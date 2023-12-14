import os
import time
from typing import Dict, List, Tuple
import numpy as np
from .configuration_factory import ConfigurationFactory
from .data import FionaPatient
from .dose import Dose
from .runner import PyFTPPRunner
from .treatment_plan import TreatmentPlan


def calculate_total_dose(target_mask: np.ndarray,
                         target_dose: float,
                         weights: List[float],
                         doses: List[Dose]) -> Dose:
    dose_shape = doses[0].data.shape
    dose_spacing = doses[0].spacing
    dose_origin = doses[0].origin

    total_dose = Dose(
        np.zeros(dose_shape, dtype=np.float32),
        dose_spacing,
        dose_origin,
    )

    for weight, dose in zip(weights, doses):
        target_mean_dose = np.sum(target_mask * dose.data) / np.sum(target_mask)
        total_dose = total_dose + (dose / target_mean_dose * target_dose * weight)

    return total_dose


class PlanEvaluator:
    
    def __init__(self,
                 patient: FionaPatient,
                 optimisation_target_name: str,
                 target_name: str, # used for normalisation
                 output_dir: str,
                 plan: TreatmentPlan,
                 constraint_builder,
                 dose_grid_resolution=0.35,
                 spot_spacing=0.4,
                 n_iterations=300):
        self._patient = patient
        self._optim_target_name = optimisation_target_name
        self._target_name = target_name
        self._output_dir = output_dir
        self._plan = plan
        self._constraint_builder = constraint_builder
        
        self._dose_grid_resolution = dose_grid_resolution
        self._spot_spacing = spot_spacing
        self._n_iterations = n_iterations

    def _build_runner(self):
        '''
        Builds a runner for evaluating a single treatment plan.
        
        The field arrangement is taken from the plan, but the constraints
        come from the constraint builder. This is a callable object with
        the optimisation ID as its only parameter.
        It should return an iterable object of constraints.
        The constraints are represented as tuples of the form
        (structure_name, dose, volume, type, importance).
        '''
        runner = PyFTPPRunner(self._output_dir)

        target_dose = self._plan.target_doses[self._plan.main_target]
        target_dose_at_edge = self._plan.target_doses_at_edge[self._plan.main_target]

        optimisation_weights = []

        for optim_ID in range(self._plan.n_optimisations):
            factory = ConfigurationFactory(self._patient,
                                           self._patient.ct_path,
                                           self._optim_target_name,
                                           target_dose,
                                           target_dose_at_edge,
                                           self._dose_grid_resolution)

            # Make sure that the fields used in the same spot optimisation have
            # the same weights.
            field_IDs = self._plan.fields_in_optimisation(optim_ID)
            weights = [self._plan.field_weight(f) for f in field_IDs]
            assert all([w == weights[0] for w in weights])
            optimisation_weights.append(len(field_IDs) * weights[0])

            # Add the fields.
            for field_ID in self._plan.fields_in_optimisation(optim_ID):
                factory.add_field(target_ID=0,
                                  gantry_angle=self._plan.gantry_angle(field_ID),
                                  couch_angle=self._plan.couch_angle(field_ID),
                                  nozzle_extraction=self._plan.nozzle_extraction(field_ID),
                                  preabsorber_setting=self._plan.preabsorber_setting(field_ID),
                                  spot_spacing=self._spot_spacing)

            # Add the constraints.
            for constr in self._constraint_builder(optim_ID):
                name, dose, volume, type, weight = constr
                factory.add_constraint(
                    structure_label=name,
                    dose=dose,
                    volume=volume,
                    type=type,
                    importance=weight,
                )

            config = factory.configurator
            config.set_n_iterations(self._n_iterations)
            runner.add_configuration(config, f'trial_{self._plan.trial_ID(optim_ID)}')

        return runner, optimisation_weights

    def evaluate(self, debug: bool = False):
        result_dose_path = f'{self._output_dir}/result_dose.dat'

        runner, optimisation_weights = self._build_runner()

        if debug:
            runner.run_all(debug=True)

            dose = Dose(
                np.zeros(self._patient.ct.data.shape),
                self._patient.ct.spacing,
                self._patient.ct.origin,
            )
            dose.save(result_dose_path)
            return result_dose_path

        #############################spot_placement_margin=0.5,
        # Run the configurations.
        #############################
        start = time.time()
        runner.run_all()
        end = time.time()
        duration = end - start
        print(f'The baseline evaluation for {self._patient.patient_ID} took {duration:.2f}s.')

        ###############################################################
        # Sum up the dose contributions from each spot optimisation.
        ###############################################################
        
        
        if self._plan.n_optimisations > 1:
            doses = []
            for optim_ID in range(self._plan.n_optimisations):
                dose = Dose.load(f'{self._output_dir}/trial_{self._plan.trial_ID(optim_ID)}/result_dose.dat')
                doses.append(dose)

            # Make sure the same grid is used for all calculations.
            dose_grid = doses[0].grid
            for dose in doses:
                assert dose_grid == dose.grid

            structure = self._patient.structures[self._target_name]
            target_mask = structure.mask

            total_dose = calculate_total_dose(
                target_mask, self._plan.target_dose, optimisation_weights, doses
            )
            total_dose.save(result_dose_path)
        else:
            link_target = result_dose_path
            if os.path.exists(link_target):
                os.remove(link_target)
            # Just create a symbolic link to have the same layout.
            os.symlink(
                f'trial_{self._plan.trial_ID(0)}/result_dose.dat',
                link_target,
            )
        
        return result_dose_path


class ImptEvaluator:
    
    def __init__(self,
                 patient: FionaPatient,
                 target_doses: Dict[str, float], # dose to each target
                 output_dir: str,
                 plan: TreatmentPlan,
                 constraint_builder,
                 dose_grid_resolution=0.35,
                 spot_spacing=0.4,
                 spot_placement_margin=0.5,
                 n_iterations=300,
                 additional_target_importance: float = 2.,
                 hotspot_avoidance_importance = None,
                 # 1.1 means that doses >= 110% are considered hotspots.
                 hotspot_avoidance_threshold_fraction = None,
                 # Used for spot placement.
                 optimisation_target = None,
                 optimisation_target_dose = None):
        self._patient = patient
        self._target_doses = target_doses
        self._output_dir = output_dir
        self._plan = plan
        self._constraint_builder = constraint_builder
        
        self._dose_grid_resolution = dose_grid_resolution
        self._spot_spacing = spot_spacing
        self._spot_placement_margin = spot_placement_margin
        self._n_iterations = n_iterations
        self._additional_target_importance = additional_target_importance
        self._hotspot_avoidance_importance = hotspot_avoidance_importance
        self._hotspot_avoidance_threshold_fraction = hotspot_avoidance_threshold_fraction
        self._optimisation_target = optimisation_target

    def _minimum_dose_target(self, target_doses) -> str:
        names = []
        doses = []
        for name, dose in target_doses.items():
            names.append(name)
            doses.append(dose)
        i = np.argmin(doses)
        return names[i]

    def _build_runner(self):
        '''
        Builds a runner for evaluating a single treatment plan.
        
        The field arrangement is taken from the plan, but the constraints
        come from the constraint builder. This is a callable object taking
        no parameters.
        It should return an iterable object of constraints.
        The constraints are represented as tuples of the form
        (structure_name, dose, volume, type, importance).
        '''
        runner = PyFTPPRunner(self._output_dir)

        if self._optimisation_target is None:
            optimisation_target = self._minimum_dose_target(self._plan.target_doses)
            target_dose = self._plan.target_doses[optimisation_target]
            target_dose_at_edge = self._plan.target_doses_at_edge[optimisation_target]
        else:
            optimisation_target = self._optimisation_target
            minimum_dose_target = self._minimum_dose_target(self._plan.target_doses)
            target_dose = self._plan.target_doses[minimum_dose_target]
            target_dose_at_edge = self._plan.target_doses_at_edge[minimum_dose_target]

        factory = ConfigurationFactory(self._patient,
                                       self._patient.ct_path,
                                       optimisation_target,
                                       target_dose,
                                       target_dose_at_edge,
                                       self._dose_grid_resolution)

        # Add the fields.
        for field_ID in self._plan.field_IDs:
            factory.add_field(target_ID=0,
                              gantry_angle=self._plan.gantry_angle(field_ID),
                              couch_angle=self._plan.couch_angle(field_ID),
                              nozzle_extraction=self._plan.nozzle_extraction(field_ID),
                              preabsorber_setting=self._plan.preabsorber_setting(field_ID),
                              spot_spacing=self._spot_spacing,
                              margin=self._spot_placement_margin)

        # Add the constraints for the OARs.
        for constr in self._constraint_builder():
            name, dose, volume, type, weight = constr
            factory.add_constraint(
                structure_label=name,
                dose=dose,
                volume=volume,
                type=type,
                importance=weight,
            )
        
        # Add objectives for the targets other than the main optimisation
        # target. The latter is alreadz accounted for by the optimiser, so no
        # objective is needed.
        maximum_target_dose = self._patient.prescriptions.hottest_target()[1]

        for target_name, dose in self._target_doses.items():
            if target_name == optimisation_target:
                continue

            factory.add_constraint(
                structure_label=target_name,
                dose=dose / target_dose,
                volume=1.0,
                type='DOSE_VOLUME',
                importance=self._additional_target_importance,
                limitType='LOWER',
            )
            factory.add_constraint(
                structure_label=target_name,
                dose=maximum_target_dose / target_dose,
                volume=0.0,
                type='DOSE_VOLUME',
                importance=self._additional_target_importance,
                limitType='UPPER',
            )
        
        # Punish hot spots.
        if self._hotspot_avoidance_importance is not None:
            for target_name, dose in self._plan.target_doses.items():
                factory.add_constraint(
                    structure_label=target_name,
                    # Punish hot spots >= 100%*self._hotspot_avoidance_threshold_fraction.
                    dose=dose / target_dose * self._hotspot_avoidance_threshold_fraction,
                    volume=0.01,                   # Approximate Dmax with D1%
                    type='DOSE_VOLUME',
                    importance=self._hotspot_avoidance_importance,
                    limitType='UPPER',
                )

        config = factory.configurator
        config.set_n_iterations(self._n_iterations)
        runner.add_configuration(config, f'optim_0')

        return runner

    def evaluate(self, debug: bool = False):
        result_dose_path = f'{self._output_dir}/result_dose.dat'

        runner = self._build_runner()

        if debug:
            runner.run_all(debug=True)

            dose = Dose(
                np.ones(self._patient.ct.data.shape),
                self._patient.ct.spacing,
                self._patient.ct.origin,
            )
            dose.save(result_dose_path)
            return result_dose_path

        #############################
        # Run the configurations.
        #############################
        start = time.time()
        runner.run_all()
        end = time.time()
        duration = end - start
        print(f'The baseline evaluation for {self._patient.patient_ID} took {duration:.2f}s.')

        # Just create a symbolic link to have the same result dose path as for SFUD/SFO plans.
        link_target = result_dose_path
        if os.path.exists(link_target):
            os.remove(link_target)

        os.symlink(
            'optim_0/result_dose.dat',
            link_target,
        )

        return result_dose_path
