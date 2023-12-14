from typing import Dict, Tuple
import numpy as np
from pyftpp import ConfigurationFactory, FionaPatient, PyFTPPRunner
from pyftpp.fiona_dij_handler import FionaDijHandler, mask_at_indices
from pyftpp.optimise.util import mean_dose
# from pyftpp.optimise.util import calculate_dose, mean_dose, maximum_dose, dvh_d, dvh_v
# from pyftpp.optimise.objectives import l1_objective, l2_objective, leq_objective, geq_objective


class TensorflowSpotOptimisation:

    def __init__(self,
                 patient: FionaPatient,
                 spot_placement_target_name: str,
                 normalisation_structure_name: str,
                 normalisation_dose: str,
                 ct_path: str,
                 output_dir: str):
        self._patient = patient
        self._spot_placement_target_name = spot_placement_target_name
        self._normalisation_structure_name = normalisation_structure_name
        self._normalisation_dose = normalisation_dose
        self._ct_path = ct_path
        self._output_dir = output_dir

        # D[i, j] = dose contributed by spot j to voxel i if w[j] = 1.
        self._optimisation_point_indices = None
        self._Dij = None
        self._w = None

    def build_Dij(self,
                  angles,
                  nozzle_extraction,
                  spot_spacing,
                  dose_grid_resolution=0.35,
                  optimisation_grid_resolution=0.35,
                  debugging=False):
        factory = ConfigurationFactory(
            self._patient,
            self._ct_path,
            self._spot_placement_target_name,
            self._normalisation_dose,
            0.9 * self._normalisation_dose,
            dose_grid_resolution=dose_grid_resolution,
            optimisation_grid_resolution=optimisation_grid_resolution,
        )

        for (gantry_angle, couch_angle) in angles:
            factory.add_field(
                gantry_angle,
                couch_angle,
                nozzle_extraction,
                spot_spacing,
                preabsorber_setting='auto',
                margin=0.5,
                target_ID=0,
            )

        config = factory.configurator

        if not debugging:
            runner = PyFTPPRunner(self._output_dir, prune_spots=False)
            runner.add_configuration(config, self._patient.patient_ID)
            runner.run_all(log_dij_matrix=True)

    def load_Dij(self):
        loader = FionaDijHandler(f'{self._output_dir}/{self._patient.patient_ID}', self._patient.ct.grid)
        Dij = loader.Dij.tocoo()
        spot_weights = loader.spots['weight'].values
        optimisation_point_indices = loader.optimisation_point_grid_indices
        # It can happen that there are some empty lines at the end of the
        # Dij matrix. In that case the dose to the corresponding optimisation
        # points is always zero and they can be removed during the optimisation.
        optimisation_point_indices = optimisation_point_indices.iloc[:Dij.shape[0], :]

        indices = np.stack([Dij.row, Dij.col], axis=1)

        self._loader = loader
        self._optimisation_point_indices = optimisation_point_indices
        self._Dij = Dij
        self._w = spot_weights.astype('float32').reshape((-1, 1))

    def structure_mask(self, structure_name: str):
        return mask_at_indices(
            self._optimisation_point_indices.loc[:, ['x', 'y', 'z']].values,
            self._patient.structures[structure_name].mask,
        ).reshape((-1, 1)).astype('float32')

    @property
    def normalisation_structure_name(self):
        return self._normalisation_structure_name

    @property
    def normalisation_dose(self):
        return self._normalisation_dose

    @property
    def Dij(self):
        return self._Dij

    @property
    def w(self):
        return self._w

    @w.setter
    def w(self, w):
        self._w = w
        self._loader.set_spot_weights(np.clip(w, a_min=0, a_max=None))

    @property
    def dose(self):
        return calculate_dose(self.Dij, self.w)

    def write_spots(self):
        self._loader.write_spots()

    def initialise_optimisation(self, normalisation_dose=None):
        if normalisation_dose is None:
            normalisation_dose = self.normalisation_dose
        normalisation_mask = self.structure_mask(self.normalisation_structure_name)
        
        w = np.ones((self.Dij.shape[1], 1), dtype='float32')
        dose = self.Dij.dot(w)

        # Normalise s. t. the Dmean of the target is correct.
        mean = mean_dose(dose, normalisation_mask)
        w *= normalisation_dose / mean
        self._w = w

    @property
    def optimiser(self):
        return self.optimiser

    @optimiser.setter
    def optimiser(self, optimiser):
        self.optimiser = optimiser
