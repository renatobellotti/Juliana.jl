import logging
from typing import List, Dict, Tuple, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from .ct import CT
from .data import DataSet
from .dose import Dose
from .dvh import DVH, cc_to_volume_fraction
from .plotting import DvhPlotter
from .prescriptions import Constraint, ConstraintType
from .structure import load_array_structure
from .treatment_plan import TreatmentPlan


class DoseComparison:
    '''
    Class that loads one or more dose distributions and compares them.
    '''

    def __init__(self,
                 dataset: DataSet,
                 doses: Dict[str, Dose],
                 plans: Union[None, List[TreatmentPlan]] = None,
                 use_near_max: bool = False):
        '''
        Parameters
        ----------
        dataset: DataSet
            Patient data.
        doses: Dict[str, Dose]
            Dose distributions to be compared. The keys of the dict are the
            labels of the distributions.
        plans: None or list of TreatmentPlan
            Optional. The treatment plans that lead to the dose distributions.
            The objectives for the spot optimiser will be added to the DVH
            plots.
        use_near_max: bool
            Optional. Whether to replace Dmax constraints with a more stable
            approximation. The D0.035cc is used as suggested by literature
            (https://doi.org/10.1007/s00066-018-1416-x).
            We use it regardless of the structure volume in order to
            be consistent.
        '''
        self._patient = dataset
        self._plans = plans
        self._use_near_max = use_near_max
        self._ct = dataset.ct

        # Load the dose distributions.
        self._labels = []
        self._doses = []
        for label, dose in doses.items():
            self._labels.append(label)
            self._doses.append(dose)

        # Sanity check.
        self._dose_shape = self.doses[0].data.shape
        self._dose_spacing = self.doses[0].spacing
        self._dose_origin = self.doses[0].origin

        same_shape = True
        same_spacing = True
        same_origin = True
        for dose in self._doses:
            same_shape = same_shape and (self._dose_shape == dose.data.shape)
            same_spacing = same_spacing and (self._dose_spacing == dose.spacing).all()
            same_origin = same_origin and (self._dose_origin == dose.origin).all()
        
        if not same_shape:
            logging.warning('Comparing dose distributions of different shapes.')
        if not same_spacing:
            logging.warning('Comparing dose distributions of different spacings.')
        if not same_origin:
            logging.warning('Comparing dose distributions of different origins.')

        # Load the prescriptions.
        self._prescriptions = dataset.prescriptions
        self._target_dose = self._get_main_target_dose()

        # Build a report table summarising the prescription constraints and
        # their achieved value for each dose distribution.
        self._report = self._build_prescription_report()

    def _get_main_target(self) -> str:
        target = self._plans[0].main_target

        if len(self._plans) == 1:
            return target

        for plan in self._plans[1:]:
            assert plan.main_target == target

        return target
    
    @property
    def target_doses(self) -> Dict[str, float]:
        target_doses = {}
        for plan in self._plans:
            target_doses.update(plan.target_doses)
        return target_doses

    def _get_main_target_dose(self) -> float:
        target = self._get_main_target()
        target_dose = self._plans[0].target_doses[target]

        if len(self._plans) == 1:
            return target_dose
        
        for plan in self._plans[1:]:
            assert np.isclose(
                plan.target_doses[target],
                target_dose,
                atol=1e-8,
            )
        
        return target_dose

    @property
    def doses(self) -> List[Dose]:
        return self._doses

    @property
    def constraint_report(self):
        '''
        Returns a pandas DataFrame that compares how well all doses fulfill
        the prescriptions.
        '''
        return self._report

    def _build_prescription_report(self):
        entries = []

        # Obtain the achieved values for each constraint, for every dose.
        for constr in self._prescriptions.constraints:
            if constr.kind == ConstraintType.UNKNOWN:
                continue

            if (constr.kind == ConstraintType.MAX) and self._use_near_max:
                volume_fraction = cc_to_volume_fraction(
                    0.035,
                    self._patient.structures[constr.structure_name].mask,
                    self._ct.spacing,
                )
                constr = Constraint(
                    constr.structure_name,
                    ConstraintType.DVH_D,
                    constr.dose,
                    constr.priority,
                    volume_fraction,
                    constr.direction,
                )

            entry = {
                'structure_name': constr.structure_name,
                'dose': constr.dose,
                'volume': constr.volume,
                'kind': constr.kind.name,
                'priority': constr.priority,
            }

            for label, dose in zip(self._labels, self._doses):
                achieved_dose = _evaluate_constr(
                    dose.data,
                    self._patient.structures[constr.structure_name].mask,
                    constr
                )
                entry[f'{label}_dose'] = achieved_dose
            entries.append(entry)

        df = pd.DataFrame(entries)
        if df.shape[0] == 0:
            columns = ['structure_name', 'dose', 'volume', 'kind']
            for label in self._labels:
                columns.append(f'{label}_dose')
            df = pd.DataFrame(columns=columns)

        return df

    def plot_dvh(self, structure_name: str, plot_optimiser_objectives: bool = False) -> plt.Figure:
        plotter = DvhPlotter(self._target_dose, structure_name)

        # Add DVH curves.
        mask = self._patient.structures[structure_name].mask
        for i, (label, dose) in enumerate(zip(self._labels, self.doses)):
            dvh = DVH(dose.data, mask, delta_d=0.005)
            plotter.add_curve(label, dvh, sns.color_palette()[i])

        # Add markers for the spot optimiser objectives.
        if (self._plans is not None) and plot_optimiser_objectives:
            for i, (label, plan) in enumerate(zip(self._labels, self._plans)):
                for optim_ID in range(plan.n_optimisations):
                    df = plan.constraints(optim_ID)

                    for _, row in df.query('`oar` == @structure_name').iterrows():
                        if plan.n_optimisations == 1:
                            text = f'weight={row["weight"]}'
                        elif plan.n_optimisations > 1:
                            text = f'weight={row["weight"]} (optim {optim_ID})'.strip()
                        else:
                            raise RuntimeError(
                                'Unexpected number of optimisations'
                            )

                        plotter.add_marker(
                            row['dose'] * 100,
                            row['volume'] * 100,
                            sns.color_palette()[i],
                            text=text
                        )

        # Add a marker for prescription constraints.
        structure_constraints = self.constraint_report.query(
            '`structure_name` == @structure_name'
        ).sort_values('structure_name')

        for _, constr in structure_constraints.iterrows():
            kind = constr['kind']

            dose = None
            volume = None
            color = ''

            if kind == 'MAX':
                dose = constr['dose'] / self._target_dose * 100.
                volume = 0
            elif kind == 'MEAN':
                # Noting to do, mean constraints cannot be displayed in DVHs.
                pass
            elif (kind == 'DVH_D') or (kind == 'DVH_V'):
                volume = constr['volume'] * 100.
                dose = constr['dose'] / self._target_dose * 100.
            else:
                logging.error(f'Unknown constraint: {constr}')

            if dose is not None:
                plotter.add_marker(dose, volume, 'black', text='Prescription')
        
        if structure_name in self.target_doses:
            plotter.add_vertical_line(
                self.target_doses[structure_name],
                color='red'
            )

        return plotter.plot()


# Helper functions.
def _load_dose(path: str) -> Dose:
    '''
    Parameters
    ----------
    path: str
        Path to a .dat file containing a dose distribution.
    '''
    path = Path(path).resolve()
    if path.is_file() and str(path).endswith('.dat'):
        return Dose.load(str(path))
    else:
        raise RuntimeError(f'Expected the path to be a .dat file, but got: {str(path)}')


def _evaluate_constr(dose_matrix, mask, constraint: Constraint, delta_d=0.005):
    '''
    If kind in ["MEAN", "MAX", "DVH_D"]:
        returns dose to the constraint quantity in Gy.
    Else:
        returns volume as a fraction in [0, 1].
    '''
    if constraint.kind == ConstraintType.MEAN:
        return np.sum(dose_matrix * mask) / np.sum(mask)
    elif constraint.kind == ConstraintType.MAX:
        return np.max(dose_matrix * mask)
    elif constraint.kind == ConstraintType.DVH_D:
        dvh = DVH(dose_matrix, mask, delta_d)
        return dvh.D(constraint.volume)
    elif constraint.kind == ConstraintType.DVH_V:
        dvh = DVH(dose_matrix, mask, delta_d)
        return dvh.V(constraint.dose)
    else:
        raise RuntimeError(f'Unsupported constraint kind: {kind}')
