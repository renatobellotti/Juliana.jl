import json
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from .data import FionaPatient


class TreatmentPlan:
    '''
    Abstraction layer representing a treatment plan.
    We consider the ensemble of all decisions made by human planners as the
    treatment plan.

    This class loads the prescriptions, beam angles and constraints.
    '''

    def __init__(self, patient: FionaPatient, plan_file: str, override_preabsorber=None):
        with open(plan_file, 'r') as file:
            self._plan = json.load(file)

        self._target_doses = patient.prescriptions.target_doses.copy()
        self._target_doses_at_edge = {k: 0.9 * v for k, v in self._target_doses.items()}
        self._fields = pd.DataFrame(self._plan['fields']).set_index('field_ID')

        self._override_preabsorber = override_preabsorber

    @property
    def n_optimisations(self) -> int:
        return len(self._plan['optimisations'])

    @property
    def field_IDs(self) -> List[str]:
        return [f['field_ID'] for f in self._plan['fields']]

    def gantry_angle(self, field_ID) -> float:
        return self._fields.loc[field_ID, 'gantry_angle']

    def couch_angle(self, field_ID) -> float:
        return self._fields.loc[field_ID, 'couch_angle']

    def nozzle_extraction(self, field_ID) -> float:
        return self._fields.loc[field_ID, 'nozzle_extraction']

    def preabsorber_setting(self, field_ID) -> str:
        if self._override_preabsorber is not None:
            assert field_ID in self._fields.index
            return self._override_preabsorber
        return self._fields.loc[field_ID, 'preabsorber']

    def target_structure(self, field_ID) -> str:
        '''
        Target structure for the given field. This is where the spots are placed.
        '''
        return self._fields.loc[field_ID, 'target_structure']

    def field_weight(self, field_ID) -> float:
        return self._fields.loc[field_ID, 'weight']

    def fields_in_optimisation(self, optim_ID: int) -> List[str]:
        return self._plan['optimisations'][optim_ID]['fields']

    @property
    def main_target(self) -> str:
        '''
        Returns the name of the target that receives the most dose.
        '''
        names = []
        doses = []
        for name, dose in self.target_doses.items():
            names.append(name)
            doses.append(dose)
        i = np.argmax(doses)
        return names[i]

    @property
    def target_doses(self) -> float:
        return self._target_doses

    @property
    def target_doses_at_edge(self) -> float:
        return self._target_doses_at_edge

    def constraints(self, optim_ID: int) -> pd.DataFrame:
        '''
        Returns a pd.DataFrame containing the settings used to optimise the
        spots in PSIPlan.

        The following columns are present:
        - 'oar':    str
          To which organ-at-risk (OAR) the constraint applies.
        - 'dose':   float
          Dose threshold relative to the target dose. A value of 0.5 means that
          no more than 50% of the target dose should be received by the OAR DVH
          indicated by 'volume'.
        - 'volume': float
          Which DVH metric is constrained. The value is given relative to the
          total volume.
        - 'weight': float
          How important the constraint is relative to the target coverage.
          The target coverage has an importance of 1. Values < 1 assign the
          constraint less importance than the target coverage, values > 1 more.
          This is the value that is probably most different between PSIPlan
          and Fiona?

        Example for interpretation of constraints:
        dose = 0.8, volume = 0.2
        D20 < 80%, i. e. "at least 20% of the OAR receive no more dose than 80%
        of the target dose"

        WARNING:
        The optimisation settings, optimisation grid and beam data are
        different than in PSIPlan. Therefore, the constraints might not
        always be useful in Fiona. Be careful when using these values to
        generate new plans!
        '''
        optim = self._plan['optimisations'][optim_ID]
        objectives = _parse_objectives(optim['objectives'])
        constraints = _parse_constraints(optim['constraints'])

        objectives.query('`weight` > 0', inplace=True)
        objectives.rename(columns={'target': 'oar'}, inplace=True)
        # We constrain the V95.
        objectives['volume'] = 0.95

        constraints['limitType'] = 'UPPER'
        objectives['limitType'] = 'LOWER'

        constraints = pd.concat([constraints, objectives], ignore_index=True)

        constraints.query('`weight` > 0', inplace=True)
        constraints.sort_values(['oar', 'weight'], inplace=True)

        return constraints

    def n_iterations(self, optim_ID: int) -> int:
        optim = self._plan['optimisations'][optim_ID]
        return optim['n_iterations']

    def trial_ID(self, optim_ID: int) -> int:
        optim = self._plan['optimisations'][optim_ID]
        return optim["trial_ID"]


######################
# Helper functions.
######################
def _parse_objectives(objectives):
    if len(objectives) == 0:
        return pd.DataFrame([], columns=['target', 'weight', 'dose'])

    objectives = pd.DataFrame(objectives).astype({
        'target': str,
        'weight': int,
        'dose': int,
    })
    objectives['dose'] /= 1000
    objectives['weight'] /= 10_000
    return objectives


def _parse_constraints(constraints):
    if len(constraints) == 0:
        constraints = pd.DataFrame(columns=[
            'oar',
            'weight',
            'dose',
            'volume',
        ])
    else:
        constraints = pd.DataFrame(constraints)

    constraints = constraints.astype({
        'oar': str,
        'weight': int,
        'dose': int,
        'volume': int,
    })
    constraints[['dose', 'volume']] /= 1000
    constraints['weight'] /= 10_000
    return constraints
