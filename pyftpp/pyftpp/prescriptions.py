from dataclasses import dataclass
from enum import Enum
import logging
import json
import re
from typing import Dict, List, Tuple, Union
import numpy as np
import pandas as pd
from .structure import StructureSet

class ConstraintType(Enum):
    DVH_D = 'DOSE_VOLUME'
    DVH_V = 'DOSE_VOLUME'
    MEAN = 'D_mean'
    MAX = 'D_max'
    UNKNOWN = 'unknown'


class ThresholdType(Enum):
    UPPER = 'UPPER'
    LOWER = 'LOWER'


class ConstraintPriority(Enum):
    SOFT = 'soft'
    HARD = 'hard'
    UNKNOWN = ''


@dataclass
class Constraint:
    structure_name: str
    kind: ConstraintType
    dose: float
    priority: ConstraintPriority
    volume: Union[None, float] = None
    direction: ThresholdType = ThresholdType.UPPER

    def __post_init__(self):
        if (self.kind in [ConstraintType.DVH_D, ConstraintType.DVH_V]) \
            and (self.volume is None):
            raise RuntimeError(
                'Constraint type is DVH, but no volume is specified.'
            )


@dataclass
class Prescriptions:
    # Dose to each target structure by name. Unit: Gy.
    target_doses: Dict[str, float]
    constraints: List[Constraint]

    def hottest_target(self) -> Tuple[str, float]:
        '''
        Get the name and the dose of the target with the highest dose
        prescription.

        Returns
        -------
        target_name: str
        target_dose: float
        '''
        target_name = None
        maximum = None

        for name, target_dose in self.target_doses.items():
            if (maximum is None) or (target_dose > maximum):
                target_name = name
                maximum = target_dose

        return target_name, maximum

    def coldest_target(self) -> Tuple[str, float]:
        '''
        Get the name and the dose of the target with the lowest dose
        prescription.

        Returns
        -------
        target_name: str
        target_dose: float
        '''
        target_name = None
        minimum = None

        for name, target_dose in self.target_doses.items():
            if (minimum is None) or (target_dose < minimum):
                target_name = name
                minimum = target_dose

        return target_name, minimum


def load_fiona_prescriptions(patient_ID: str,
                             target_doses_path: str,
                             oar_presc_path: str,
                             grid_spacing: np.ndarray,
                             masks: Dict[str, np.ndarray]) -> Prescriptions:
    target_doses = _load_target_doses(target_doses_path, patient_ID)

    # The 100% dose level.
    target_dose = max(target_doses.values())

    constraints = parse_fiona_constraint_file(
        oar_presc_path,
        target_dose,
        grid_spacing,
        masks
    )

    return Prescriptions(
        target_doses,
        constraints
    )


def _load_target_doses(path: str, patient_ID: str) -> Dict[str, float]:
    '''
    Loads the dose to each target from a CSV file at the given path.

    Parameters
    ----------
    path: str
        Path to a CSV file containing the dose to each target.
        The CSV file must contain the following columns:
        patient_ID, target_name, total_dose_gy, dose_gy_per_frac
    
    Returns
    -------
    Dict[str, float]:
        Dose in Gy for each target. The target name is the dictionary key.
    '''
    df = pd.read_csv(path)
    df.query('`patient_ID` == @patient_ID', inplace=True)

    doses = {}
    for _, row in df.iterrows():
        doses[row['target_name']] = row['total_dose_gy']
    return doses


def parse_fiona_constraint_file(path: str,
                                target_dose: float,
                                grid_spacing: np.ndarray,
                                structures: StructureSet) -> List[Constraint]:
    '''
    The parameters are needed to parse relative constraints.
    '''
    df = pd.read_csv(path)
    constraints = []
    for i, row in df.iterrows():
        structure_name = row['structure_name']
        if structure_name not in structures.structure_names:
            logging.error(f'Could not find mask for {structure_name}')
            continue
        mask = structures[structure_name].mask
        structure_volume = np.sum(mask) * np.prod(grid_spacing)

        priority = row['soft_hard']
        if priority == 'hard':
            priority = ConstraintPriority.HARD
        elif priority == 'soft':
            priority = ConstraintPriority.SOFT
        elif (priority == '') or (priority == 'unknown') or (np.isnan(priority)):
            priority == ConstraintPriority.UNKNOWN
        else:
            raise RuntimeError(f'Unexpected constraint priority: {priority}')

        kind, dose, volume = _parse_constraint(
            structure_name,
            row['constraint_quantity'].strip(),
            row['constraint_threshold'].strip(),
            target_dose,
            structure_volume,
        )
        constraint = Constraint(structure_name, kind, dose, priority, volume)
        constraints.append(constraint)
    return constraints


def _parse_constraint(structure_name: str,
                      quantity: str,
                      threshold: str,
                      target_dose: float,
                      structure_volume: float) -> Constraint:
    '''
    Parse input of the form:
        quantity < threshold
    
    The target dose is used to parse relative constraints like D2 <= 103%.

    Returns
    -------
    kind: ConstraintType
        Kind of constraint.
    dose: float
        For "MEAN", "MAX", "DVH_D" constraints: Dose threshold not to be
        exceeded. Constraint dose ("V<dose>") for "DVH_V". Unit: Gy
    volume: None or float
        For "MEAN", "MAX": None.
        For "DVH_D": Constraint volume ("D<volume>"), as a percentage in [0, 100].
        For "DVH_V": Volume threshold not to be exceeded, as a percentage in [0, 100].
    '''
    dvh_D_pattern = re.compile(r'D\d+(\.\d+)?(cc)?(%)?\Z', re.IGNORECASE) # Ex.: D2, D4cc
    dvh_V_pattern = re.compile(r'V\d+(\.\d+)?(gy|%)?\Z', re.IGNORECASE) # Ex.: V60gy, V30

    # Get kind, dose and volume of the constraint.
    # If dose and/or volume are not set: Return -1 for the unset quantities.
    # MEAN and MAX constraints don't have an associated volume values.
    dose = -1
    volume = None
    if structure_name == 'unknown':
        kind = ConstraintType.UNKNOWN
        dose = -1
        volume = -1
    if quantity == 'D_mean':
        kind = ConstraintType.MEAN
        dose = _parse_dose(threshold, target_dose)
    elif quantity == 'D_max':
        kind = ConstraintType.MAX
        dose = _parse_dose(threshold, target_dose)
    elif dvh_D_pattern.match(quantity):
        # D<volume> < threshold
        kind = ConstraintType.DVH_D
        # Remove the "D" at the beginning.
        volume = _parse_volume(quantity[1:], structure_volume)
        dose = _parse_dose(threshold, target_dose)
    elif dvh_V_pattern.match(quantity):
        # V<dose> < threshold
        kind = ConstraintType.DVH_V
        # Remove the "V" at the beginning.
        dose = _parse_dose(quantity[1:], target_dose)
        volume = _parse_volume(threshold, structure_volume)
    elif quantity == 'unknown':
        logging.warning(
            f'Skipping unknown constraint: ' \
            f'{quantity} <= {threshold} for {structure_name}'
        )
        kind = ConstraintType.UNKNOWN
        dose = -1
        volume = -1
    else:
        raise RuntimeError(f'Unsupported constraint quantity: {quantity}')

    return kind, dose, volume


def is_percentage(query: str):
    percentage_pattern = re.compile(r'\d+.?\d*%')
    return percentage_pattern.match(query)


def _parse_dose(dose: str, target_dose_gy: float) -> float:
    '''
    Convert a dose string to a float representing the dose in Gy.
    '''
    if dose.lower().endswith('gy'):
        dose = float(dose[:-2])
    elif is_percentage(dose):
        dose = float(dose[:-1]) * target_dose_gy / 100.
    else:
        raise RuntimeError(f'Unknown dose format: {dose}')
    return dose


def _parse_volume(volume: str, structure_volume: float) -> float:
    '''
    Convert volume string to a float representing the volume relative to the
    target volume.

    Possible formats (case insensitive):
    30:   30% of the structure volume; return 0.3
    10cc: 10cm^3;                      return 10 / structure_volume

    Parameters
    ----------
    volume: str
        Volume string to be converted to a float in [0, 1]. Ex.: V60Gy.
        I. e. the volume of the structure that receives 60Gy or less.
    structure_volume: float
        Volume of the structure in cm^3 (== cc).

    Returns
    -------
    float
        Volume as a fraction of the total structure volume in [0, 1].
    '''
    number_only_pattern = re.compile(r'\d+\Z')
    if is_percentage(volume):
        return float(volume[:-1]) / 100.
    elif number_only_pattern.match(volume):
        return float(volume) / 100.
    elif volume.lower().endswith('cc'):
        return float(volume[:-2]) / structure_volume
    else:
        raise RuntimeError(f'Malformatted volume string: {volume}')
