from typing import Dict, Tuple
import numpy as np
from .dose import Dose
from .structure import StructureSet


def plan_quality_value(structures: StructureSet,
                       dose: Dose,
                       target_doses: Dict[str, float]) -> Tuple[Dose, np.ndarray]:
    '''
    Calculates plan quality Q according to [1].
    
    [1] https://doi.org/10.3109/0284186X.2014.906747

    Parameters
    ----------
    structures: StructureSet
        The strcuture set to be used.
    dose: Dose
        The dose for which we want to calculate the Q distribution.
    target_doses: Dict[str, float]
        Dictionary containing the dose to each target structure. The key
        is the target name and the value is its prescribed dose.

    Returns
    -------
    Q_dsitribution: Dose
        The Q-value at each voxel. The Q-value for non-target voxels is zero.
    total_target_mask: np.ndarray
        The union of the binary masks of all targets (logical or).
    '''
    if len(target_doses) == 0:
        raise RuntimeError('Need at least one target to calculate Q_RMS!')
    elif len(target_doses) == 1:
        target_dose = list(target_doses.values())[0]
        target_name = list(target_doses.keys())[0]
        Q = dose / (target_dose + 1e-6)
        return Q, structures[target_name].mask
    else:
        prescription_matrices = []
        total_target = []
        for target_name, target_dose in target_doses.items():
            mask = structures[target_name].mask
            presc = mask.astype(dose.data.dtype) * target_dose
            prescription_matrices.append(presc)
            total_target.append(mask)
        prescription_matrices = np.stack(prescription_matrices, -1)
        prescription_tensor = np.max(prescription_matrices, axis=-1)
        total_target = np.logical_or(*total_target)
        
        Q = dose.data / (prescription_tensor + 1e-6)
        # Avoid blowups due to dividing by "zero".
        # Typical Q values are around 1, so this is by far big enough.
        # The blowups are in regions receiving no dose, anways.
        # Therefore, they will not affect the calculation of Q_RMS,
        # but they will cause an out of memory error for the DVH calculation.
        Q = Q * total_target

        return Dose(Q, dose.spacing.copy(), dose.origin.copy()), total_target

def q_rms(Q: Dose, full_target_mask: np.ndarray) -> float:
    return np.sqrt(np.sum((Q.data - 1)**2 * full_target_mask) / np.sum(full_target_mask))
