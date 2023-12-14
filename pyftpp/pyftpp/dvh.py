import numpy as np
from numba import jit
from .dose import Dose


@jit(nopython=True)
def _build_dvh(dose_tensor: np.ndarray,
               d_max: float,
               mask: np.ndarray,
               delta_d: float):
    # Make sure there is a bin also for the maximum dose value.
    bins = np.arange(0., d_max+delta_d, delta_d)
    cnt = np.zeros_like(bins)

    for i in range(dose_tensor.shape[0]):
        for j in range(dose_tensor.shape[1]):
            for k in range(dose_tensor.shape[2]):
                # Only add the dose to the histogram if it is part of the
                # structure of interest.
                if mask[i, j, k]:
                    ind = int(dose_tensor[i, j, k] // delta_d)
                    cnt[ind] += 1

    # So far we just have a histogram. Now we want to make it cumulative,
    # i. e. the bin for the lowest dose should also contain the counts of
    # the higher dose bins.
    cnt = np.cumsum(cnt[::-1])[::-1]

    return bins, cnt


class DVH:

    def __init__(self, dose: Dose, mask: np.ndarray, delta_d: float = 0.005):
        '''
        Spacing in Gy.
        '''
        self._dose_tensor = dose.data
        self._mask = mask
        self._delta_d = delta_d

        # These two arrays contain the DVH curve at discrete points.
        # Between the points, the curve is interpolated linearly.
        n_voxels = np.sum(mask)

        bins, cnt = _build_dvh(self._dose_tensor,
                               np.max(self._dose_tensor),
                               self._mask,
                               self._delta_d)
        self._doses = bins
        self._volume_fractions = cnt / n_voxels

    def V(self, dose):
        '''
        Returns the volume fraction receiving at least dose_fraction [Gy] dose.

        Linear interpolation is used.

        Ex.: V(2.4) returns V2.4Gy dose.

        Parameters
        ----------
        dose: array-like
            A float or a numpy array or list of evaluation positions.
            Dose in Gy at which to evaluate the DVH.
        '''
        return np.interp(dose, self._doses, self._volume_fractions)

    def D(self, volume_fraction):
        '''
        Returns the dose of the hottest volume_fraction voxels in Gy.

        Linear interpolation is used.

        Ex.: D(0.02) returns D2 dose.

        Parameters
        ----------
        volume_fraction: array-like
            A float or a numpy array or list of evaluation positions.
        '''
        return np.interp(volume_fraction,
                         self._volume_fractions[::-1],
                         self._doses[::-1])


def cc_to_volume_fraction(volume_cc: float,
                          mask: np.ndarray,
                          grid_spacing: np.ndarray) -> float:
    '''
    Calculate things like D0.035cc for a given dose volume histogram.
    This is the "dose received by the hottest 0.035cc of the structure".

    Parameters
    ----------
    volume: float
        The volume in cc (== cm^3) at which to evaluate the DVH.
    dvh: DVH
        The dose volume histogram to evaluate.
    grid_spacing: np.ndarray
        Grid used to compute the volume of a voxel. Shape: (3,)
    '''
    voxel_volume = grid_spacing[0] * grid_spacing[1] * grid_spacing[2]

    # Number of voxels in the structure.
    n_voxels_structure = np.sum(mask)
    # Number of voxels in the given volume.
    n_voxels_volume = volume_cc / voxel_volume

    volume_fraction = n_voxels_volume / n_voxels_structure

    return volume_fraction