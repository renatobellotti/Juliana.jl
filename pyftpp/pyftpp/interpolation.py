import numpy as np
from scipy.interpolate import interp2d
from .grid import Grid


def interpolate_data(coarse_data: np.ndarray,
                     coarse_grid: Grid,
                     fine_grid: Grid) -> np.ndarray:
    '''
    Interpolate the coarse dose distribution onto the fine grid.
    Return the resulting distribution.
    '''
    # Check preconditions:
    # - Same spacing in x- and y-direction.
    # - Same slice spacing.
    # - The finer grid is fully contained within the coarser grid.
    same_res_coarse = np.isclose(
        coarse_grid.spacing[0],
        coarse_grid.spacing[1],
        atol=1e-6
    )
    same_res_fine = np.isclose(
        fine_grid.spacing[0],
        fine_grid.spacing[1],
        atol=1e-6
    )
    if not (same_res_coarse and same_res_fine):
        raise NotImplementedError(
            'Cannot interpolate grids of different spacings ' \
            'when the grid spacing is different in the x and y ' \
            'direction.'
        )
    if not np.isclose(coarse_grid.spacing[2], fine_grid.spacing[2], atol=1e-6):
        raise NotImplementedError(
            'Cannot interpolate grids of different spacings when the ' \
            'slice spacing is not the same.'
        )

    is_contained = (fine_grid.bounding_box.x_min >= coarse_grid.bounding_box.x_min) \
        and (fine_grid.bounding_box.y_min >= coarse_grid.bounding_box.y_min) \
        and (fine_grid.bounding_box.z_min >= coarse_grid.bounding_box.z_min) \
        and (fine_grid.bounding_box.x_max <= coarse_grid.bounding_box.x_max) \
        and (fine_grid.bounding_box.y_max <= coarse_grid.bounding_box.y_max) \
        and (fine_grid.bounding_box.z_max <= coarse_grid.bounding_box.z_max)

    if not is_contained:
        raise NotImplementedError(
            'Cannot add dose distributions if the finer grid is ' \
            'not contained in the coarser grid.'
        )
    
    # Not sure if the implementation is correct if the origins are
    # not the same...
    if not np.isclose(fine_grid.origin, coarse_grid.origin, atol=1e-6).all():
        raise NotImplementedError(
            'Can only add Dose objects if the origins are the same.'
        )
    
    # Interpolate the values of the coarse grid on the fine grid,
    # slice by slice.
    interpolated = np.empty(fine_grid.shape, dtype=coarse_data.dtype)
    for slice_ind in range(fine_grid.shape[2]):
        m = _interpolate_slice(
            coarse_data, coarse_grid, fine_grid, slice_ind
        )
        interpolated[:, :, slice_ind] = m
    
    return interpolated

    
def _interpolate_slice(coarse_data: np.ndarray,
                       coarse_grid: Grid,
                       fine_grid: Grid,
                       slice_ind: int) -> np.ndarray:
    # The voxel positions are in the center of the voxel, so we have to shift
    # the interpolation points by 0.5 times the spacing.
    coarse_spacing = np.arange(
        0,
        coarse_grid.shape[0] * coarse_grid.spacing[0], 
        coarse_grid.spacing[0],
    ) + 0.5 * coarse_grid.spacing[0]
    fine_spacing = np.arange(
        0,
        fine_grid.shape[0] * fine_grid.spacing[0], 
        fine_grid.spacing[0],
    ) + 0.5 * fine_grid.spacing[0]

    f = interp2d(
        coarse_spacing + coarse_grid.origin[0],
        coarse_spacing + coarse_grid.origin[1],
        coarse_data[:, :, slice_ind]
    )

    interpolated = f(
        fine_spacing + fine_grid.origin[0],
        fine_spacing + fine_grid.origin[1],
    )
    interpolated = interpolated.reshape(fine_grid.shape[:2])

    return interpolated
