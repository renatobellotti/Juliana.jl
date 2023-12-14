from collections import Counter
import cv2 as cv
from numba import jit
import numpy as np
from .structure import Structure


@jit
def count_neighbours(binary_img):
    '''
    Count the number of active (1) neighbours in a binary
    2D image. 4-connectivity is used.

    Returns
    -------
    np.ndarray
        Image containing the number of active neighbours for
        each pixel. Same shape as binary_img.
    '''
    count = np.empty(binary_img.shape, dtype=np.int64)
    for i in range(binary_img.shape[0]):
        for j in range(binary_img.shape[1]):
            if i == 0:
                top = 0
            else:
                top = binary_img[i-1, j]
            if i == binary_img.shape[1]-1:
                bottom = 0
            else:
                bottom = binary_img[i+1, j]
            if j == 0:
                left = 0
            else:
                left = binary_img[i, j-1]
            if j == binary_img.shape[1]-1:
                right = 0
            else:
                right = binary_img[i, j+1]
            count[i, j] = top + bottom + left + right
    return count


def extract_biggest_connected_component(binary_img):
    kernel = np.ones((1, 1), dtype=np.uint8)
    binary_img = cv.dilate(cv.erode(binary_img.astype(np.int8), kernel), kernel)
    retval, component_labels = cv.connectedComponents(binary_img, connectivity=8)
    counter = Counter(component_labels.flatten().tolist())
    # Ignore background.
    counter[0] = -1
    return (component_labels == counter.most_common()[0][0])


def find_whole_body_contour(ct):
    normal_tissue_mask = np.logical_and(ct.data >= -850, ct.data < 2000)

    # Find the border cells (called "outside") and interior cells (called "inside").
    neighbour_count = []
    for iz in range(normal_tissue_mask.shape[2]):
        slice_neighbour_count = count_neighbours(normal_tissue_mask[:, :, iz])
        neighbour_count.append(slice_neighbour_count)
    neighbour_count = np.stack(neighbour_count, axis=-1)

    inside_mask = (neighbour_count == 4)
    outside_mask = np.logical_xor(normal_tissue_mask, inside_mask)

    # Clean up: Remove isolated voxels, e. g. artifacts from the CT.
    is_isolated = (neighbour_count == 0)
    outside_mask = np.logical_and(outside_mask, np.logical_not(is_isolated))

    # We do not care about interior contours, so we focus only on the biggest
    # connected component.
    extracted_outside_contour = []
    for iz in range(outside_mask.shape[2]):
        extracted_outside_contour.append(
            extract_biggest_connected_component(outside_mask[:, :, iz])
        )
    outside_mask = np.stack(extracted_outside_contour, axis=-1)

    # Use OpenCV to convert the binary mask of the whole body contour
    # to polygons.
    # For some reason, the last slice contains the entire image...
    # Remove the last slice. (-1 in the range)
    all_points = []
    for iz in range(outside_mask.shape[2]-1):
        z = ct.origin[2] + iz*ct.spacing[2]

        img = outside_mask[:, :, iz].astype(np.uint8)
        contours = cv.findContours(
            img,
            mode=cv.RETR_EXTERNAL,
            method=cv.CHAIN_APPROX_NONE,
        )
        assert len(contours[0]) == 1
        # I don't know exactly why [0][0] is needed, but it seems to lead
        # to good results...
        points = contours[0][0][:, 0, :]

        points = ct.origin[:2] + points * ct.spacing[:2]
        tmp = points[:, 0].copy()
        points[:, 0] = points[:, 1]
        points[:, 1] = tmp

        points_3d = np.empty((points.shape[0], 3), dtype=points.dtype)
        points_3d[:, :2] = points
        points_3d[:, 2] = z
        all_points.append(points_3d)

    all_points = np.vstack(all_points)

    # Build a structure from the polygons.
    whole_body = Structure('WHOLE_BODY', all_points, ct.grid)

    return whole_body
