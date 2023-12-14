import copy
import os
from numba import jit
import numpy as np
from PIL import Image, ImageDraw
from typing import Dict, Iterator, List, Tuple, Union
from .ct import CT
from .grid import Grid
from .rounding import drop_precision


PyFTPP_slice = Dict[str, List[Dict[str, float]]]


class Structure:

    def __init__(self, name: str, points: np.ndarray, grid: Union[None, Grid] = None):
        self._name = name
        self._points = points
        self._grid = grid
        self._contours = None
        self._mask = None

    @property
    def fiona_slices(self) -> List[PyFTPP_slice]:
        '''
        Returns
        -------
        List[PyFTPP_slice]
            The structure in the format of the Fiona standalone executable.
        '''
        if self._grid is None:
            raise RuntimeError(
                f'Need to set a grid to convert structure {self.name} '
                'to Fiona slices.'
            )
        if self._contours is None:
            self._contours = numpy_to_pyftpp_contour(
                self.points, self._grid.spacing[2], self._grid.origin[2]
            )
        return self._contours

    @property
    def mask(self) -> np.ndarray:
        '''
        A voxel is considered part of the structure if its center lies
        on or within a contour.

        Voxels that are intersected are not considered part of the structure
        unless their center lies within a contour.
        '''
        if self._mask is None:
            self._mask = binary_mask_from_contours(
                self._points,
                self.grid.shape,
                self.grid.spacing,
                self.grid.origin
            )
        return self._mask.copy()

    @property
    def name(self) -> str:
        return copy.copy(self._name)
    
    @property
    def points(self) -> np.ndarray:
        return self._points.copy()
    
    @property
    def grid(self) -> Grid:
        return self._grid
    
    @grid.setter
    def grid(self, grid: Grid):
        self._grid = grid
    
    #@jit(nopython=False)
    def expand(self, distance: float):
        raise NotImplementedError()
        # Super slow and not tested on a real-world structure!!!
        assert False
        if self._grid is None:
            raise RuntimeError(
                f'Need to set a grid to expand structure {self.name}.'
            )

        expanded_points = []
        for points in split_slices(self._points, self.grid.spacing[2], self.grid.origin[2]):
            lines = _polygon_to_lines(self._points[:, :2])
            _fill_normal_vectors(lines)
            new_lines = _translate_lines_along_normal(lines, distance)
            new_points = _lines_to_polygon(new_lines)
            new_points = np.hstack([
                new_points,
                np.ones(new_points.shape[0]).reshape(-1, 1) * points[0, 2],
            ])
            expanded_points.append(new_points)
        expanded_points = np.vstack(expanded_points)
        return Structure(self.name, expanded_points, self.grid)


class StructureSet:

    def __init__(self, structures: Iterator[Structure], ct: CT = None):
        self._structures = {}
        for s in structures:
            self._structures[s.name] = s
        self._ct = ct
    
    def add(self, structure: Structure):
        self._structures[structure.name] = structure

    def __getitem__(self, key: str) -> Structure:
        return self._structures[key]
    
    def __iter__(self) -> Structure:
        for name in self.structure_names:
            yield self[name]

    @property
    def structure_names(self) -> Iterator[str]:
        for name in self._structures.keys():
            yield name

    @property
    def ct(self) -> Union[None, CT]:
        return self._ct

    def to_dicom(
        self,
        study_instance_UID: str,
        frame_of_reference_UID: str,
        base_ct_series_instance_UID: str,
        structureset_series_instance_UID: str,
        patient_ID: str,
        patient_name: str,
        ct_datasets = None,
        drop_precision: bool = True,
    ):
        # pydicom is an optional module: Import it only when it is needed.
        import pydicom

        if ct_datasets is None:
            ct_datasets = self.ct.to_dicom(
                study_instance_UID,
                frame_of_reference_UID,
                ct_series_instance_UID,
                patient_ID,
                patient_name,
            )

        ds = pydicom.Dataset()
        ds.Modality = 'RTSTRUCT'
        # RT Structure Set Storage
        ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
        ds.StructureSetLabel = 'StructureSet'
        ds.StructureSetName = 'StructureSet'
        ds.ApprovalStatus = 'UNAPPROVED'
        ds.AccessionNumber = ''
        ds.StudyID = ''
        ds.SeriesNumber = '0'

        referenced_frame_of_reference = pydicom.Dataset()
        referenced_frame_of_reference.FrameOfReferenceUID = frame_of_reference_UID
        referenced_frame_of_reference.RTReferencedStudySequence = build_referenced_study_sequence(
            study_instance_UID,
            base_ct_series_instance_UID,
            ct_datasets,
        )

        ds.ReferencedFrameOfReferenceSequence = [referenced_frame_of_reference]

        roi_seq = []
        for i, structure in enumerate(self):
            roi_ds = pydicom.Dataset()
            roi_ds.ROINumber = i
            roi_ds.ReferencedFrameOfReferenceUID = frame_of_reference_UID
            roi_ds.ROIName = structure.name
            roi_ds.ROIGenerationAlgorithm = 'MANUAL'
            roi_seq.append(roi_ds)
        ds.StructureSetROISequence = roi_seq

        roi_contour_sequence = []
        for i, structure in enumerate(self):
            contours = build_contour_sequence(structure, self.ct, ct_datasets, drop_precision)
            roi_contour = pydicom.Dataset()
            roi_contour.ReferencedROINumber = i
            roi_contour.ContourSequence = contours
            roi_contour_sequence.append(roi_contour)
        ds.ROIContourSequence = roi_contour_sequence

        observations = []
        for i, structure in enumerate(self):
            observation = pydicom.Dataset()
            observation.ObservationNumber = i
            observation.ReferencedROINumber = i
            observation.ROIObservationLabel = structure.name
            observation.ROIInterpreter = ''
            
            n = structure.name.lower()
            is_target = ('ctv' in n) or ('gtv' in n) or ('ptv' in n)
            observation.RTROIInterpretedType = 'PTV' if is_target else 'ORGAN'

            observations.append(observation)
        ds.RTROIObservationsSequence = observations

        ds.file_meta = pydicom.dataset.FileMetaDataset()
        ds.file_meta.magic = b'DICM'
        # Use the default "DICOM Implicit VR Little Endian Transfer Syntax"
        # https://dicom.nema.org/dicom/2013/output/chtml/part05/chapter_10.html#sect_10.1
        ds.file_meta.TransferSyntaxUID = '1.2.840.10008.1.2'
        # RT Structure Set Storage
        ds.file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
        ds.file_meta.MediaStorageSOPInstanceUID = structureset_series_instance_UID

        ds.PatientName = patient_name
        ds.PatientID = patient_ID
        ds.StudyInstanceUID = study_instance_UID
        ds.SeriesInstanceUID = structureset_series_instance_UID
        ds.SOPInstanceUID = structureset_series_instance_UID
        ds.FrameOfReferenceUID = frame_of_reference_UID

        return ds


def build_referenced_series_sequence(base_ct_series_instance_UID, ct_datasets):
    # pydicom is an optional module: Import it only when it is needed.
    import pydicom

    ref_series = pydicom.Dataset()
    ref_series.SeriesInstanceUID = base_ct_series_instance_UID

    contour_image_sequence = []
    for ct_ds in ct_datasets:
        ds = pydicom.Dataset()
        ds.ReferencedSOPClassUID = ct_ds.SOPClassUID
        ds.ReferencedSOPInstanceUID = ct_ds.SOPInstanceUID
        contour_image_sequence.append(ds)
    ref_series.ContourImageSequence = contour_image_sequence
    
    return [ref_series]


def build_referenced_study_sequence(study_instance_UID, base_ct_series_instance_UID, ct_datasets):
    # pydicom is an optional module: Import it only when it is needed.
    import pydicom

    ref_study = pydicom.Dataset()
    ref_study.ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3' # RT Structure Set Storage
    ref_study.ReferencedSOPInstanceUID = study_instance_UID
    ref_study.RTReferencedSeriesSequence = build_referenced_series_sequence(base_ct_series_instance_UID, ct_datasets)
    
    return [ref_study]


def build_contour_sequence(structure, ct, ct_datasets, decrease_precision: bool = True):
    # pydicom is an optional module: Import it only when it is needed.
    import pydicom

    slices = split_slices(structure.points, ct.spacing[2], ct.origin[2])

    contours = []
    for points in slices:
        z = points[0, 2]
        assert np.isclose(points[:, 2], z, atol=1e-6).all()
        slice_ind = int(np.round((z - ct.origin[2] - ct.spacing[2] / 2) / ct.spacing[2]))

        # Avoid explicitly closing the contour: DICOM takes care of this because
        # of the CLOSED_PLANAER geometric type (see below).
        if np.isclose(points[0, :], points[-1, :], atol=1e-6).all():
            points = points[:-1, :]

        ct_dataset = ct_datasets[slice_ind]

        contour = pydicom.Dataset()

        contour_image = pydicom.Dataset()
        contour_image.ReferencedSOPClassUID = ct_dataset.SOPClassUID
        contour_image.ReferencedSOPInstanceUID = ct_dataset.SOPInstanceUID

        contour.ContourImageSequence = [contour_image]
        contour.ContourGeometricType = 'CLOSED_PLANAR'
        contour.NumberOfContourPoints = points.shape[0]
        # We use cm because the Fiona standalone does, but DICOM expects mm.
        points = points * 10
        # DICOM export to Fiona makes problems if the numbers are not exactly
        # the slice spacing.
        if decrease_precision:
            points = drop_precision(points)
        contour.ContourData = points.flatten().tolist()
        contours.append(contour)

    return contours


class FionaStructureSet(StructureSet):

    def __init__(self, structure_dir: str, ct: CT = None):
        self._structure_dir = structure_dir
        self._ct = ct

        self._structure_names = []
        for file in os.listdir(self._structure_dir):
            if not file.endswith('.npy'):
                continue
            name = file.replace('.npy', '')
            self._structure_names.append(name)
        self._structures = {}

    def add(self, structure: Structure):
        self._structure_names.append(structure.name)
        self._structures[structure.name] = structure

    @property
    def structure_names(self) -> Iterator[str]:
        for name in self._structure_names:
            yield name

    def __getitem__(self, key: str) -> Structure:
        # Lazy loading the structures: Only load them when they are requested,
        # then keep them cached in memory.
        if key not in self._structures.keys():
            path = f'{self._structure_dir}/{key}.npy'
            grid = Grid(
                self.ct.data.shape, self.ct.spacing, self.ct.origin
            )
            self._structures[key] = load_array_structure(path, key, grid)
        return self._structures[key]


def load_array_structure(path: str, name: str, grid: Union[None, Grid] = None) -> Structure:
    if os.path.exists(path):
        points = np.load(path)
        return Structure(name, points, grid)
    else:
        raise RuntimeError(
            f'Could not load structure {name} from path "{path}".'
        )


################################
# Actual conversion functions.
################################
def numpy_to_pyftpp_contour(arr: np.ndarray,
                            slice_spacing: float,
                            z_origin: float) -> List[PyFTPP_slice]:
    '''
    Convert a numpy array of points for a structure into the PyFTPP format.

    The PyFTPP format is documented in the docstring for the class
    OptimizationSettings on pyftpp.config.

    Parameters
    ----------
    arr: np.ndarray
        Array containing all the points in the structure.
        Shape: (n_points, 3)
    slice_spacing: float
        Spacing between CT slices.
    z_origin: float
        z-coordinate of the CT origin.

    Returns
    -------
    list of PyFTPP_slice
        The structure representation in the PyFTPP format.
    '''
    slices = split_slices(arr, slice_spacing, z_origin)

    to_ret = []
    for sl in slices:
        points_in_slice = [{'x': p[0], 'y': p[1], 'z': p[2]} for p in sl]
        to_ret.append({'points': points_in_slice})

    return to_ret


def split_slices(points: np.ndarray, slice_spacing: float, z_origin: float):
    '''
    Take an array and split it into a list of sub-arrays based on the
    z-coordinate.

    Parameters
    ----------
    points: np.ndarray
        Array representing a set of points. The points are assumed to represent
        a structure contour, i. e. subsequent points within the same slice are
        supposed to be connected by lines. Within a slice, all points share the
        same z-coordinate. The points are supposed to be ordered such that
        the points belonging to the same slice form a contiguous subarray.
        Shape: (n_points, 3)
    slice_spacing: float
        Spacing between slices. This is needed when there was a registration
        and the points need to be mapped to slices first.
    z_origin: float
        z-component of the coordinate system origin. The slice with index 0
        has this z-coordinate.

    Returns
    -------
    list of np.ndarray
        Each element in the list represents the points in a single slice.
        Empty slices are omitted.
    '''
    slices = []

    z_indices = np.rint((points[:, 2] - z_origin) / slice_spacing).astype(np.int64)

    i_start = 0
    i_end = 0
    while i_end < points.shape[0]:
        if z_indices[i_end] != z_indices[i_start]:
            # Save the old slice and start a new one.
            slices.append(points[i_start:i_end, :])
            i_start = i_end
        i_end += 1

    assert i_end == points.shape[0]

    if i_end - i_start > 1:
        # Add the last slice.
        slices.append(points[i_start:i_end, :])

    return slices


@jit(nopython=True)
def shortest_distance_to_line(x, a, b):
    d_ab = b - a
    d_xa = x - a
    
    #assert not np.isclose(d_ab, 0.).all()
    
    # TODO: Why squared???
    t = np.dot(d_xa, d_ab) / np.linalg.norm(d_ab)**2
    
    if t < 0:
        return np.linalg.norm(x - a)
    elif t < 1.:
        x_proj = a + t * d_ab
        return np.linalg.norm(x - x_proj)
    else:
        return np.linalg.norm(x - b)


@jit(nopython=True)
def shortest_distance_to_polygon(x, contour):
    N = contour.shape[0]
    minimum_distance = np.inf

    for i in range(N):
        p0 = contour[i, :2]
        p1 = contour[(i+1)%N, :2]

        d = shortest_distance_to_line(x[:2], p0, p1)
        minimum_distance = min(minimum_distance, d)

    return minimum_distance


@jit(nopython=True)
def contour_line_intersections(line_value, contour, axis=0):
    '''
    Calculate intersection points of the contour with an axis-parallel line.

    Returns a sorted list of intersection points in ascending order along
    the axis.
    '''
    N = contour.shape[0]
    intersection_points = []
    for i in range(N):
        p0 = contour[i, :2]
        p1 = contour[(i+1)%N, :2]

        if (abs(p0[axis] - p1[axis]) <= 1e-6):
            # Axis parallel line.
            if abs(p0[axis] - line_value) <= 1e-6:
                # The line is close to the scan line.
                intersection_points.append(p0)
                intersection_points.append(p1)
            # Avoid division by zero for parallel lines.
            continue

        t = (line_value - p0[axis]) / (p1[axis] - p0[axis])

        if 0. <= t and t <= 1.:
            intersection_points.append(p0 + t * (p1 - p0))

    return intersection_points


@jit(nopython=True)
def calculate_distance_mask_for_contour(shape, spacing, contour, x_indices, y_indices):
    '''Assumes zero origin.'''
    distance = np.zeros((len(x_indices), len(y_indices)), dtype=np.float64)

    for ni, i in enumerate(x_indices):
        x = i*spacing[0]
        intersection_points = contour_line_intersections(x, contour, axis=0)

        for nj, j in enumerate(y_indices):
            y = j*spacing[1]

            crossed_intersections = [p for p in intersection_points if p[1] <= y]
            is_inside = len(crossed_intersections) % 2 == 1

            p = np.array([x, y])
            d = shortest_distance_to_polygon(p, contour)
            if is_inside:
                d = -d
            distance[ni, nj] = d
    return distance


@jit(nopython=True)
def calculate_distance_mask(shape, spacing, contours):
    '''Assumes zero origin.'''
    mask = np.zeros(shape)
    for contour in contours:
        z = contour[0, 2]
        iz = int(z // spacing[2])
        slice_mask = calculate_distance_mask_for_contour(
            shape[:2],
            spacing,
            contour,
            np.arange(shape[0], dtype=np.int32),
            np.arange(shape[1], dtype=np.int32),
        )
        mask[:, :, iz] = slice_mask
    return mask


def calculate_binary_mask(shape, spacing, contours, max_distance=0.1):
    '''
    Calculate the binary mask for the given contours and grid.

    Assumes zero origin.

    Parameters
    ----------
    shape:
        Shape of the grid on which to calculate the binary mask. Shape (3,).
    spacing:
        Spacing of the grid on which to calculate the binary mask. Shape (3,).
    contours:
        Iterable of nx3 arrays representing a structure.
    max_distance: float
        A voxel is considered to be a part of the structure if its center lies
        within this distance from any contour. (in-slice distance; cm)
    '''
    mask = np.zeros(shape)
    for contour in contours:
        z = contour[0, 2]
        iz = int(z // spacing[2])
        
        # Don't need to check for all grid positions, just those close
        # to the structure. Add a safety margin to make sure we don't
        # miss any voxel.
        safety_margin = 2 + int(np.ceil(np.max(max_distance / spacing)))

        x_min = int(np.round(np.min(contour[:, 0]) / spacing[0])) - safety_margin
        x_max = int(np.round(np.max(contour[:, 0]) / spacing[0])) + safety_margin
        x_min, x_max = np.clip([x_min, x_max], 0, shape[0]-1)

        y_min = int(np.round(np.min(contour[:, 1]) / spacing[1])) - safety_margin
        y_max = int(np.round(np.max(contour[:, 1]) / spacing[1])) + safety_margin
        y_min, y_max = np.clip([y_min, y_max], 0, shape[1]-1)

        slice_distance_mask = calculate_distance_mask_for_contour(
            shape[:2],
            spacing,
            contour,
            np.arange(x_min, x_max+1, dtype=np.int32),
            np.arange(y_min, y_max+1, dtype=np.int32),
        )
        mask[x_min:x_max+1, y_min:y_max+1, iz] = slice_distance_mask <= max_distance
    return mask


def binary_mask_from_contours(points: np.ndarray,
                              grid_size: Tuple[int, int, int],
                              grid_spacing: np.ndarray,
                              grid_origin: np.ndarray) -> np.ndarray:
    '''
    Calculates the binary 3D mask for an entire structure.

    Parameters
    ----------
    points: np.ndarray
        Contours to convert to a binary mask
        Shape: (n, 3)
    grid_size: tuple of int
        Shape: (3,)
    grid_spacing: np.ndarray
        Shape: (3,)
    grid_origin: np.ndarray
        Shape: (3,)

    Returns
    -------
    np.ndarray
        Binary mask of shape grid_size.
    '''
    points = points - grid_origin
    contours = split_slices(points, grid_spacing[2], grid_origin[2])
    return calculate_binary_mask(
        grid_size,
        grid_spacing,
        contours,
        max_distance=0.,
    )



class _Line2D:
    '''
    Represents a 2D line of the following form: x(t) = x1 + t * (x2 - x1),
    where x, x1, x2 are 2D arrays and t is a scalar.
    
    The line has a start (t=0) and an end (t=1) point.
    '''
    
    def __init__(self, x1, x2):
        self.x1 = x1
        self.x2 = x2
        self.d = x2 - x1
        self._normal_vector = None
    
    def intersection_t(self, other):
        '''
        We obtain the intersection point between both lines by setting the
        equations equal.
        '''
        A = np.hstack([
            self.d.reshape(-1, 1),
            -other.d.reshape(-1, 1)
        ])
        b = other.x1 - self.x1
        return np.linalg.solve(A, b)[0]

    def __eq__(self, other):
        return all((
            np.isclose(self.x1, other.x1).all(),
            np.isclose(self.x2, other.x2).all(),
        ))

    def __call__(self, t):
        return self.x1 + t * self.d

    @property
    def normal_vector(self):
        return self._normal_vector
    
    @normal_vector.setter
    def normal_vector(self, n):
        self._normal_vector = n
    
    def __str__(self):
        n = 'None' if self.normal_vector is None else self.normal_vector.tolist()
        return f'{self.x1.tolist()} -- {self.x2.tolist()}; normal vector {n}'


def _polygon_to_lines(points):
    lines = []
    for i in range(points.shape[0]):
        x1 = points[i, :]
        x2 = points[(i+1)%points.shape[0], :]
        lines.append(_Line2D(x1, x2))
    return lines


def _lines_to_polygon(lines):
    new_points = []
    for i in range(len(lines)):
        l1 = lines[i]
        l2 = lines[(i+1)%len(lines)]
        try:
            t = l1.intersection_t(l2)
            new_points.append(l1(t))
        except np.linalg.LinAlgError:
            pass
    return np.vstack(new_points)


@jit(nopython=True)
def _intersecting_points(line_of_interest, lines):
    x1 = np.array([
        line_of_interest.x1[0],
        line_of_interest.x1[1] + line_of_interest.d[1] / 2,
    ])
    x2 = x1 + np.array([1, 0])
    scan_line = _Line2D(x1, x2)

    intersection_points = []
    for line in lines:
        if line is line_of_interest:
            continue
        try:
            if np.isclose(line.d, scan_line.d).all():
                # Skip lines parallel to the horizontal line.
                continue
            t = line.intersection_t(scan_line)
            if 0 <= t and t <= 1:
                intersection_points.append(line(t))
        except np.linalg.LinAlgError:
            continue

    return sorted(intersection_points, key=lambda p: p[1])


def _is_left_halfspace_inside(line, intersection_points):
    cnt = 0
    for p in intersection_points:
        if np.isclose(p, line.x1).all():
            break
        cnt += 1
    left_inside = (cnt % 2) == 0
    return left_inside


def _calculate_normal_vector(line, left_inside):
    normal_vector = np.abs(np.array([line.d[1], line.d[0]]))
    normal_vector /= (np.linalg.norm(normal_vector) + 1e-6)
    if left_inside:
#         print('left inside')
        # See documentation and use that d[0] == d_2 and d[1] == d_1
        normal_vector[1] *= -1 * np.sign(line.d[0])
        normal_vector[0] *= -1 * np.sign(line.d[1])
    else:
#         print('left outside')
        normal_vector[1] *= -1 * np.sign(line.d[0])
        normal_vector[0] *= np.sign(line.d[1])
    return normal_vector


def _fill_normal_vectors(lines):
    for line_of_interest in lines:
        intersection_points = _intersecting_points(line_of_interest, lines)
        left_inside = _is_left_halfspace_inside(
            line_of_interest,
            intersection_points
        )
        line_of_interest.normal_vector = _calculate_normal_vector(line_of_interest, left_inside)


def _translate_lines_along_normal(lines, delta):
    new_lines = []
    for line in lines:
        new_lines.append(_Line2D(
            line.x1 + delta * line.normal_vector,
            line.x2 + delta * line.normal_vector,
        ))
    return new_lines


class CompositeStructure(Structure):

    def __init__(self, name: str, structures: List[Structure]):
        if len(structures) < 2:
            raise RuntimeError('Need at least two structures to compose!')

        self._name = name
        self._structures = structures
        for s in structures[1:]:
            if s.grid != self._structures[0].grid:
                raise RuntimeError('Structures must have the same grid to be composed!')
        self._points = None
        self._mask = None
        self._grid = None

    @property
    def name(self) -> str:
        return copy.copy(self._name)

    @property
    def points(self) -> np.ndarray:
        if self._points is None:
            self._points = self._build_points()
        return self._points

    @property
    def fiona_slices(self) -> List[PyFTPP_slice]:
        all_contours = []
        for s in self._structures:
            contours = s.fiona_slices
            all_contours.extend(contours)
        return all_contours

    @property
    def mask(self) -> np.ndarray:
        if self._mask is None:
            self._mask = self._build_mask()
        return self._mask
    
    @property
    def grid(self) -> Grid:
        return self._grid
    
    @grid.setter
    def grid(self, grid: Grid):
        self._grid = grid

    def _build_points(self) -> np.ndarray:
        points = []
        for s in self._structures:
            points.append(s.points)
        points = sorted(points, key=lambda p: p[0, 2])
        return np.vstack(points)

    def _build_mask(self) -> np.ndarray:
        mask = self._structures[0].mask
        for s in self._structures[1:]:
            mask = np.logical_or(mask, s.mask)
        return mask


@jit(forceobj=True)
def expand_mask(mask):
    new_mask = np.zeros_like(mask)
    for i, cell in enumerate(mask.flatten()):
        if not cell:
            continue

        ix, iy, iz = np.unravel_index(i, shape=mask.shape)

        for x_delta in [-1, 0, 1]:
            for y_delta in [-1, 0, 1]:
                for z_delta in [-1, 0, 1]:
                    # Clamp indices to be valid, i. e. 0 <= i <= size-1.
                    index = (
                        min(max(ix + x_delta, 0), mask.shape[0]-1),
                        min(max(iy + y_delta, 0), mask.shape[1]-1),
                        min(max(iz + z_delta, 0), mask.shape[2]-1),
                    )
                    new_mask[index] = True
    return new_mask


class MaskExpansion(Structure):
    '''
    Defines a new structure whose binary mask is the one of
    a preexisting structure plus a dilation.

    The contour points are NOT expanded!!!
    '''

    def __init__(self, structure, n_voxels, new_name=None):
        self._structure = structure
        if new_name is None:
            self._name = copy.copy(structure.name)
        else:
            self._name = new_name
        self._grid = copy.copy(structure.grid)
        self._original_mask = copy.copy(structure.mask)
        self._mask = None
        self._n_voxels = n_voxels

    def _build_mask(self):
        self._mask = self._original_mask
        for _ in range(self._n_voxels):
            self._mask = expand_mask(self._mask)

    @property
    def fiona_slices(self) -> List[PyFTPP_slice]:
        return self._structure.fiona_slices

    @property
    def mask(self) -> np.ndarray:
        if self._mask is None:
            self._build_mask()
        return self._mask

    @property
    def name(self) -> str:
        return self._name

    @property
    def points(self) -> np.ndarray:
        return self._structure.points

    @property
    def grid(self) -> Grid:
        return self._grid

    @grid.setter
    def grid(self, grid: Grid):
        self._grid = grid
