import copy
from typing import List, Set, Union
import numpy as np
import pydicom
from pydicom.pixel_data_handlers.util import apply_modality_lut
from .ct import CT
from .dose import Dose
from .structure import Structure, StructureSet


# Save DICOM files.
def export_to_dicom(ct,
                    structureset,
                    output_dir,
                    study_instance_UID,
                    new_patient_ID,
                    dose_distributions={},
                    decrease_precision: bool = True):
    frame_of_reference_UID = f'{study_instance_UID}.0'
    ct_series_instance_UID = f'{study_instance_UID}.1'
    structureset_series_instance_UID = f'{study_instance_UID}.2'
    patient_name = f'{new_patient_ID}^{new_patient_ID}'

    ct_datasets = ct.to_dicom(
        study_instance_UID,
        frame_of_reference_UID,
        ct_series_instance_UID,
        new_patient_ID,
        patient_name,
        decrease_precision,
    )
    for ds in ct_datasets:
        ds.SeriesDescription = 'CT Images'
        ds.save_as(
            f'{output_dir}/CT.{ds.InstanceNumber}.dcm',
            write_like_original=False,
        )

    ds = structureset.to_dicom(
        study_instance_UID,
        frame_of_reference_UID,
        ct_series_instance_UID,
        structureset_series_instance_UID,
        new_patient_ID,
        patient_name,
        ct_datasets,
        decrease_precision,
    )
    ds.SeriesDescription = 'StructureSet'
    ds.save_as(
        f'{output_dir}/RS.{structureset_series_instance_UID}.dcm',
        write_like_original=False,
    )

    for i, (label, dose) in enumerate(dose_distributions.items()):
        dose_series_instance_UID = f'{study_instance_UID}.{3+i}'
        dicom_dose = dose.to_dicom_dataset(
            study_instance_UID,
            frame_of_reference_UID,
            dose_series_instance_UID,
            new_patient_ID,
            patient_name,
        )
        dicom_dose.SeriesDescription = label
        dicom_dose.save_as(
            f'{output_dir}/RD_plan_{label}.dcm',
            write_like_original=False,
        )


###########################
# Load DICOM files.
###########################

def load_structureset(path: str, ct: Union[None, CT] = None) -> StructureSet:
    dicom_structure_set = DicomStructureSet(path)

    structures = []
    for s in dicom_structure_set.structure_names:
        if len(dicom_structure_set.point_sets[s]) == 0:
            points = np.empty((0, 3))
        else:
            points = np.vstack(dicom_structure_set.point_sets[s])

        # Convert mm -> cm.
        points = points / 10.

        # Avoid problems when saving the structure.
        s = s.replace('/', '_')

        grid = None
        if ct is not None:
            grid = ct.grid

        structures.append(Structure(s, points, grid))
    return StructureSet(structures, ct)



def dicom_to_dict(dataset, extract_pixel_array=False, recursive=False):
    '''
    Parameters
    ----------
    dataset: pydicom.Dataset
        The dataset to convert into a Python dict.
    '''
    keys = dataset.keys()
    names = [dataset[k].name for k in keys]

    entry = {}

    for name, key in zip(names, keys):
        if name == 'Pixel Data':
            continue
        val = dataset[key].value

        if recursive and type(val) == pydicom.sequence.Sequence:
            val_converted = []
            for v in val:
                v_conv = dicom_to_dict(v, False, recursive)
                val_converted.append(v_conv)
            val = val_converted

        entry[name] = val

    if extract_pixel_array:
        entry['Pixel Data'] = dataset.pixel_array

    return entry


class DicomCTImage:
    '''
    Class representing a DICOM CT volume.

    This class assumes that all the slices have the same spacing
    and physical image position in the Reference Coordinate System.
    (see DICOM standard for how the coordinate system is defined)
    '''

    def __init__(self, paths, use_cm: bool = True):
        '''
        Parameters
        ----------
        paths: list of str
            Files in which the CT slices are stored. Each file is a DICOM file.
            The list is assumed to be ordered along the z-axis.
        use_cm: bool
            Whether to express the origin and spacing in cm.
            Otherwise mm is used.
            Default: True
        '''
        self._use_cm = use_cm
        self._slices = []

        image_positions = []
        image_orientations_X = []
        image_orientations_Y = []
        row_spacings = []
        col_spacings = []

        # load the slices
        for path in paths:
            CT_ds = pydicom.read_file(path)
            CT_d = dicom_to_dict(CT_ds, extract_pixel_array=True)
            image_positions.append(CT_d['Image Position (Patient)'])
            arr = CT_d['Pixel Data']
            arr = apply_modality_lut(arr, CT_ds)
            self._slices.append(arr)

            X = [int(x) for x in CT_d['Image Orientation (Patient)'][:3]]
            X = np.array(X)

            Y = [int(y) for y in CT_d['Image Orientation (Patient)'][3:]]
            Y = np.array(Y)

            image_orientations_X.append(X)
            image_orientations_Y.append(Y)

            # Unit: mm
            row_spacings.append(CT_d['Pixel Spacing'][0])
            col_spacings.append(CT_d['Pixel Spacing'][1])

        image_positions = np.vstack(image_positions)
        image_orientations_X = np.vstack(image_orientations_X)
        image_orientations_Y = np.vstack(image_orientations_Y)
        row_spacings = np.array(row_spacings)
        col_spacings = np.array(col_spacings)
        self._slices = np.stack(self._slices, axis=-1)
        self._slices = np.swapaxes(self._slices, 0, 1)

        # Make sure that all the CT images have the same x and y positions.
        assert (image_positions[:, 0] == image_positions[0, 0]).all()
        assert (image_positions[:, 1] == image_positions[0, 1]).all()

        self._image_position = image_positions[0, :].copy()
        self._slice_depths = image_positions[:, 2]

        # There should be no depth-duplicates.
        assert len(self._slice_depths) == len(set(self._slice_depths))

        # Make sure that all the slices have the same spacing.
        assert (row_spacings == row_spacings[0]).all()
        assert (col_spacings == col_spacings[0]).all()

        self._row_spacing = row_spacings[0]
        self._col_spacing = col_spacings[0]

        # Make sure that all the image orientations are the same.
        assert (image_orientations_X[0, :] == image_orientations_X).all()
        assert (image_orientations_Y[0, :] == image_orientations_Y).all()

        self._image_orientation_X = image_orientations_X[0, :]
        self._image_orientation_Y = image_orientations_Y[0, :]

        if self._use_cm:
            self._image_position /= 10.
            self._row_spacing /= 10.
            self._col_spacing /= 10.
            self._slice_depths /= 10.

    @property
    def voxel_matrix(self):
        '''
        Returns
        -------
        m: array
            Shape: (n_pixels_x, n_pixels_y, n_pixels_z)
        '''
        return np.copy(self._slices)

    @property
    def row_spacing(self):
        '''Spacing between rows.'''
        return self._row_spacing

    @property
    def col_spacing(self):
        '''Spacing between columns.'''
        return self._col_spacing

    @property
    def slice_spacing(self):
        '''Spacing between slices.'''
        # We have one separation less than slices
        separations = np.empty(self._slice_depths.shape[0] - 1)

        for i in range(1, self._slice_depths.shape[0]):
            separations[i-1] = self._slice_depths[i] - self._slice_depths[i-1]

        assert np.isclose(separations, separations[0], atol=1e-6).all()

        return separations[0]

    @property
    def image_position(self):
        '''
        Returns the position of the top left corner of the CT in units of mm.

        The position of the middle point of the top left pixel is returned.
        '''
        return self._image_position

    @property
    def image_orientation_x(self):
        '''
        First half of the "Image Orientation (Patient)" DICOM attribute,
        which indicates how the coordinate system x-axis is oriented.
        '''
        return self._image_orientation_X

    @property
    def image_orientation_y(self):
        '''
        First half of the "Image Orientation (Patient)" DICOM attribute,
        which indicates how the coordinate system y-axis is oriented.
        '''
        return self._image_orientation_Y

    @property
    def ct(self):
        return CT(
            self.voxel_matrix,
            np.array([self.row_spacing, self.col_spacing, self.slice_spacing]),
            self.image_position,
        )


class DicomStructureSet:
    '''
    Class representing a DICOM RT Structure Set.
    '''

    def __init__(self, file: str, use_cm=True):
        '''
        Parameters
        ----------
        file: str
            Path to a DICOM file containing the structure set.
        use_cm: bool
            Whether to use cm as the unit for the contour point coordinates.
            Otherwise, mm is used.
            The unit must match the unit for the CT origin and spacing.
            Default: True
        '''
        self._use_cm = use_cm
        ds = pydicom.read_file(file)
        self._d = dicom_to_dict(ds, recursive=True)

        # Collect names of the structures (regions of interest).
        self._structure_names = []
        for roi in self._d['Structure Set ROI Sequence']:
            self._structure_names.append(roi['ROI Name'])

        # Collect the points belonging to each structure and their color.
        self._pointsets = {}
        self._structure_colors = {}
        for ID, name in enumerate(self._structure_names):
            color = None
            all_points = []
            if 'Contour Sequence' in self._d['ROI Contour Sequence'][ID]:
                contour_sequence = self._d['ROI Contour Sequence'][ID]['Contour Sequence']

                color = self._d['ROI Contour Sequence'][ID].get('ROI Display Color', '')

                # Collect the points in each slice.
                for contour in contour_sequence:
                    n_points = contour['Number of Contour Points']

                    points = np.empty((n_points, 3))
                    for i, val in enumerate(contour['Contour Data']):
                        points[i // 3, i % 3] = val

                    all_points.append(points)

            self._structure_colors[name] = color
            self._pointsets[name] = all_points

    @property
    def n_structures(self) -> int:
        '''
        Returns
        -------
        n_structures: int
        '''
        return len(self._d['ROI Contour Sequence'])

    @property
    def structure_names(self) -> Set[str]:
        '''
        Returns
        -------
        structure_names: set of str
            alphabetically sorted list of structure names
        '''
        return set(copy.deepcopy(self._structure_names))

    @property
    def structure_colors(self) -> dict:
        return copy.deepcopy(self._structure_colors)

    @property
    def point_sets(self) -> dict:
        '''
        Point sets.

        Returns
        -------
        point_sets: dict
            The keys are the names of the structures. The values are list of
            arrays of shape n_points_per_slicex3, where the columns indicate
            the x, y and z coordinates of the contour point w.r.t. the
            Reference Coordinate System. The unit for all coordinates is mm.
            Each element in the list contains the points of a single slice,
            i. e. all the points in one array share the same z-coordinate.
        '''
        return copy.deepcopy(self._pointsets)

    def structure_set(self, ct):
        '''
        Return the structure set.

        Parameters
        ----------
        ct: CT
            The CT on which the binary masks of the structures will be
            calculated.
        '''
        structures = []

        for name, points in self.point_sets.items():
            if len(points) > 0:
                points = np.vstack(points)
            if self._use_cm:
                points = points / 10
            structure = Structure(name, points, ct.grid)
            structures.append(structure)

        return StructureSet(structures, ct)


class DicomDose:
    '''
    Class representing a DICOM RT Dose.
    '''

    def __init__(self, file, use_cm: bool = True):
        '''
        Parameters
        ----------
        file: str
            Path to a DICOM file containing the dose.
        use_cm: bool
            Whether to use cm as the unit for the contour point coordinates.
            Otherwise, mm is used.
            The unit must match the unit for the CT origin and spacing.
            Default: True
        '''
        self._use_cm = use_cm
        ds = pydicom.read_file(file)
        d = dicom_to_dict(ds, extract_pixel_array=True, recursive=True)

        self._matrix = d['Pixel Data'] * float(d['Dose Grid Scaling'])
        # column-major zyx -> row major xyz
        self._matrix = np.swapaxes(self._matrix, 0, 2)
        self._units = d['Dose Units']
        self._summation_type = d['Dose Summation Type']

        # Make sure the dose layers are equidistant along the z-axis.
        offsets = d['Grid Frame Offset Vector']
        # We have one difference less than slices
        differences = np.empty(len(offsets) - 1)

        for i in range(1, len(offsets)):
            differences[i-1] = offsets[i] - offsets[i-1]

        assert np.isclose(differences, differences[0]).all()

        self._grid_spacing = np.array([*d['Pixel Spacing'], differences[0]])

        self._image_position = np.array(d['Image Position (Patient)'])
        self._orientation = np.array(d['Image Orientation (Patient)'])
        # We support only HFS (head first suppine) for now.
        assert np.isclose(
            self._orientation,
            np.array([1, 0, 0, 0, 1, 0]),
            atol=1e-6,
        ).all()

        if self._use_cm:
            self._image_position /= 10.
            self._grid_spacing /= 10.

    @property
    def voxel_matrix(self):
        '''
        Returns
        -------
        m: array
            Shape: (n_pixels_x, n_pixels_y, n_pixels_z)
        '''
        return np.copy(self._matrix)

    @property
    def row_spacing(self):
        '''Spacing between rows.'''
        return self._grid_spacing[0]

    @property
    def col_spacing(self):
        '''Spacing between columns.'''
        return self._grid_spacing[1]

    @property
    def slice_spacing(self):
        '''Spacing between slices'''
        return self._grid_spacing[2]

    @property
    def image_position(self):
        '''
        Returns the position of the top left corner of the dose.

        The position of the middle point of the top left pixel is returned.
        '''
        return self._image_position

    @property
    def units(self):
        return self._units

    @property
    def summation_type(self):
        return self._summation_type

    @property
    def dose(self):
        return Dose(
            self.voxel_matrix,
            np.array([self.row_spacing, self.col_spacing, self.slice_spacing]),
            self.image_position,
        )
