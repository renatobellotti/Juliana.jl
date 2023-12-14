import logging
from pathlib import Path
import re
import struct
from typing import Dict, List, Union
import numpy as np
import pandas as pd
from pyftpp import CT, Dose, Grid, Structure, StructureSet


def read_three_longs(file):
    triple = []
    for i in range(3):
        triple.append(struct.unpack('<l', file.read(4))[0])
    return tuple(triple)


def parse_psiplan_ct_or_dose(path: str, dtype='int16'):
    with open(path, 'rb') as file:
        # I don't really understand why it's only 4 bytes
        # if we use little-endian...
        n_bytes = struct.unpack('<l', file.read(4))[0]

        size = read_three_longs(file)
        n_voxels = read_three_longs(file)
        # What is the difference between size and n_voxels???
        assert size == n_voxels
        # Index of bottom lower left corner.
        origin_index = read_three_longs(file)
        assert origin_index == (0, 0, 0)
        # Origin in real-world coordinates in header units (below).
        origin = read_three_longs(file)
        # Spacing in real-world coordinates in header units (below).
        spacing = read_three_longs(file)
        # ???
        data_type = struct.unpack('<l', file.read(4))[0]
        n_bits_per_voxel = struct.unpack('<l', file.read(4))[0]

        header_units = struct.unpack('<l', file.read(4))[0]
        data_units = struct.unpack('<l', file.read(4))[0]

        # ???
        n_files = struct.unpack('<l', file.read(4))[0]
        file_num = struct.unpack('<l', file.read(4))[0]
        assert n_files == 1
        assert file_num == 0
        couch = struct.unpack('<l', file.read(4))[0]

        # Unused values
        for i in range(4):
            _ = struct.unpack('<l', file.read(4))[0]

        patient_weight = struct.unpack('<l', file.read(4))[0]
        # "Additional weight of equipment on patient couch"
        extra_weight = struct.unpack('<l', file.read(4))[0]
        # "Total weight on couch"
        system_weight = struct.unpack('<l', file.read(4))[0]

        ct_type = struct.unpack('<l', file.read(4))[0]
        gantry_type = struct.unpack('<l', file.read(4))[0]

        study_id = struct.unpack('<l', file.read(4))[0]
        series_id = struct.unpack('<l', file.read(4))[0]

        couch_offset_y = struct.unpack('<l', file.read(4))[0]

        for _ in range(116):
            blank = struct.unpack('c', file.read(1))

        comment = ''.join([struct.unpack('c', file.read(1))[0].decode('ascii') for _ in range(256)])

        if n_bits_per_voxel != 16:
            logging.warning('CTCTFile ignores the data_size field in the CTCT file and assumes 2 bytes per voxel.')

        slices = []
        for z in range(size[2]):
            plane = np.frombuffer(
                # Each pixel in the slice needs two bytes.
                file.read(size[0]*size[1]*2),
                dtype=dtype,
            ).reshape(size[:2])
            slices.append(plane)
        data = np.stack(slices, axis=2)
        data = np.flip(data, axis=0)

    return origin, spacing, data, header_units, data_units


class CTCTFile:

    def __init__(self, path: str):
        self._path = path
        self._parse_file()

    def _parse_file(self):
        origin, spacing, data, header_units, data_units = parse_psiplan_ct_or_dose(self._path)

        self._data = data - 1000
        self._origin = np.array(origin)
        self._size = data.shape
        # The spacing is in [mm].
        self._spacing = np.array(spacing) / 10**(header_units-2)
    
    @property
    def ct(self) -> CT:
        '''
        Returns the CT. Unit of the spacing: cm (needed by Fiona standalone).
        '''
        if not np.isclose(np.array(self._origin), 0.).all():
            logging.warning('Flipping the CT with non-zero origin is not tested. Please check whether CT and contours are correctly displayed!!!')
        # Flipping the z-axis is needed so that the beam angles are correct
        # in Fiona GUI and Fiona standalone.
        return CT(self._data[:, :, ::-1], self._spacing, np.array(self._origin))


class PTDTFile:

    def __init__(self, path: str):
        self._path = path
        self._parse_file()

    def _parse_file(self):
        origin, spacing, data, header_units, data_units = parse_psiplan_ct_or_dose(self._path)

        # Division by 10 yields percentage points.
        self._data = data.astype('float32') / 10**(data_units) / 10.
        self._origin = np.array(origin)
        self._size = data.shape
        # The spacing is in [mm].
        self._spacing = np.array(spacing) / 10**(header_units-2)

    def dose(self, target_dose: Union[float, None] = None) -> Dose:
        '''
        Returns the dose distribution.
        Unit of the spacing: cm (needed by Fiona standalone).
        Unit of the dose: Percentage if target_dose is None (default), else Gy.

        Parameters
        ----------
        target_dose: float or None
            The dose corresponding to 100%.

        Returns
        -------
        Dose
        '''
        # Flipping the z-axis is needed so that the beam angles are correct
        # in Fiona GUI and Fiona standalone.
        if not np.isclose(np.array(self._origin), 0.).all():
            logging.warning('Flipping a dose distribution with non-zero origin is not tested. Please check whether CT and contours are correctly displayed!!!')
        data = self._data[:, :, ::-1]
        if target_dose is not None:
            data = data / 100. * target_dose
        return Dose(data, self._spacing, self._origin)


class VDXFile:

    def __init__(self, path: str, grid: Union[Grid, None] = None):
        self._grid = grid
        with open(path, 'r') as file:
            lines = file.readlines()
            lines = [l.strip() for l in lines]

        # Remove double empty lines to facilitate parsing.
        cleaned_lines = []
        for i, line in enumerate(lines):
            if (i < len(lines)-1) and (line == '') and (lines[i+1] == ''):
                continue
            cleaned_lines.append(line)
        del lines
        lines = cleaned_lines

        paragraph_start = 0
        paragraph_end = find_next_paragraph_end(lines)
        header_lines = lines[:paragraph_end]

        n_vois = parse_file_header(header_lines)

        info = []

        for _ in range(n_vois):
            paragraph_start = paragraph_end + 1
            voi_name, voi_key, voi_type, points, paragraph_end = parse_voi_lines(lines[paragraph_start:])
            paragraph_end += paragraph_start
            info.append({
                'name': voi_name.replace('/', '_'),
                'key': voi_key,
                'type': voi_type,
                'points': points,
            })

        self._info = info

    @property
    def info(self):
        return self._info

    def structureset(self, ct: CT) -> StructureSet:
        assert np.isclose(ct.origin, 0., atol=1e-6).all()
        structures = []
        for structure_info in self.info:
            name = structure_info['name']
            points = sorted(
                structure_info['points'],
                key=lambda p: p[0, 2],
            )
            # Drop duplicate rows. The DicomConverter seems to do it, too.
            points = np.vstack(points)
            # Taken from https://stackoverflow.com/a/15637512
            _, indices = np.unique(points, axis=0, return_index=True)
            points = points[np.sort(indices)]
            # Shift of half a slice spacing is needed because the VDX file
            # contains middler points of the slice (I think...)
            # The shift along the y-direction is also done in the DicomConverter
            # and is needed to align the target structure and the PSIPlan
            # dose distribution... Perhaps that was a PSIPlan bug...?
            points += np.array([0, -ct.spacing[1], -0.5*ct.spacing[2]])
            structure = Structure(name, points, self._grid)
            # Flipping the z-axis is needed so that the beam angles are correct
            # in Fiona GUI and Fiona standalone.
            structure = invert_z_for_structure(structure, ct)
            structures.append(structure)
        return StructureSet(structures, ct)


def find_next_paragraph_end(lines):
    '''
    Get all lines until the next empty line appears.
    '''
    for i, line in enumerate(lines):
        if line == '':
            return i

    return len(lines)



def parse_file_header(header_lines):
    assert len(header_lines) == 3
    assert header_lines[0] == 'vdx_file_version   2.0'
    assert header_lines[1] == 'all_indices_zero_based'
    assert header_lines[2].startswith('number_of_vois')

    n_vois = int(header_lines[2].split()[1].strip())

    return n_vois


def parse_structure_header(structure_header_paragraph):
    assert len(structure_header_paragraph) == 3
    assert structure_header_paragraph[0].startswith('voi')
    assert structure_header_paragraph[1].startswith('key')
    assert structure_header_paragraph[2].startswith('type')

    # Structure name, key and type (OAR or target (?)).
    voi_name = structure_header_paragraph[0].split()[1].strip()
    voi_key = structure_header_paragraph[1].split()[1].strip()
    voi_type = structure_header_paragraph[2].split()[1].strip()

    return voi_name, voi_key, voi_type


def line_to_triple(line):
    # [1:] is needed because we want to remove the line "header".
    return np.array([float(x) for x in line.strip().split()[1:]])


def parse_voi_header(contour_header_paragraph):
    assert contour_header_paragraph[0] == 'contours'
    assert contour_header_paragraph[1] == 'reference_frame'
    assert contour_header_paragraph[2].startswith('origin')
    assert contour_header_paragraph[3].startswith('point_on_x_axis')
    assert contour_header_paragraph[4].startswith('point_on_y_axis')
    assert contour_header_paragraph[5].startswith('point_on_z_axis')
    assert contour_header_paragraph[6].startswith('number_of_slices')

    # We don't need this information for now, but we check it to ensure consistency.
    reference_frame_origin = line_to_triple(contour_header_paragraph[2])
    reference_frame_point_on_x = line_to_triple(contour_header_paragraph[3])
    reference_frame_point_on_y = line_to_triple(contour_header_paragraph[4])
    reference_frame_point_on_z = line_to_triple(contour_header_paragraph[5])
    matrix = np.vstack([
        reference_frame_point_on_x,
        reference_frame_point_on_y,
        reference_frame_point_on_z,
    ])

    assert np.isclose(reference_frame_origin, np.zeros(3)).all()
    assert np.isclose(matrix, np.eye(3)).all()
    

    n_slices = int(contour_header_paragraph[6].split()[1].strip())

    return n_slices


def parse_contour(lines):
    assert lines[0].startswith('contour')
    assert lines[1].startswith('internal')
    assert lines[2].startswith('number_of_points')
    n_points = int(lines[2].replace('number_of_points', ''))
    assert len(lines) == (n_points + 3)
    
    points = []
    for line in lines[3:]:
        parts = line.strip().split()
        x, y, z, a, b, c = (float(p) for p in parts)
        assert np.isclose(np.array([a, b, c]), np.zeros(3)).all()
        points.append([x, y, z])
    return np.array(points) / 10.


def parse_voi_slice(lines):
    assert lines[0].startswith('slice ')
    assert lines[1].startswith('slice_in_frame ')
    assert lines[2].startswith('thickness')
    assert lines[3].startswith('number_of_contours')
    
    n_contours = int(lines[3].replace('number_of_contours', ''))
    assert n_contours == 1
    
    return parse_contour(lines[4:])


def parse_voi_lines(lines):
    # Header of the structure.
    paragraph_start = 0
    paragraph_end = paragraph_start + find_next_paragraph_end(lines)
    structure_header_paragraph = lines[paragraph_start:paragraph_end]

    voi_name, voi_key, voi_type = parse_structure_header(structure_header_paragraph)

    # Meta-information about the contours.
    paragraph_start = paragraph_end + 1
    paragraph_end = paragraph_start + find_next_paragraph_end(lines[paragraph_start:])
    contour_header_paragraph = lines[paragraph_start:paragraph_end]

    n_slices = parse_voi_header(contour_header_paragraph)

    # Slices of the contour.
    if n_slices == 0:
        logging.warning(f'Trying to load VDX structure with zero slices: {voi_name}!')

    points = []
    for _ in range(n_slices):
        paragraph_start = paragraph_end + 1
        paragraph_end = paragraph_start + find_next_paragraph_end(lines[paragraph_start:])
        points.append(parse_voi_slice(lines[paragraph_start:paragraph_end]))

    return voi_name, voi_key, voi_type, points, paragraph_end


def invert_z_for_structure(structure: Structure, ct: CT) -> Structure:
    new_points = structure.points.copy()
    for i in range(new_points.shape[0]):
        z = new_points[i, 2]
        iz = np.round(z / ct.spacing[2])
        assert np.isclose(z - iz * ct.spacing[2], 0, atol=1e-6)
        iz_new = ct.data.shape[2] - iz - 1
        z_new = iz_new * ct.spacing[2]
        new_points[i, 2] = z_new
    return Structure(structure.name, new_points, ct.grid)


class TrialDef:

    def __init__(self, path):
        with open(path, 'r') as file:
            lines = file.readlines()
            lines = [l.strip() for l in lines]

        n_iterations, fields, constraints, objectives = self._parse_lines(lines)

        self._path = path
        self._n_iterations = n_iterations
        self._fields = fields
        self._constraints = constraints
        self._objectives = objectives

    @property
    def path(self) -> str:
        return self._path

    @property
    def n_iterations(self) -> int:
        return self._n_iterations

    @property
    def fields(self) -> List[str]:
        return self._fields

    @property
    def constraints(self) -> List[Dict[str, str]]:
        return pd.DataFrame(self._constraints)

    @property
    def objectives(self) -> List[Dict[str, str]]:
        return pd.DataFrame(self._objectives)

    @staticmethod
    def _parse_lines(lines):
        assert (lines[0] == 'iterations:') \
            or (lines[0] == 'ititerations:')
        assert lines[1] == '# Define number of iterations'

        n_iterations = int(lines[2])

        assert lines[3] == ''
        assert lines[4] == 'fields:'
        assert (lines[5] == '# Define fields to be included') \
            or (lines[5] == '# define fields to be included')

        fields = []
        pattern = re.compile(r'F\d')
        i = 6
        # Skip empty lines.
        while lines[i].strip() == '':
            i += 1

        while pattern.match(lines[i]):
            fields.append(lines[i])
            i += 1

        # Skip empty lines.
        while lines[i].strip() == '':
            i += 1

        assert lines[i] == 'constraints:'
        assert lines[i+1] == '#Define constraints'
        parts = lines[i+2].split()
        assert len(parts) == 4
        assert parts[0] == '#OAR'
        assert parts[1] == 'Wgt.'
        assert parts[2] == 'Dose'
        assert parts[3] == 'Vol.'

        constraints = []
        i = i + 3
        while lines[i] not in  ['', 'prescriptions:']:
            parts = lines[i].split()
            assert len(parts) == 4
            constraints.append({
                'oar': parts[0],
                'weight': parts[1],
                'dose': parts[2],
                'volume': parts[3],
            })
            i += 1

        while lines[i] == '':
            i += 1
        assert lines[i] == 'prescriptions:'
        assert (lines[i+1] == '#Define ptv boost factors') \
            or (lines[i+1] == '#Define ptv dose prescription')
        parts = lines[i+2].split()
        assert len(parts) == 3
        assert parts[0] == '#PTV'
        assert parts[1] == 'Wgt.'
        assert parts[2] == 'Dose'

        objectives = []
        for line in lines[i+3:]:
            if line == '':
                continue

            parts = line.split()
            assert len(parts) == 3
            objectives.append({
                'target': parts[0],
                'weight': parts[1],
                'dose': parts[2],
            })

        return n_iterations, fields, constraints, objectives


class G2FieldInfo:
    '''
    Example of a field info file::

               9
        AUTO_1
              23
               0
               0
               0
               0
               0
              180.000
              180.000
            -0.240000     0.240000     -135.910
             -45.0000
             0.250000     0.400000     0.400000
             -60      60
             -40      40
             -40      40
               1       1       1
             0.500000
             0.200000
             0.800000
             0.800000
             0.800000
              0.00000
             0.900000
              10.0000
               0
               0
             -60      60
             -40      40
             -40      40
               1       1       1
               1       1       0       0       0       0       0       0       0
               0

    Documentation for FIELD_INFO2 file format:
    https://codebeamer.psi.ch/cb/doc/292094
    '''

    def __init__(self, path: str):
        '''
        Parameters
        ----------
        path: str
            path to a field info file
        '''
        tmp = Path(path).stem
        self._patient_ID = tmp.split('_')[0]
        self._ID = tmp.replace(f'{self._patient_ID}_', '')
        self._name = self._ID.split('_')[-1]

        with open(path, 'r') as file:
            lines = file.readlines()
            self._target_structure_key = lines[0].strip()
            self._preabsorber = lines[1].strip()
            self._nozzle_extraction = float(lines[2].strip())
            self._gantry_angle = float(lines[8].strip())
            self._couch_angle = float(lines[9].strip())

    @property
    def patient_ID(self) -> str:
        return self._patient_ID

    @property
    def ID(self) -> str:
        return self._ID

    @property
    def name(self) -> str:
        return self._name

    @property
    def preabsorber(self) -> str:
        '''
        WARNING: Untested code!!!
        '''
        # TODO: Write tests.
        return self._preabsorber

    @property
    def nozzle_extraction(self) -> float:
        '''
        Nozzle extraction in cm (according to PSIPlan).
        '''
        return self._nozzle_extraction

    @property
    def gantry_angle(self) -> float:
        '''
        Gantry angle in degrees.
        '''
        return self._gantry_angle

    @property
    def couch_angle(self) -> float:
        '''
        Couch angle in degrees.
        '''
        return self._couch_angle

    @property
    def target_structure_key(self) -> str:
        '''
        Key of the target structure in the corresponding VDX_PSI file.
        '''
        return self._target_structure_key


class G2PlanDefinition:
    '''
    Example of a plan definition file::

        \t  1.63636
        \t   3
        \t  1.00000      1.00000     0.333333     0.545455 F0
        \t  1.00000      1.00000     0.333333     0.545455 F1
        \t  1.00000      1.00000     0.333333     0.545455 F2
        \t  <blank line>
    '''

    def __init__(self, path: str):
        '''
        Parameters
        ----------
        path: str
            path to a plan definition file
        '''
        # Example path:
        # /asm/p2data/PRD/<patient_ID>/<patient_ID>_CT<i>_T<j>_P<k>.PLAN_DEF
        tmp = Path(path).stem
        self._patient_ID = tmp.split('_')[0]
        self._ID = tmp.replace(f'{self._patient_ID}_', '')
        self._directory = str(Path(path).parent)

        with open(path, 'r') as file:
            lines = file.readlines()

            self._fraction_dose = float(lines[0].strip())
            self._n_fields = int(lines[1].strip())

            self._field_names = []
            self._field_weights = {}

            for i in range(self._n_fields):
                columns = lines[2 + i].split()
                name = columns[-1].strip()
                weight = float(columns[2].strip())
                self._field_names.append(name)
                self._field_weights[name] = weight

    @property
    def patient_ID(self) -> str:
        return self._patient_ID

    @property
    def ID(self) -> str:
        return self._ID

    @property
    def n_fields(self) -> int:
        return self._n_fields

    @property
    def field_names(self) -> List[str]:
        return self._field_names

    @property
    def fields(self) -> List[G2FieldInfo]:
        tmp = self.ID.split('_')
        ct = tmp[0][2:]
        t = tmp[1][1:]

        fields = []
        for f in self.field_names:
            filename = f'{self.patient_ID}_CT{ct}_T{t}_{f}.FIELD_INFO2'
            path = f'{self._directory}/{filename}'
            fields.append(G2FieldInfo(path))
        return fields

    @property
    def field_weights(self) -> Dict[str, float]:
        '''
        The dict keys are the field names, e. g. "F3".
        '''
        return self._field_weights

    @property
    def dose_per_frac(self) -> float:
        '''
        Returns the dose_per_fraction * 1.1.
        '''
        return self._fraction_dose * 1.1

    def __str__(self):
        return str(self.field_names)


class G2TreatmentDefinition:
    '''
    Class representing a treatment definition.

    Example of a treatment definition file::

        TREATMENT_NUMBER: 0
        NUMBER_OF_SERIES: 2;
        CT_SET: 0;
        T_SET: 3;
        PLAN: 0;
        START_DOSE: 0.00000;
        END_DOSE: 54.0000;
        END_OF_SERIES:
        CT_SET: 0;
        T_SET: 4;
        PLAN: 0;
        START_DOSE: 54.0000;
        END_DOSE: 74.0000;
        END_OF_SERIES:
    '''

    def __init__(self, path: str):
        '''
        Parameters
        ----------
        path: str
            path to a treatment definition file
        '''
        self._ct_numbers = []
        self._t_numbers = []
        self._plans = []
        self._start_doses = []
        self._end_doses = []
        self._patient_ID = Path(path).stem.split('_')[0]
        self._patient_directory = str(Path(path).parent)

        with open(path, 'r') as file:
            lines = [line.strip() for line in file.readlines()]
            line = lines[0]
            if not line.startswith('TREATMENT_NUMBER: '):
                raise RuntimeError(f'Malformatted treatment definition file: {path}')
            self._treatment_number = int(line.split(' ')[1])

            line = lines[1]
            if not line.startswith('NUMBER_OF_SERIES: '):
                raise RuntimeError(f'Malformatted treatment definition file: {path}')
            self._n_series = int(line.split(' ')[1].replace(';', ''))

            # Process the series information.
            error = RuntimeError(f'Malformatted treatment definition file: {path}')

            lines = lines[2:]
            i = 0
            while (i < len(lines)) and (lines[i] != ''):
                line = lines[i]
                if not line.startswith('CT_SET: '):
                    raise error
                self._ct_numbers.append(int(line.split(' ')[1].replace(';', '')))

                line = lines[i + 1]
                if not line.startswith('T_SET: '):
                    raise error
                self._t_numbers.append(int(line.split(' ')[1].replace(';', '')))

                line = lines[i + 2]
                if not line.startswith('PLAN: '):
                    raise error
                self._plans.append(int(line.split(' ')[1].replace(';', '')))

                line = lines[i + 3]
                if not line.startswith('START_DOSE: '):
                    raise error
                self._start_doses.append(float(line.split(' ')[1].replace(';', '')))

                line = lines[i + 4]
                if not line.startswith('END_DOSE: '):
                    raise error
                self._end_doses.append(float(line.split(' ')[1].replace(';', '')))

                line = lines[i + 5]
                if line != 'END_OF_SERIES:':
                    raise error

                i += 6

    def ct_set(self, series: int) -> int:
        return self._ct_numbers[series]

    def t_set(self, series: int) -> int:
        return self._t_numbers[series]

    def plan_ID(self, series: int) -> int:
        return self._plans[series]

    def plan(self, series: int) -> G2PlanDefinition:
        # TODO: Write tests for this.
        ct = self.ct_set(series)
        t = self.t_set(series)
        plan = self.plan_ID(series)
        filename = f'{self._patient_ID}_CT{ct}_T{t}_P{plan}.PLAN_DEF'
        path = f'{self._patient_directory}/{filename}'
        return G2PlanDefinition(path)

    def start_dose(self, series: int) -> float:
        return self._start_doses[series]

    def end_dose(self, series: int) -> float:
        return self._end_doses[series]

    @property
    def n_series(self) -> int:
        # Could also choose any other member list, they have the same length.
        return len(self._ct_numbers)

    def __str__(self) -> str:
        msg = []
        for i in range(len(self._ct_numbers)):
            msg.append('\n'.join([
                '----------------',
                f'Series {i}',
                '----------------',
                f' CT Set: {self.ct_set(i)}',
                f'  T Set: {self.t_set(i)}',
                f'   Plan: {self.plan(i)}',
                f'D_start: {self.start_dose(i)}',
                f'  D_end: {self.end_dose(i)}',
            ]))
        return '\n'.join(msg)
