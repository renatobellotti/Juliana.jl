import json
import struct
import numpy as np
import pandas as pd
from scipy.sparse import csc_array


class FionaDijHandler:

    def __init__(self, directory, grid):
        self._dir = directory
        self._grid = grid
        self._spots = None
        self._prescribed_dose = None
        self._fields = None

    @property
    def Dij(self):
        return load_dij_matrix(f'{self._dir}/dij_matrix.dat')

    @property
    def spots(self):
        if self._spots is None:
            self._spots = load_spots(f'{self._dir}/result_plan.json')
        return self._spots.copy()

    @property
    def optimisation_point_grid_indices(self):
        optim_point_file = f'{self._dir}/optimization_points.json'
        optimisation_points = load_optimisation_points(optim_point_file)
        # Map the indices of the optimisation points to flattened (ravelled)
        # CT grid indices.
        return points_to_indices(optimisation_points, self._grid)

    def write_spots(self):
        if self._fields is None:
            self._prescribed_dose, self._fields = load_fields(f'{self._dir}/result_plan.json')
        for field in self._fields:
            spots = self.spots.query('`field` == @field["id"]').copy()
            spots.drop(columns='field', inplace=True)
            s = spots.reset_index().to_json(orient='records')
            field['spots'] = json.loads(s)

        with open(f'{self._dir}/result_plan.json', 'w') as file:
            plan = {
                'prescribedDose': self._prescribed_dose,
                'fields': self._fields
            }
            json.dump(plan, file)

    def set_spot_weights(self, w):
        if self._spots is None:
            self._spots = load_spots(f'{self._dir}/result_plan.json')
        self._spots['weight'] = w


def load_dij_matrix(dij_matrix_path):
    with open(dij_matrix_path, 'rb') as file:
        n_entries = struct.unpack('<i', file.read(4))[0]
        n_spots = struct.unpack('<i', file.read(4))[0]
        # Each entry needs 4 bytes.
        entries = np.frombuffer(file.read(n_entries*4), dtype='float32')
        index = np.frombuffer(file.read(n_entries*4), dtype='int32')
        spot_stop = np.frombuffer(file.read((n_spots+1)*4), dtype='int32')
        return csc_array((entries, index, spot_stop))


def load_spots(plan_file):
    with open(plan_file, 'r') as file:
        data = json.load(file)

    spots = []
    for field in data['fields']:
        for spot in field['spots']:
            spot['field'] = field['id']
        spots.extend(field['spots'])
    spots = pd.DataFrame(spots)
    spots.set_index('id', inplace=True)
    assert (spots.index.value_counts() == 1).all()
    return spots


def load_fields(plan_file):
    with open(plan_file, 'r') as file:
        data = json.load(file)
    return data['prescribedDose'], data['fields']


def load_optimisation_points(path: str) -> pd.DataFrame:
    with open(path) as file:
        raw_data = json.load(file)
        data = {
            'optim_grid_index': raw_data['indices'],
            'x': raw_data['xCoordinates'],
            'y': raw_data['yCoordinates'],
            'z': raw_data['zCoordinates'],
        }
        df = pd.DataFrame(data)
        # Make sure the index column is unique.
        assert (df['optim_grid_index'].value_counts() == 1).all()
        df.set_index('optim_grid_index', inplace=True)
        return df


def mask_at_indices(indices: np.array, mask: np.array) -> np.array:
    '''
    Parameters
    ----------
    indices: np.array
        Array of shape (n, 3) of indices. The mask will be evaluated at these
        positions.
    mask: np.array
        The binary mask that should be evaluated at the given indices.
    Returns
    -------
    np.array
        Array of shape (n,) containing the values of the binary mask at the
        given indices.
    '''
    return np.array([mask[x, y, z] for x, y, z in indices])


def points_to_indices(points, grid) -> np.ndarray:
    return (points / grid.spacing).astype(int)
