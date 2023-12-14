import os
from pathlib import Path
from typing import List
import numpy as np
from .config import PyFTPPConfigurator, PyFTPP_slice
from .ct import CT
from .data import DataSet
from .grid import Grid
from .structure import load_array_structure


class ConfigurationFactory:
    '''
    This class aims to provide an easy-to-use interface for loading all
    relevant patient information in order to avoid code duplication.
    '''

    def __init__(self,
                 patient_data: DataSet,
                 ct_path: str,
                 target_structure: str,
                 target_dose: float,
                 target_dose_at_edge: float,
                 dose_grid_resolution: float,
                 optimisation_grid_resolution: float = 0.35):
        '''
        Parameters
        ----------
        patient_data: DataSet
            All the data about the patient.
        ct_path: str
            Path to the CT file in the Fiona standalone format.
        target_structure: str
            Name of the target structure.
            This is where the spots will be placed.
        target_dose: float
            Target dose in Gy.
        target_dose_at_edge: float
            Target dose at edge in Gy.
        dose_grid_resolution: float
            Resolution of the dose calculation grid in cm.
        optimisation_grid_resolution: float
            Resolution of the optimisation grid in cm.
        '''
        self._data = patient_data
        self._ct_path = ct_path
        self._target_structure = target_structure
        self._bin_dir = Path(f'{__file__}/../bin')

        # We need to get a boundary box for all structures
        # to know where to put the optimisation grid.
        self._p_min_x = np.inf
        self._p_min_y = np.inf
        self._p_min_z = np.inf
        self._p_max_x = -np.inf
        self._p_max_y = -np.inf
        self._p_max_z = -np.inf

        # Configuration parameters.
        self._target_dose = target_dose
        self._target_dose_at_edge = target_dose_at_edge
        self._dose_grid_resolution = dose_grid_resolution
        self._optim_grid_resolution = optimisation_grid_resolution

        self._constraints = []
        self._fields = []

    def _find_ROI(self):
        '''
        Finds the region of interest (ROI), i. e. a bounding box around the
        target structure and all OARs with at least one constraint.
        The bounding box will confine the spot optimisation grid.
        '''
        # This implementation is an attempt at copying Fiona's behaviour
        # in order to have the same (or at least very similar) results as if
        # one were using the GUI.
        # The optimisation grid plays an important role because
        # Siddon's algorithm is very sensitive w. r. t. even sub-voxel shifts,
        # which leads to jumps in the calculated dose distribution in regions
        # where the stopping power changes abruptly, e. g. at the edges of
        # bones.
        # The relevant code in the Fiona repository can be found at:
        # shared/ch.psi.ftpp.planning.impl.optimization/src/ch/psi/ftpp/planning/impl/optimization/OptimizationGridProvider.java

        # Get list of relevant structures.
        structure_names_for_bounding_box = [self._target_structure]
        for constraint in self._constraints:
            structure_names_for_bounding_box.append(constraint['structure_label'])

        # Calculate the extreme points that define the bounding box around
        # these structures.
        p_min = np.array([np.inf, np.inf, np.inf])
        p_max = np.array([-np.inf, -np.inf, -np.inf])
        for structure_name in structure_names_for_bounding_box:
            structure = self._data.structures[structure_name]
            structure_min = np.min(structure.points, axis=0)
            structure_max = np.max(structure.points, axis=0)
            p_min = np.minimum(p_min, structure_min)
            p_max = np.maximum(p_max, structure_max)

        # Extend the bounding box.
        EXTENSION = 3.
        p_min = p_min - EXTENSION
        p_max = p_max + EXTENSION

        # Round the origin in order to make them lie on CT grid nodes.
        origin = self._data.ct.origin
        spacing = self._data.ct.spacing
        p_min = np.floor((p_min - origin) / spacing) * spacing + origin

        # Needed in Fiona because of the way in which voxels are assigned to
        # structures? I don't understand it fully...
        p_min[2] += 0.1

        self._p_min_x = p_min[0]
        self._p_min_y = p_min[1]
        self._p_min_z = p_min[2]
        self._p_max_x = p_max[0]
        self._p_max_y = p_max[1]
        self._p_max_z = p_max[2]

    @property
    def configurator(self) -> PyFTPPConfigurator:
        self._find_ROI()

        # Configure the target and its prescription.
        conf = PyFTPPConfigurator(self._ct_path,
                                  self._target_dose,
                                  self._target_dose_at_edge)
        conf.set_target(0, self._data.structures[self._target_structure])

        # Set optimisation grid.
        d = tuple([
            self._p_max_x - self._p_min_x,
            self._p_max_y - self._p_min_y,
            self._p_max_z - self._p_min_z,
        ])
        optim_grid_origin = tuple((
            self._p_min_x,
            self._p_min_y,
            self._p_min_z,
        ))
        optim_grid_spacing = [
            self._optim_grid_resolution,
            self._optim_grid_resolution,
            self._optim_grid_resolution,
        ]
        optim_grid_size = tuple([np.ceil(d[i] / optim_grid_spacing[i]) for i in range(3)])
        optim_grid = Grid(
            optim_grid_size, optim_grid_spacing, optim_grid_origin
        )
        conf.set_optimization_grid(optim_grid)

        # Set the dose calculation grid.
        conf.set_dose_resolution(self._dose_grid_resolution)

        # Add the constraints.
        # We order the structures with constraints alphabetically to ensure
        # reproducibility: The same structure should have the same index
        # over subsequent runs. Sets eliminate duplicates, but are not ordered.
        structure_names = set([c['structure_label'] for c in self._constraints])
        structure_names = sorted(list(structure_names))
        for i, constr in enumerate(self._constraints):
            label = constr['structure_label']
            conf.add_constraint(structure_ID=structure_names.index(label),
                                constraint_ID=i,
                                structure=self._data.structures[label],
                                importance=constr['importance'],
                                dose=constr['dose'],
                                volume=constr['volume'],
                                type=constr['type'],
                                enforced=True,
                                limitType=constr['limitType'])
        # Add the fields.
        for field in self._fields:
            conf.add_field(**field)

        return conf

    def add_constraint(self,
                       structure_label: str,
                       dose: float,
                       volume: float,
                       type: str,
                       importance: float,
                       limitType: str = 'UPPER'):
        '''
        Parameters
        ----------
        structure_label: str
            Name of the structure. Must be consistent with the naming of the
            structures file.
        dose: float
            Dose threshold for mean/max dose constraints. Also used for DVH
            constraints:
                dose = 96, volume = 2 means D2 < 96% * target_dose
                ("hottest 2% of the volume receive no more than 96% of the target dose")
                (https://www.carlosjanderson.com/understanding-the-meaning-of-dvh-metrics/)
                (viewed 2021-09-10)
                ???
        volume: float
            Only used for DVH constriants. See documentation for the dose
            parameter.
        '''
        self._constraints.append({
            'structure_label': structure_label,
            'dose': dose,
            'volume': volume,
            'type': type,
            'importance': importance,
            'limitType': limitType,
        })

    def add_field(self,
                  gantry_angle: float,
                  couch_angle: float,
                  nozzle_extraction: float,
                  spot_spacing: float,
                  preabsorber_setting: str = 'auto',
                  margin: float = 0.5,
                  target_ID: int = 0):
        if preabsorber_setting.lower().startswith('preab'):
            preabsorber_setting = 'IN'
        elif preabsorber_setting.lower().startswith('auto'):
            preabsorber_setting = 'AUTO'
        elif preabsorber_setting.lower().startswith('nopreab'):
            preabsorber_setting = 'OUT'
        else:
            raise RuntimeError(f'Unexpected preabsorber setting {preabsorber_setting}')

        self._fields.append({
            'gantry_angle': gantry_angle,
            'couch_angle': couch_angle,
            'nozzle_extraction': nozzle_extraction,
            'margin': margin,
            'preabsorberSetting': preabsorber_setting,
            'xSpacing': spot_spacing,
            'ySpacing': spot_spacing,
            'target_ID': target_ID,
        })
