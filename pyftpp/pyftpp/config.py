from collections import defaultdict
import json
import logging
from typing import List, Tuple, Dict
import numpy as np
from dataclasses import dataclass, field, asdict
from .ct import CT
from .dose import Dose
from .grid import Grid
from .structure import Structure, PyFTPP_slice


# The contour in list-of-slices format, i. e. every slice is represented
# as a dictionary with a single key "points" whose value is a list of dict
# with keys "x", "y", "z". The z-coordinate is the same for all points
# in a slice.
# Example:
# {
#     'points': [
#         {'x': 1, 'y': 2, 'z': 3},
#         {'x': 4, 'y': 5, 'z': 3},
#         ...
#     ]
# }


class PyFTPPConfigurator:

    def __init__(self,
                 ct_path: str,
                 target_dose: float,
                 target_dose_at_edge: float):
        self._ct_path = ct_path
        self._target_dose = target_dose
        self._target_dose_at_edge = target_dose_at_edge
        self._target = None
        self._fields: List[FieldDefinition] = []
        # The following dictionaries have the OAR ID as keys.
        self._oar_labels: Dict[int, str] = {}
        self._constraints: Dict[int, List[Constraint]] = defaultdict(list)
        self._oar_slices: Dict[int, List[Dict]] = {}
        self._optimization_grid: OptimizationGrid = None
        self._dose_resolution: float = None
        self._n_iterations: int = 100

    def set_target(self, ID: int, structure: Structure):
        '''
        Set the target structure for the spot optimization.

        Parameters
        ----------
        ID: int
            ID of the target structure.
        structure: Structure
            The structure to use as the target.
        '''
        self._target = StructureDefinition(structure.name, ID, structure.fiona_slices)

    def add_field(self,
                  target_ID: int,
                  gantry_angle: int,
                  couch_angle: int,
                  nozzle_extraction: int = 15,
                  preabsorberSetting: str = 'OUT',
                  margin: float = 0.5,
                  xSpacing: float = 0.4,
                  ySpacing: float = 0.4):
        '''
        Add a field to this plan configuration.

        Parameters
        ----------
        target_ID: int
            The centre of the structure corresponding to the target_ID is used
            as the isocentre of the field.
        gantry_angle: int
            Gantry angle in degrees. Range: [-30, 180] for G2.
        couch_angle: int
            Couch angle in degrees. Range: [0, 180] for G2.
        nozzle_extraction: int
            Nozzle extraction in cm.
        preabsorber_setting: str
            Can be "IN", "OUT", or "AUTO".
        margin: float
            The spots or placed within the interval
            [minimum_coord-margin, maximum_coord+margin] for voxels inside the
            target. Unit: cm
        xSpacing: float
            Spacing of the spots in x direction. Unit: cm
        ySpacing: float
            Spacing of the spots in y direction. Unit: cm
        '''
        self._fields.append(FieldDefinition(target_ID,
                                            gantry_angle,
                                            couch_angle,
                                            nozzle_extraction,
                                            preabsorberSetting,
                                            margin,
                                            xSpacing=xSpacing,
                                            ySpacing=ySpacing))

    def add_constraint(self,
                       structure_ID: int,
                       constraint_ID: int,
                       structure: Structure,
                       importance: float,
                       dose: float,
                       volume: float,
                       type: str,
                       enforced: bool = True,
                       limitType: str = 'UPPER'):
        '''
        Add a constraint for the spot optimiser.

        Parameters
        ----------
        structure_label: str
            Name of the structure to constrain.
        structure_ID: int
            ID of the structure to constrain.
        constraint_ID: int
            ID of the constraint being added.
        structure: Structure
            The structure to constrain.
        importance: float
            Importance of the constraint. An importance of 1 means that the
            constraint is as important as covering the target, zero means
            that the constraint is deactivated.
        dose: float
            Dose value of the constraint. A value of 1 means
            100% of the target dose, a value of 1.2 means
            120% of the target dose.
            This means that dose% of the target dose should be limited. E. g.:
            D_volume < dose% * target_dose
        volume: float
            Volume for the DVH to be constraint. This number is in [0, 1].
            A value of zero means D0, i. e. the maximum dose.
        type: str
            Can be "DOSE_VOLUME" or "MEAN".
            If it is "MEAN", then the volume parameter is not used.
        enforced: bool
            Whether to enforce the constraint.
            TODO: Check if this is even needed...Perhaps it is related to the
            optimiser by Michael Matter?
        limitType: str
            Either "UPPER" or "LOWER".
            If "UPPER": D_volume < dose% * target_dose
            If "LOWER": D_volume > dose% * target_dose
            Lower limits can be useful to boost some target subvolumes.
        '''
        self._oar_labels[structure_ID] = structure.name
        self._oar_slices[structure_ID] = structure.fiona_slices
        self._constraints[structure_ID].append(Constraint(
            constraint_ID,
            importance,
            dose,
            volume,
            type,
            enforced,
            limitType
        ))

    def set_optimization_grid(self, grid: Grid):
        '''
        Define the grid for the optimisation.

        By default (i. e. if this function is not called), the CT grid is used.

        Parameters
        ----------
        grid: Grid
            The new optimisation grid.
        '''
        self._optimization_grid = OptimizationGrid(grid.shape, grid.spacing, grid.origin)

    def set_dose_resolution(self, res: float):
        '''
        Set the resolution of the dose calculation grid.

        Parameters
        ----------
        res: float
            Spacing of the dose calculation grid in cm.
        '''
        self._dose_resolution = res

    def set_n_iterations(self, n):
        '''Set the number of iterations for the spot optimiser.'''
        self._n_iterations = n

    def write_config(self, output_dir: str, path_to_bin: str):
        '''
        Write the configuration files for the standalone executable.

        Parameters
        ----------
        output_dir: str
            Directory that will contain all config files and the optimized plan
            and calculated dose distribution.
        path_to_bin: str
            Directory where the standalone executable and its static config
            files are stored.
        '''
        if not self._fields:
            logging.warning('No fields have been defined!')
        if not self._target:
            logging.warning('No target has been defined!')

        spot_config_path = f'{output_dir}/spot_config.json'
        optim_config_path = f'{output_dir}/optim_config.json'
        result_dose_path = f'{output_dir}/result_dose.dat'
        result_plan_path = f'{output_dir}/result_plan.json'
        properties_path = f'{output_dir}/config.properties'

        # Build the constraints in the correct format.
        structure_constraints = []
        for structure_ID, list_of_constr in self._constraints.items():
            structure_constraints.append(StructureConstraints(
                StructureDefinition(
                    self._oar_labels[structure_ID],
                    structure_ID,
                    self._oar_slices[structure_ID]
                ), list_of_constr
            ))

        # Write the spot config JSON file.
        spot_config = SpotConfig(self._fields, [self._target])
        with open(spot_config_path, 'w') as file:
            json.dump(asdict(spot_config), file, indent=4)

        # Write the optimisation config JSON file.
        if self._optimization_grid is None:
            ct = CT.load(self._ct_path)
            self._optimization_grid = OptimizationGrid(
                ct.data.shape,
                (
                    float(ct.spacing[0]),
                    float(ct.spacing[1]),
                    float(ct.spacing[2]),
                ),
                (
                    float(ct.origin[0]),
                    float(ct.origin[1]),
                    float(ct.origin[2]),
                ),
            )

        optim_config = OptimizationSettings(self._target_dose,
                                            self._target_dose_at_edge,
                                            self._target,
                                            structure_constraints,
                                            self._optimization_grid,
                                            numberOfIterations=self._n_iterations)
        with open(optim_config_path, 'w') as file:
            json.dump(asdict(optim_config), file, indent=4)

        # Write the properties file.
        config = PyFtppConfig(self._ct_path,
                              optim_config_path,
                              spot_config_path,
                              result_dose_path,
                              result_plan_path,
                              self._target_dose,
                              path_to_bin,
                              plan=result_plan_path)
        if self._dose_resolution is not None:
            config.doseResolution = self._dose_resolution
        with open(properties_path, 'w') as file:
            for key, val in asdict(config).items():
                if (val is None) or (key == 'binPath'):
                    continue
                file.write(f'{key}={val}\n')


@dataclass
class PyFtppConfig:
    '''
    Main config file for the standalone program.

    Contains paths to various config files and the results, as well as some
    default settings that don't need to be changed by the user.
    '''
    # Paths to data and config files.
    ct: str
    optimizationSettings: str
    spotConfiguration: str
    doseResults: str
    # Path to which the plan JSON file will be written. (output!)
    planResults: str
    prescribedDose: float
    binPath: str
    # Path to a plan JSON (input). When using the RUN command,
    # this is overwritten during thefield definition step. However, useful if
    # only a dose calculation is needed.
    plan: str = None
    preabsorberWed: float = 4.173
    cutoffFactor: float = 5.0
    doseResolution: float = 0.35
    huToSp: str = field(init=False)
    beamdata: str = field(init=False)
    beamlinePhasespace: str = field(init=False)
    nozzlePhasespace: str = field(init=False)
    absorberPhasespace: str = field(init=False)

    def __post_init__(self):
        self.huToSp = f'{self.binPath}/huToSp.json'
        self.beamdata = f'{self.binPath}/beamdata.xml'
        self.beamlinePhasespace = f'{self.binPath}/beamline.xml'
        self.nozzlePhasespace = f'{self.binPath}/nozzle.xml'
        self.absorberPhasespace = f'{self.binPath}/absorber.xml'


@dataclass
class StructureDefinition:
    '''
    Description of a structure (valid for both target and OAR).

    Some identifiers and the contour information.
    '''
    label: str
    id: int
    # The contour in list-of-slices format, i. e. every slice is represented
    # as a dictionary with a single key "points" whose value is a list of dict
    # with keys "x", "y", "z". The z-coordinate is the same for all points
    # in a slice.
    contours:  List[PyFTPP_slice]
    zMin: float = field(init=False)
    zMax: float = field(init=False)

    def __post_init__(self):
        z_coords = np.array([slice['points'][0]['z'] for slice in self.contours])
        self.zMin = np.min(z_coords)
        self.zMax = np.max(z_coords)


@dataclass
class Constraint:
    '''
    Constraint to a target or OAR.
    '''
    id: int
    importance: float
    dose: float
    volume: float
    # Type can be either 'DOSE_VOLUME' or 'MEAN'.
    type: str
    enforced: bool = True
    # limitType is either "UPPER" or "LOWER".
    limitType: str = 'UPPER'


@dataclass
class StructureConstraints:
    structure: StructureDefinition
    constraints: List[Constraint]


@dataclass(frozen=True)
class OptimizationGrid:
    size: Tuple[int]
    spacing: Tuple[float]
    origin: Tuple[float]


@dataclass
class OptimizationSettings:
    '''
    Representation of an optimization settings JSON.

    This config file contains hardcoded values as well as the contours and
    their associated constraints.

    For the target and every constraint, we have a dictionary "contours"
    that contains all points of the contours. It is implemented as a list of
    slices.
    Every slice is a dictionary with a single entry "points", which is a list
    of dicts with the keys "x", "y", and "z". The z-coordinate is the same for
    all points within the slice.

    Apart from the contours, we have also "zMin", "zMax", "label"
    (name of the structure), and "id", which is an arbitrarily chosen integer
    number (has to be the same as in the spot config file).
    '''
    targetDose: float
    targetDoseAtEdge: float
    target: StructureDefinition
    structureConstraints: List[StructureConstraints]
    optimizationGrid: OptimizationGrid
    numberOfIterations: int = 100
    sigmaFalloff: float = 0.8
    importanceExponent: float = 10.0
    falloffDistance: float = -0.1672349
    inTargetDistance: float = 3.0
    minWeight: float = 70.0
    hotSpotEnabled: bool = False
    hotSpotDose: float = 107.0
    hotSpotImportance: float = 3.0
    rrEnabled: bool = False
    rrScaleValue: float = 3.0
    rrScenarioImportance: float = 1.0


@dataclass
class FieldDefinition:
    targetStructureId: int
    gantryAngle: int
    couchAngle: int
    nozzleExtraction: int = 15
    preabsorberSetting: str = 'OUT'
    margin: float = 0.5
    xSpacing: float = 0.4
    ySpacing: float = 0.4


@dataclass
class SpotConfig:
    '''
    Contains the field specification and a list of targets.
    '''
    fields: List[FieldDefinition]
    # Contains a list of target structures.
    structures: List[StructureDefinition]
