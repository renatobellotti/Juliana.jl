'''
This module provides classes taking care of running the standalone spot
optimisation and dose calculation program.
'''
import json
import os
import shutil
from typing import Dict
from .config import PyFTPPConfigurator
from .plan import Plan
from .standalone import Standalone, get_path_to_bin, get_path_to_jar


class PyFTPPRunner:
    '''
    This class is a simple wrapper and scheduler for the standalone FTPP
    program.

    It allows to perform a spot optimisation and dose calculation as a function
    call.
    '''

    def __init__(self, output_dir: str, prune_spots: bool = True):
        '''
        Parameters
        ----------
        output_dir: str
            Output directory. All the config files and results will be written
            to this directory. Every configuration is stored in a separate
            subfolder.
        path_to_bin: str
            Path to the directory that contains the JAR and the static config
            files for the stanbdalone FTPP program.
        '''
        self._output_dir = output_dir
        self._configs: Dict[str, PyFTPPConfigurator] = {}
        self._path_to_bin = get_path_to_bin()
        self._prune_spots = prune_spots

    def add_configuration(self, config: PyFTPPConfigurator, label: str):
        '''
        Add another configuration to be run.

        The default implementation will just evaluate all configurations
        as a subprocrocess in a sequential way. More sophisticated
        implementations could realise running on clusters, e. g. using SLURM,
        and/or implementing load balancing, including resource handling
        (using CPUs for some configurations and GPUs for others to fully
        utilise the available resources).

        Parameters
        ----------
        config: PyFTPPConfigurator
            The plan configuration for which to perform a spot optimisation and
            calculate the dose distribution.
        label: str
            Label that is used to refer to this configuration. This will be
            used as the name of the subdirectory containing the configuration
            and the results of this evaluation, so it must be a valid directory
            name.
            The label should be unique.
        '''
        self._configs[label] = config

    def _run_config(self, label: str, debug: bool = True, log_dij_matrix: bool = False, log_wed: bool = False):
        '''
        Generates the config files and runs the configuration with the given
        label.

        Parameters
        ----------
        label: str
            Label of the configuration to run.
        debug: bool
            If True, only write the configuration files.
        log_dij_matrix: bool
            If True, write the Dij matrix and the optimisation points to the
            output directory.
        log_wed: bool
            If True, write the WED vectors to the output directory.

        Returns
        -------
        str
            The directory containing the configuration files and the results.
        '''
        workdir = f'{self._output_dir}/{label}'
        if not os.path.exists(workdir):
            os.makedirs(workdir)

        # Log the version of the standalone program.
        shutil.copyfile(
            f'{self._path_to_bin}/commit.txt',
            f'{workdir}/commit.txt'
        )

        # Write the config files.
        self._configs[label].write_config(workdir, self._path_to_bin)

        if debug:
            return workdir

        # Run the standalone program.
        property_file = f'{workdir}/config.properties'
        jar_file = get_path_to_jar()


        # We need to create this empty file so we don't get errors when
        # calling the standalone executable.
        plan_file = f'{workdir}/result_plan.json'
        with open(plan_file, 'w') as file:
            json.dump({}, file)
        
        dij_dir = workdir if log_dij_matrix else None
        wed_dir = workdir if log_wed else None

        # Place spots and optimise them, then remove unnecessary spots to
        # speed up the dose calculation.
        standalone = Standalone(jar_file, property_file, dij_dir, wed_dir)
        standalone.generate_plan()
        standalone.optimise_plan()
        if self._prune_spots:
            plan = Plan(plan_file)
            plan.prune_spots()
            plan.write_json()
        standalone.calculate_dose()

        # Log the run times.
        times = standalone.times
        with open(f'{workdir}/times.json', 'w') as file:
            json.dump(times, file)

        return workdir

    def run_all(self, debug: bool = False, log_dij_matrix: bool = False, log_wed: bool = False):
        '''
        Run all configurations that have been added to this runner.

        This function call is blocking.

        Returns
        -------
        Dict[str, str]
            A dictionary containing (label, result_dir) pairs.
            The final plan can be found in <result_dir>/result_plan.json, and
            the final dose distribution in <result_dir>/result_dose.dat.

        Raises
        ------
        RuntimeError
            If the standalone program finished with a non-zero exit code.
        '''
        result_dirs = {}
        for label in self._configs.keys():
            result_dir = self._run_config(label, debug=debug, log_dij_matrix=log_dij_matrix, log_wed=log_wed)
            result_dirs[label] = result_dir

        return result_dirs
