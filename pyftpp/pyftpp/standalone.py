import logging
from pathlib import Path
import subprocess
import time


def get_path_to_bin():
    return str(Path(__file__).parent / '..' / 'bin')


def get_path_to_jar():
    return f'{get_path_to_bin()}/ch.psi.ftpp.standalone.planner-1.0.6.jar'


class Standalone:
    '''
    Class representing the standalone executable.

    This class does not handle the generation of the configuration files,
    but assumes that they have already been written.

    The execution time per step is measured.
    '''

    def __init__(self, jar_file: str, properties_file: str, dij_dir: str = None, wed_dir: str = None):
        '''
        Parameters
        ----------
        jar_file: str
            Path to the JAR file of the standalone executable.
        properties_file: str
            Path to the .properties file containing all the information.
        dij_dir: str or None
            Path to the directory where the Dij matrix should be written.
            No Dij information is written if the value is None.
        wed_dir: str or None
            Path to the directory where the WED vectors for each field should
            be written. No WED information is written if the value is None.
        '''
        self._jar_file = str(Path(jar_file).resolve())
        self._properties_file = str(Path(properties_file))
        self._log_dij_matrix_argument = f'-dijLogDir {dij_dir}' if dij_dir else ''
        self._log_wed_matrix_argument = f'-wedLogDir {wed_dir}' if wed_dir else ''


        self._time_run = None
        self._time_generation = None
        self._time_optimisation = None
        self._time_dose_calculation = None

    def run(self):
        start = time.time()
        cmd = f'java -jar {self._jar_file} {self._log_dij_matrix_argument} {self._log_wed_matrix_argument} {self._properties_file}'
        logging.info(f'PyFTPP Run: {cmd}')
        subprocess.run(cmd.split(), check=True, capture_output=True)
        end = time.time()
        self._time_run = end - start

    def generate_plan(self):
        start = time.time()
        cmd = f'java -jar {self._jar_file} -c GENERATE_PLAN {self._properties_file}'
        logging.info(f'PyFTPP Generate plan: {cmd}')
        subprocess.run(cmd.split(), check=True, capture_output=True)
        end = time.time()
        self._time_generation = end - start

    def optimise_plan(self):
        start = time.time()
        cmd = f'java -jar {self._jar_file} -c OPTIMIZE {self._log_dij_matrix_argument} {self._log_wed_matrix_argument} {self._properties_file}'
        logging.info(f'PyFTPP Optimise plan: {cmd}')
        subprocess.run(cmd.split(), check=True, capture_output=True)
        end = time.time()
        self._time_optimisation = end - start

    def calculate_dose(self):
        start = time.time()
        cmd = f'java -jar {self._jar_file} -c CALCULATE_DOSE {self._properties_file}'
        logging.info(f'PyFTPP calculate dose: {cmd}')
        subprocess.run(cmd.split(), check=True, capture_output=True)
        end = time.time()
        self._time_dose_calculation = end - start

    @property
    def times(self):
        return {
            'time_run': self._time_run,
            'time_generation': self._time_generation,
            'time_optimisation': self._time_optimisation,
            'time_dose_calculation': self._time_dose_calculation,
        }
