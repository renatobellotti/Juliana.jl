from abc import ABC, abstractmethod
import os
from typing import Iterator, List, Tuple, Union
from .ct import CT
from .grid import Grid
from .structure import Structure, StructureSet, FionaStructureSet
from .prescriptions import Prescriptions, load_fiona_prescriptions


class DataSet(ABC):

    @property
    @abstractmethod
    def structures(self) -> StructureSet:
        pass

    @property
    @abstractmethod
    def prescriptions(self) -> Prescriptions:
        pass

    @property
    @abstractmethod
    def ct(self) -> CT:
        pass


# list of (patient_ID, original_name, standardised_name)
LookupType = List[Tuple[str, str, str]]


class StandardisedNameStructureSet(StructureSet):
    '''
    Wrapper for another dataset that allows to access structures by another
    name, e. g. a standardised name.
    There will be an exception if the given name does not match a unique
    name in the underlying StructureSet.
    '''

    def __init__(self, structure_set: StructureSet, lookup: LookupType):
        self._other = structure_set
        self._lookup = lookup

        self._name_cache = {}
    
    def _lookup_name(self, standardised_name: str) -> str:
        if standardised_name in self._name_cache:
            return self._name_cache[standardised_name]

        matches = [entry[1] for entry in self._lookup if entry[2] == standardised_name]

        if len(matches) == 0:
            raise RuntimeError(f'No standardised name found for {standardised_name}')
        elif len(matches) > 1:
            raise RuntimeError(f'Multiple standardised names found for {standardised_name}!')
        else:
            self._name_cache[standardised_name] = matches[0]
            return matches[0]
    
    def __getitem__(self, key: str) -> Structure:
        original_name = self._lookup_name(key)
        return self._other[original_name]
    
    @property
    def structure_names(self) -> Iterator[str]:
        standard_names = set(standard_name for patient_ID, original, standard_name in self._lookup)
        for name in sorted(list(standard_names)): 
            yield name

    @property
    def ct(self) -> Union[None, CT]:
        return self._other.ct


class FionaPatient(DataSet):
    '''
    Class for loading data for the Fiona standalone executable.

    The data directory contains all the patient data. It is supposed to have
    the following structure:

    data_dir/
        dose_per_target.csv
        CTs/
            P000000_0.dat
            P000000_1.dat
            ...
            P000001_0.dat
            ...
        prescriptions/
            P000000_constraints.csv
            ...
        structures/
            P000000/
                0/  # Series 0
                    structure1.npy
                    structure2.npy
                    ...
                1/  # Series 1
                    ...
            ...

    The CT file name follows the scheme "<patient_ID>_<series_ID>.dat".

    A constraints CSV file must have the following columns:
        structure_name
        soft_hard            (can take values "soft" or "hard")
        constraint_quantity  ("D_max", "D_mean" or "Dn", where n is an integer)
        constraint_threshold (a dose threshold not to be exceeded; can be given in units of Gy or % of the target dose)

    The file dose_per_target.csv must have the following columns:
        patient_ID
        target_name
        total_dose_gy
        dose_gy_per_frac
    '''

    def __init__(self,
                 data_dir: str,
                 patient_ID: str,
                 series_ID: int):
        self._data_dir = data_dir
        self._patient_ID = patient_ID
        self._series_ID = series_ID

        self._ct_path = f'{data_dir}/CTs/{patient_ID}_{series_ID}.dat'
        self._structure_dir = f'{data_dir}/structures/{patient_ID}/{series_ID}'

        self._ct = None
        self._structureset = None
        self._prescriptions = None

    @property
    def patient_ID(self) -> str:
        return self._patient_ID
    
    @property
    def series_ID(self) -> int:
        return self._series_ID

    @property
    def prescriptions(self) -> Prescriptions:
        if self._prescriptions is None:
            self._prescriptions = self._load_prescriptions()
        return self._prescriptions

    @property
    def structures(self) -> StructureSet:
        if self._structureset is None:
            self._structureset = FionaStructureSet(
                self._structure_dir, self.ct
            )
        return self._structureset
    
    @property
    def ct(self) -> CT:
        if self._ct is None:
            self._ct = CT.load(self._ct_path)
        return self._ct
    
    @property
    def ct_path(self) -> str:
        return self._ct_path

    def _load_prescriptions(self) -> Prescriptions:
        masks = {}

        return load_fiona_prescriptions(
            self.patient_ID,
            f'{self._data_dir}/dose_per_target.csv',
            f'{self._data_dir}/prescriptions/{self._patient_ID}_constraints.csv',
            self.ct.spacing,
            self.structures,
        )
