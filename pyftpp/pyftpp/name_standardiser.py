import numpy as np
import pandas as pd


class NameStandardiser:

    def __init__(self, structure_name_file: str, patient_ID: str):
        '''
        Parameters
        ----------
        structure_name_file: str
            Path to the CSV file containing the structure name lookup table.
            The CSV file must contain at least the following columns:
            "patient_ID", "original_name", "standardised_name"
        '''
        lookup = pd.read_csv(structure_name_file).loc[:, ['patient_ID', 'original_name', 'standardised_name']]
        lookup.query('`patient_ID` == @patient_ID', inplace=True)
        lookup.query('`standardised_name` != ["unknown", "irrelevant"]', inplace=True)
        lookup.drop_duplicates(inplace=True)
        lookup = [tuple(r) for r in lookup.to_numpy()]
        self._lookup = lookup

    def __call__(self, original_name: str) -> str:
        '''
        Returns
        -------
        str
            The standardised structure name.
        '''
        for patient_ID, unstandardised_name, standardised_name in self._lookup:
            if original_name == unstandardised_name:
                return standardised_name
        #raise RuntimeError(f'No standardised name found for structure {original_name}')
        return np.nan
