import os
import numpy as np


# Drop everything beyond 6 digits of precision by default.
MEDICAL_FLOAT_DECIMALS = os.environ.get('MEDICAL_FLOAT_PRECISION', 6)


def drop_precision(a, decimals=MEDICAL_FLOAT_DECIMALS):
    return np.around(a, decimals=decimals)
