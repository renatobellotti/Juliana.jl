# PyFTPP

Python Wrapper for the FTPP standalone JAR.

# Reading, plotting and saving DICOM CTs and structure sets

The following code assumes that your DICOM files are stored in the directory `input_dir`.
The CT slice files follow the pattern `CT.<something>.<i>.dcm`, where `i` is the slice index.
The structure set file follows the pattern `RS.<something else>.dcm`. The code below assumes that
there is exactly one structure set file, but it can easily be adjusted to other cases as well.

```python
import os
import numpy as np
import pydicom
from pyftpp.dicom import DicomCTImage, DicomStructureSet, export_to_dicom
from pyftpp.ct import CT
from pyftpp.structure import Structure, StructureSet
from pyftpp.plotting import CTPlotter

input_dir = '/some/path'
output_dir = '/my/output/directory'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Load the CT from the DICOM files.
ct_files = [f for f in os.listdir(input_dir) if f.startswith('CT.')]
ct_files = sorted(ct_files, key=lambda f: int(f.split('.')[-2]))
ct_files = [f'{input_dir}/{f}' for f in ct_files]

img = DicomCTImage(ct_files)
ct = img.ct

# Load the structure set from the DICOM files.
ss = DicomStructureSet(structureset_file)
structure_set = ss.structure_set(ct)

# Plot (works only if executed as the only command in a Jupyter notebook cell).
%matplotlib widget
CTPlotter(ct)

# (Optional: Modify the CT and/or structures.)

# (optional) Export to DICOM.
# You could change the study instance ID to whatever you like, but make sure to be reproducible!
patient_ID = 'P12345'
study_instance_UID = pydicom.uid.generate_uid(entropy_srcs=[patient_ID])

export_to_dicom(ct, structure_set, output_dir, study_instance_UID, 'some_test')
```

# Generating UML diagrams
Run the following commands from within the repo root directory:

```bash
pyreverse -d doc/uml pyftpp/
dot -Tpdf doc/uml/classes.dot -o doc/uml/classes.pdf
```

# Running tests
You can run the unit tests by executing the following command from within the
repo top-level directory:

```bash
coverage run --omit tests/ -m unittest discover tests/
coverage html
```
This will also generate an HTML report that measures the test coverage.
Open ``htmlcov/index.html`` in a browser to view the report.

# Example code
Generate some example config files for a water phantom by running the following
commands from within ``pyftpp``:

```bash
python write_test_ct.py
python runner.py
```
This creates a directory ``example`` containing the config files for
a simple water phantom, the final optimised plan and the corresponding
dose distribution. (same format as the CT?)


# File formats
## CT file format
- Little endian; one float is 4 bytes; one int is 4 bytes; one short is 2 bytes
- First 3 bytes: ASCII characters denoting the patient orientation. See [1].
  - HFS: Head First Supine (most common [2])
  - HFP: Head First Prone
  - FFS: Feet First Supine
  - FFP: Feet first Prone

  The DICOM standard mentions more possibilities, but FTPP only implements
  these four as of 2021-07-06.
- Next 3 floats: origin
- Next 3 floats: grid spacing
- Next 3 ints: grid size = size[]
- Next size[0] * size[1] * size[2] shorts: CT values


[1] https://dicom.innolitics.com/ciods/mr-image/general-series/00185100<br />
[2] https://deepai.org/publication/a-tool-for-automatic-estimation-of-patient-position-in-spinal-ct-data


## Dose file format
- The almost the same as the CT, but:
  - There is no patient orientation.
  - The data entries are float32 values.
