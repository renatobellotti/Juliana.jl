import logging
import struct
import numpy as np
from .grid import Grid
from .rounding import drop_precision


class CT:
    '''
    Class representing a CT image.

    Notice
    ------
    Array values are clipped to the CT range [-1024, 3071]
    '''
    allowed_orientations = ['HFS', 'HFP', 'FFS', 'FFP']

    def __init__(self,
                 data: np.ndarray,
                 spacing: np.ndarray,
                 origin: np.ndarray = np.zeros(3, dtype=np.float32),
                 orientation: str = 'HFS'):
        '''
        Parameters
        ----------
        data: np.ndarray
            The values of the CT image.
        spacing: np.ndarray
            Grid spacing.
        origin: np.ndarray
            Shape (3,). Origin of the CT. Default: (0, 0, 0)
        orientation: str
            Orientation of the patient. Can be "HFS", "HFP", "FFS",
            "FFP". Default: "HFS"

        Raises
        ------
        ValueError
            If the orientation is not in one of the following:
            "HFS", "HFP", "FFS", "FFP".
        '''
        # Clip values to the CT range.
        if (np.any(data < -1024) or np.any(data > 3071)):
            logging.warning('WARNING: Have to clip CT values to [-1024, 3071].')

        self._data = np.clip(data, -1024, 3071).astype('<i2')
        self._spacing = spacing.astype(np.float32)
        self._origin = origin.astype(np.float32)
        self._orientation = orientation

        if orientation not in self.allowed_orientations:
            raise ValueError(
                f'Orientation {orientation} is not in {self.allowed_orientations}'
            )

    def save(self, path: str):
        '''
        Saves this CT in the format readable by the standalone executable.

        Parameters
        ----------
        path: str
            Where to save the file to. A typical file extension is .dat.
        '''
        with open(path, 'wb') as file:
            file.write(bytes(self._orientation, 'ascii'))
            file.write(self._origin.tobytes())
            file.write(self._spacing.tobytes())
            file.write(np.array(self._data.shape, dtype=np.int32).tobytes())
            size = self._data.shape
            for z in range(size[2]):
                file.write(self._data[:, :, z].astype('<i2').tobytes(order='F'))

    @staticmethod
    def load(path: str):
        with open(path, 'rb') as file:
            orientation = file.read(3).decode('ascii')

            origin = np.empty(3, dtype=np.float32)
            for i in range(3):
                # little-endian 4-byte float
                origin[i] = struct.unpack('<f', file.read(4))[0]

            spacing = np.empty(3, dtype=np.float32)
            for i in range(3):
                spacing[i] = struct.unpack('<f', file.read(4))[0]

            size = np.empty(3, dtype=np.int32)
            for i in range(3):
                size[i] = struct.unpack('<i', file.read(4))[0]

            slices = []
            for z in range(size[2]):
                plane = np.frombuffer(
                    # Each pixel in the slice needs two bytes.
                    file.read(size[0]*size[1]*2),
                    dtype='<i2'
                ).reshape(size[:2], order='F')
                slices.append(plane)
            data = np.stack(slices, axis=2)

            return CT(data, spacing, origin, orientation)

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def spacing(self) -> np.ndarray:
        return self._spacing

    @property
    def origin(self) -> np.ndarray:
        return self._origin

    @property
    def orientation(self) -> str:
        return self._orientation
    
    @property
    def grid(self) -> Grid:
        return Grid(self.data.shape, self.spacing, self.origin)

    def __eq__(self, other):
        if not isinstance(other, CT):
            raise RuntimeError('Can only compare CT to ther CT instances.')

        if not np.isclose(self.spacing, other.spacing).all():
            return False
        if not np.isclose(self.origin, other.origin).all():
            return False
        if self.orientation != other.orientation:
            return False
        if not np.isclose(self.data, other.data).all():
            return False

        return True

    def _slice_to_dicom(self, slice_ind: int, decrease_precision: bool = True):
        # pydicom is an optional module: Import it only when it is needed.
        import pydicom

        # Monochrome2 is only 16bit.
        # https://pydicom.github.io/pydicom/dev/guides/encoding/rle_lossless.html
        array = (self.data[:, :, slice_ind] + 1000).astype('uint16')

        dicom_ct = pydicom.Dataset()
        dicom_ct.Modality = 'CT'
        dicom_ct.SOPClassUID = '1.2.840.10008.5.1.4.1.1.2' # CT Image Storage
        dicom_ct.ImageType = ['ORIGINAL', 'PRIMARY', 'AXIAL']
        
        dicom_ct.StudyID = '0'
        dicom_ct.SeriesNumber = 0
        dicom_ct.AcquisitionNumber = slice_ind
        dicom_ct.InstanceNumber = slice_ind

        dicom_ct.PixelData = array.tobytes(order='F')
        dicom_ct.RescaleIntercept = -1000
        dicom_ct.RescaleSlope = 1.

        # Explanation of PixelRepresentation:
        # https://dicom.innolitics.com/ciods/ct-image/image-pixel/00280103
        dicom_ct.PixelRepresentation = 0
        dicom_ct.HighBit = 15
        dicom_ct.Rows = array.shape[1]
        dicom_ct.Columns = array.shape[0]
        dicom_ct.NumberOfFrames = 1
        dicom_ct.SamplesPerPixel = 1
        dicom_ct.PixelSpacing = [
            drop_precision(self.spacing[0] * 10),
            drop_precision(self.spacing[1] * 10),
        ]

        dicom_ct.SliceThickness = self.spacing[2] * 10
        dicom_ct.SliceLocation = (self._origin[2] + self.spacing[2] * slice_ind) * 10.
        if decrease_precision:
            dicom_ct.SliceThickness = drop_precision(dicom_ct.SliceThickness)
            dicom_ct.SliceLocation = drop_precision(dicom_ct.SliceLocation)

        if decrease_precision:
            dicom_ct.ImagePositionPatient = [
                drop_precision(self._origin[0] * 10.),
                drop_precision(self._origin[1] * 10.),
                drop_precision((self._origin[2] + self.spacing[2] * slice_ind) * 10.),
            ]
        else:
            dicom_ct.ImagePositionPatient = [
                self._origin[0] * 10.,
                self._origin[1] * 10.,
                (self._origin[2] + self.spacing[2] * slice_ind) * 10.,
            ]
        # Assumption about the patient orientation!!!
        dicom_ct.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        dicom_ct.PatientPosition = self.orientation
        assert self.orientation == 'HFS' # other orientations are not tested!!!
        dicom_ct.PhotometricInterpretation = 'MONOCHROME2'
        dicom_ct.BitsAllocated = 16
        dicom_ct.BitsStored = 16
        dicom_ct.PixelRepresentation = 0 # Unsigned integer.

        dicom_ct.file_meta = pydicom.dataset.FileMetaDataset()
        dicom_ct.file_meta.magic = b'DICM'
        # Use the default "DICOM Implicit VR Little Endian Transfer Syntax"
        # https://dicom.nema.org/dicom/2013/output/chtml/part05/chapter_10.html#sect_10.1
        dicom_ct.file_meta.TransferSyntaxUID = '1.2.840.10008.1.2'
        # CT Image Storage
        dicom_ct.file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'

        return dicom_ct

    def to_dicom(
        self,
        study_instance_UID: str,
        frame_of_reference_UID: str,
        ct_series_instance_UID: str,
        patient_ID: str,
        patient_name: str,
        decrease_precision: bool = True,
    ):
        datasets = []
        for slice_ind in range(self.data.shape[2]):
            slice_ID = f'{ct_series_instance_UID}.{slice_ind}'
            ds = self._slice_to_dicom(slice_ind, decrease_precision)
            
            ds.file_meta.MediaStorageSOPInstanceUID = slice_ID
            ds.PatientName = patient_name
            ds.PatientID = patient_ID
            ds.SOPInstanceUID = slice_ID
            ds.StudyInstanceUID = study_instance_UID
            ds.SeriesInstanceUID = ct_series_instance_UID
            ds.FrameOfReferenceUID = frame_of_reference_UID

            datasets.append(ds)
        return datasets
