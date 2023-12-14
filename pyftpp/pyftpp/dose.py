import numbers
import struct
import numpy as np
from .grid import Grid
from .interpolation import interpolate_data


class Dose:
    '''
    Class representing a dose distribution.
    '''

    def __init__(self,
                 data: np.ndarray,
                 spacing: np.ndarray,
                 origin: np.ndarray = np.zeros(3, dtype=np.float32)):
        '''
        Parameters
        ----------
        data: np.ndarray
            The values of the dose distibution.
        spacing: np.ndarray
            Grid spacing.
        origin: np.ndarray
            Shape (3,). Origin of the CT. Default: (0, 0, 0)
        '''
        # Clip values to the CT range.
        self._data = data
        self._spacing = spacing.astype(np.float32)
        self._origin = origin.astype(np.float32)

    def save(self, path: str):
        '''
        Saves this dose in the format readable by the standalone executable.

        Parameters
        ----------
        path: str
            Where to save the file to. A typical file extension is .dat.
        '''
        with open(path, 'wb') as file:
            file.write(self._origin.tobytes())
            file.write(self._spacing.tobytes())
            file.write(np.array(self._data.shape, dtype=np.int32).tobytes())
            size = self._data.shape
            for z in range(size[2]):
                file.write(self._data[:, :, z].tobytes(order='F'))

    @staticmethod
    def load(path: str):
        with open(path, 'rb') as file:
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
                    file.read(size[0]*size[1]*4),
                    dtype='float32',
                ).reshape(size[:2], order='F')
                slices.append(plane)
            data = np.stack(slices, axis=2)

            return Dose(data, spacing, origin)

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
    def grid(self) -> Grid:
        return Grid(self.data.shape, self.spacing, self.origin)

    def normalise(self, normalisation_factor: float):
        self._data /= normalisation_factor

    def __add__(self, other):
        if not isinstance(other, Dose):
            raise NotImplementedError(
                'Addition is only defined for two instances of Dose.'
            )

        same_grid = np.isclose(self.origin, other.origin).all() \
            and np.isclose(self.spacing, other.spacing).all()

        if same_grid:
            return Dose(self.data + other.data, self.spacing, self.origin)

        if self.spacing[0] <= other.spacing[0]:
            finer = self
            coarser = other
        else:
            finer = other
            coarser = self
        
        # Interpolate the values of the coarse grid on the fine grid
        # bilinearly, then add the grids normally.
        # Do it slice by slice.
        interpolated = interpolate_data(
            coarser.data,
            coarser.grid,
            finer.grid
        )
        interpolated = Dose(
            interpolated, finer.grid.spacing, finer.grid.origin
        )

        return finer + interpolated

    def __mul__(self, other):
        # Implement dose * 2.
        if not isinstance(other, numbers.Real):
            raise NotImplementedError('Scalar multiplication for Dose objects is only valid for real numbers.')

        return Dose(self.data * other, self.spacing, self.origin)

    def __rmul__(self, other):
        # Implement 2. * dose.
        return self * other

    def __truediv__(self, other):
        if not isinstance(other, numbers.Real):
            raise NotImplementedError('Scalar division for Dose objects is only valid for real numbers.')

        return Dose(self.data / other, self.spacing, self.origin)

    def to_dicom_dataset(
            self, 
            study_instance_UID: str,
            frame_of_reference_UID: str,
            series_instance_UID: str,
            patient_ID: str,
            patient_name: str,
        ):
        # pydicom is an optional module: Import it only when it is needed.
        import pydicom

        # Physical dose.
        ###################################################
        # Assumption: This dose object contains dose RBE.
        ###################################################
        data = self.data / 1.1

        # Make sure we can really fit the entire range into our ints.
        MAX_16BIT_NUMBER = 65504
        dose_scaling = (data.max() + 1) / MAX_16BIT_NUMBER
        
        # Monochrome2 is only 16bit.
        # https://pydicom.github.io/pydicom/dev/guides/encoding/rle_lossless.html
        array = (data / dose_scaling).astype('uint16')

        dicom_dose = pydicom.Dataset()
        dicom_dose.Modality = 'RTDOSE'
        dicom_dose.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.2' # RT Dose Storage
        dicom_dose.PixelData = array.tobytes(order='F')
        dicom_dose.DoseGridScaling = dose_scaling
        # Explanation of PixelRepresentation:
        # https://dicom.innolitics.com/ciods/ct-image/image-pixel/00280103
        dicom_dose.PixelRepresentation = 0
        dicom_dose.HighBit = 15
        dicom_dose.Rows = data.shape[0]
        dicom_dose.Columns = data.shape[1]
        dicom_dose.NumberOfFrames = self.data.shape[2]
        dicom_dose.PixelSpacing = [self.spacing[0] * 10, self.spacing[1] * 10]
        dicom_dose.SliceThickness = self.spacing[2] * 10
        dicom_dose.GridFrameOffsetVector = [i*self.spacing[2]*10 for i in range(self.data.shape[2])]

        dicom_dose.ImagePositionPatient = (self.origin * 10.).tolist()
        # Assumption about the patient position!!!
        dicom_dose.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        dicom_dose.PhotometricInterpretation = 'MONOCHROME2'
        dicom_dose.DoseUnits = 'GY'
        dicom_dose.DoseType = 'PHYSICAL'
        dicom_dose.DoseSummationType = 'PLAN'
        dicom_dose.BitsAllocated = 16
        dicom_dose.BitsStored = 16
        dicom_dose.SamplesPerPixel = 1
        dicom_dose.file_meta = pydicom.dataset.FileMetaDataset()
        dicom_dose.file_meta.magic = b'DICM'
        # Use the default "DICOM Implicit VR Little Endian Transfer Syntax"
        # https://dicom.nema.org/dicom/2013/output/chtml/part05/chapter_10.html#sect_10.1
        dicom_dose.file_meta.TransferSyntaxUID = '1.2.840.10008.1.2'
        dicom_dose.file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.2' # RT Dose Storage
        dicom_dose.file_meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid()

        dicom_dose.PatientName = patient_name
        dicom_dose.PatientID = patient_ID
        dicom_dose.StudyInstanceUID = study_instance_UID
        dicom_dose.SeriesInstanceUID = series_instance_UID
        dicom_dose.SOPInstanceUID = series_instance_UID
        dicom_dose.FrameOfReferenceUID = frame_of_reference_UID

        return dicom_dose


def normalise_dose(dose: Dose, target_mask: np.ndarray, target_dose: float) -> Dose:
    '''
    Normalise the dose distribution s. t. the mean dose to the given binary
    mask is 100% of the given target dose.
    '''
    D_mean = (dose.data * target_mask).sum() / np.sum(target_mask)
    return dose * target_dose / D_mean
