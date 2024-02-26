"""
A module that contains classes for loading in CT and PET scans, and RTDOSE and 
RTSTRUCT files and preparing them for further processing.
Parts of this code are a refectorisation of code written originally by 
Henriëtte Kuipers

By: Jacob Menzinga
Version: 1.6
"""
# Standard imports
import glob
import os
import logging
from datetime import datetime

# External imports
import numpy as np
import pydicom as pdcm
import SimpleITK as sitk
from scipy.ndimage.interpolation import shift

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

class RtdoseCtArrayNotEqual(Exception):
    'raised when the rtdose and ct arrays are not equal in size'

class CTLoader:
    """
    A class for loading CT-scans and preparing them for further processing

    Attributes:
    - ct_path (str): The path to the CT scan.
    - slices (list): A list of DICOM slices.
    - array (numpy.ndarray): A NumPy array representation of the CT scan.
    - metadata (dict): A dictionary containing the direction, origin, and 
        spacing of the DICOM series.
    - image (SimpleITK.Image): The converted image object.

    Args:
    - ct_path (str): The path to the CT scan.

    Methods:
    - make_array(): Reads DICOM files from the specified directory and returns a
        tuple containing the DICOM slices and a NumPy array representation of 
        the CT scan.
    - collect_metadata(): Collects metadata for the loaded DICOM series.
    - convert_to_hu(): Convert the pixel intensity of the CT scan to Hounsfield
        Units (HU) using the RescaleIntercept and RescaleSlope values of the 
        first slice in the scan. The HU values are returned as a numpy array of
        dtype int16.
    - array_to_image(): Converts a numpy array to a SimpleITK image object.
    """

    def __init__(self, ct_path) -> None:
        """
        Initializes the CTLoader class with the path to the CT scan.
        """
        self.ct_path = ct_path
        self.slices, self.array = self.make_array()
        self.metadata = self.collect_metadata()
        self.array = self.convert_to_hu()
        self.image = self.array_to_image()

    def make_array(self):
        """
        Reads DICOM files from the specified directory and returns a tuple
        containing the DICOM slices and a NumPy array representation of the CT 
        scan.

        Returns:
        - slices: A list of DICOM slices.
        - ct_array: A NumPy array representation of the CT scan.
        """

        slices = [os.path.join(self.ct_path, dcm)
                  for dcm in os.listdir(self.ct_path)]
        slices = [pdcm.read_file(dcm) for dcm in slices]
        slices.sort(key=lambda x: float(x.ImagePositionPatient[2]),
                    reverse=True
                    )

        ct_array = np.stack([s.pixel_array for s in slices], axis=0)
        ct_array = np.flip(ct_array, axis=0)

        return slices, ct_array

    def collect_metadata(self):
        """
        Collects metadata for the loaded DICOM series.

        Returns:
        - metadata (dict): A dictionary containing the direction, origin, and
            spacing of the DICOM series.
        """
        metadata = dict()

        # Calculated metadata, used for converting pixel intensity to HU

        spacing = [float(i) for i in self.slices[0].PixelSpacing] \
                   + [float(self.slices[0].SliceThickness)]
        
        z_indices = [round(s.ImagePositionPatient[2], 1) for s in self.slices]

        direction = list(self.slices[0].ImageOrientationPatient) + [0, 0, 1]

        metadata['direction'] = direction
        metadata['origin'] = self.slices[0].ImagePositionPatient
        metadata['spacing'] = spacing
        metadata['patient_id'] = self.slices[0].PatientID
        metadata['z_indices'] = z_indices
        metadata['UID'] = self.slices[0].FrameOfReferenceUID

        # TODO: build in some checks here, look in A1_RTdose for inspiration

        return metadata

    def convert_to_hu(self):
        """
        Convert the pixel intensity of the CT scan to Hounsfield Units (HU) 
        using the RescaleIntercept and RescaleSlope values of the first slice in
        the scan. The HU values are returned as a numpy array of dtype int16.

        Returns:
        - hu_ct_array (numpy.ndarray): 
            A 3D numpy array of Hounsfield Units (HU) repdresenting the CT scan.
        """
        intercept = self.slices[0].RescaleIntercept
        slope = self.slices[0].RescaleSlope
        hu_ct_array = (self.array * slope + intercept).astype(np.int16)

            
        return hu_ct_array

    def array_to_image(self):
        """
        Converts a numpy array to a SimpleITK image object.

        Returns:
        - ct_image (SimpleITK.Image): The converted image object.
        """
        ct_image = sitk.GetImageFromArray(self.array)
        ct_image.SetSpacing(self.metadata['spacing'])
        ct_image.SetOrigin(self.metadata['origin'])
        ct_image.SetDirection(self.metadata['direction'])

        return ct_image


class RTDOSELoader:

    """
    A class for loading and processing CT and RTDOSE images.

    Args:
        - rtdose_path (str): The file path to the RTDOSE image.
        - ct_img (SimpleITK.Image): The SimpleITK image object for the CT image.
        - ct_arr (np.ndarray): The numpy array for the CT image.

    Attributes:
        - rtdose_path (str): The file path to the RTDOSE image.
        - ct_img (SimpleITK.Image): The SimpleITK image object for the CT image.
        - ct_arr (np.ndarray): The numpy array for the CT image.
        - image (SimpleITK.Image): The SimpleITK image object for the RTDOSE 
            image.
        - metadata (dict): A dictionary containing metadata for the RTDOSE 
            image.
        - resampled_rtdose_img (SimpleITK.Image): The resampled RTDOSE image.
        - resampled_rtdose_arr (np.ndarray): The numpy array for the 
            resampled RTDOSE image.
    """

    def __init__(self, rtdose_path: str, ct:CTLoader) -> None:
        self.rtdose_path = rtdose_path
        self.ct = ct
        self.dcm = self.load_rtdose()
        self.metadata = self.collect_metadata()
        self.array, self.img = self.create_array_and_image()
        self.resampled_array, self.resampled_img = self.resample_rtdose_to_ct()

    def load_rtdose(self):
        """
        Loads the RTDOSE file specified by the rtdose_path attribute.

        If the file has a .DCM extension, it is loaded using the SimpleITK 
        ImageFileReader. Otherwise, the file is assumed to be a DICOM series
        and is loaded using the SimpleITK ImageSeriesReader. The last slice
        of the loaded image is returned.

        Returns:
            - image (SimpleITK.Image): The loaded RTDOSE image.
        """

        if self.rtdose_path.upper().endswith('.DCM'):
            dose_dcm = pdcm.read_file(self.rtdose_path)
            
            return dose_dcm
        
        else:
            dose_path = os.path.join(self.rtdose_path, 
                                    os.listdir(self.rtdose_path)[0])
            dose_dcm = pdcm.read_file(dose_path)
            
            return dose_dcm

    def collect_metadata(self):
        """
        Collects metadata from the DICOM file at the specified RTDOSE path.

        Returns:
            - metadata (dict): A dictionary containing the following 
                metadata:
                - PatientID (str): The ID of the patient associated with 
                    the DICOM file.
                - DoseGridScaling (float): The scaling factor for the dose
                    grid.
                - DoseUnits (str): The units of the dose values.
        """
        metadata = dict()
        metadata['PatientID'] = self.dcm.PatientID
        metadata['DoseGridScaling'] = float(self.dcm.DoseGridScaling)
        metadata['DoseUnits'] = self.dcm.DoseUnits
        metadata['ImagePositionPatient'] = [float(pos) for pos in 
                                            self.dcm.ImagePositionPatient]
        metadata['Spacing']  = [f for f in self.dcm.PixelSpacing] +\
                               [self.dcm.SliceThickness]

        return metadata
    
    def create_array_and_image(self):
        dose_array = self.dcm.pixel_array * self.metadata['DoseGridScaling']
        dose_img = sitk.GetImageFromArray(dose_array)
        dose_img.SetSpacing(self.metadata['Spacing'])

        return dose_array, dose_img

    def resample_rtdose_to_ct(self):
        """
        Resamples the RTDOSE image to the same grid as the CT image and converts
        the pixel values to dose values.

        Returns:
        - resampled_rtdose_img (SimpleITK.Image): The resampled RTDOSE image.
        - resampled_rtdose_arr (numpy.ndarray): The resampled RTDOSE array with 
            pixel values converted to dose values.
        """
        resampler = sitk.ResampleImageFilter()
        resampler.SetOutputSpacing(self.ct.metadata['spacing'])
        resampler.SetSize(self.ct.array.shape[::-1])
        resampled_img = resampler.Execute(self.img)
        resampled_array = sitk.GetArrayFromImage(resampled_img)

        if self.ct.array.shape != resampled_array.shape:
            raise RtdoseCtArrayNotEqual(
                'The CT and RTDOSE arrays are not equal in size')
            
        # Shifting the rtdose array to match the ct array
        ct_imgpos_patMax = [float(each) for each in self.ct.slices[-1].ImagePositionPatient]
        dy, dx, dz = (np.array(self.metadata['ImagePositionPatient']) - np.array(ct_imgpos_patMax)) / np.array(self.ct.metadata['spacing'])
        
        resampled_array = shift(resampled_array, (dz, dx, dy))
        resampled_img = sitk.GetImageFromArray(resampled_array)
        
        return resampled_array, resampled_img
        
class RTSTRUCTLoader:
    """
    A class to load and extract information from RTSTRUCT files.
    """

    def __init__(self, rtdose_path: str, ct:CTLoader) -> None:
        """
        Initializes an RTSTRUCTLoader object.

        Args:
            - rtdose_path (str): The path to the RTSTRUCT file.
        """
        self.rtdose_path = rtdose_path
        self.dcm = self.load_rtstruct()
        self.ct = ct        
        self.contours = self.get_contours()
        
    def load_rtstruct(self):
        """
        Loads an RTSTRUCT file and returns a pydicom object.

        If the file path ends with '.DCM', it is assumed to be a single file 
        and is loaded directly. Otherwise, the function searches for all '.dcm' 
        files in the directory and loads the first one found.

        Returns:
            - A pydicom object representing the loaded RTSTRUCT file.
        """
        if self.rtdose_path.upper().endswith('.DCM'):
            dcm = pdcm.read_file(self.rtdose_path)
            return dcm

        else:
            dcm = [os.path.join(self.rtdose_path, dcm)
                   for dcm in os.listdir(self.rtdose_path)]
            dcm = pdcm.read_file(dcm[0])
            return dcm

    def get_contours(self):
    
        """
        Extracts contour information from the loaded RTSTRUCT file.

        Returns:
            - A nested dictionary containing the following information for each
                Contour:
                - name (str): The name of the contour.
                - number (int): The number of the contour.
                - contour_sequence (numpy.ndarray): A numpy array containing the
                    contour coordinates for each slice.
        """
        contours =  dict()
        
        # Creating dicts for referencing ROI's and ContourSequences
        rois = dict()
        for ss in self.dcm.StructureSetROISequence:
            try:
                rois[ss.ROINumber] = ss.ROIName

            except AttributeError:
                continue
            
        contour_sequences = dict()
        
        for rcs in self.dcm.ROIContourSequence:
            try:
                contour_sequences[rcs.ReferencedROINumber] = rcs.ContourSequence
                
            except Exception:
                pass

        target = dict([
            ('lung_left', ['Lung_L', 'longli', 'long_li', 'Lung_Li', 'long-li',
                            'long li', 'Lung (Left)']),
            ('lung_right', ['Lung_R', 'longre', 'long_re', 'Lung_Re',
                            'long-re', 'long re','Lung (Right)']),
            ('lungs', ['Lungs', 'longen', 'Lung']),
            ('gtv_lung', ['Longen-GTV', 'Lung minus GTV', 'longen-gtv',
                            'Lung-GTV', 'long-gtv', 'longen_gtv', 
                            'Lung minus GTV']),
            ('heart', ['Heart', 'hart']),
            ('itv', ['ITV', '2ITV-tumor', 'ITV-tumor', 'ITV_klieren',
                        'ITVtumor', 'ITVklieren', 'ITVtu', 'IGTV']),
            ('gtv', ['GTV', 'GTVp','2GTV-tumor', 'GTV-tumor', 'GTV_klieren',
                        'GTVtumor', 'GTVklieren', 'GTVpostop'])
            ])

        rois_of_interest = {'lung_left': None, 'lung_right': None,
                    'lungs': None, 'gtv_lung': None, 'heart': None,
                    'itv': None, 'gtv': None}

        # Collecting the contour sequences of the ROIs of interest
        for num, roi in rois.items():
            for target_name, variations in target.items():
                
                if roi in variations:
                    if rois_of_interest[target_name] is not None:
                        pass
                    else:
                        rois_of_interest[target_name] = num
                        break
        
        for name, num in rois_of_interest.items():
            if num is not None:
                contour = dict()
                contour['name'] = name
                contour['number'] = num

                con_sequence = tuple()
                try:
                    for con_per_slice in contour_sequences[num]:           
                    
                        condata_slice = np.asarray(con_per_slice.ContourData)
                        x = condata_slice[range(0, len(condata_slice), 3)]
                        y = condata_slice[range(1, len(condata_slice), 3)]
                        z = condata_slice[range(2, len(condata_slice), 3)]

                        condata_pts_slice = np.array([x, y, z]).T
                        
                        cont_coord_perslice = (
                            (condata_pts_slice - self.ct.metadata['origin'])
                            / self.ct.metadata['spacing'])
                        cont_coord_perslice = np.round(cont_coord_perslice).astype(int)
                        
                        
                        con_sequence += (cont_coord_perslice,)
                        
                    contour['contour_sequence'] = np.array(con_sequence, 
                                                           dtype='object')
                    contours[name] = contour
                except Exception as e:
                    logger.error(e)
                    logger.warning(f"{self.ct.metadata['patient_id']} - No contour found for {name}")
                    continue
                
            if num is None:
                logger.warning(f"{self.ct.metadata['patient_id']} - No contour found for {name}")
                continue    
            
        return contours
    
class PETLoader:
    def __init__(self, pet_path) -> None:
        self.pet_path = pet_path
        self.slices, self.array = self.make_array()
        self.metadata = self.collect_metadata()

    def make_array(self):
        slices = glob.glob(self.pet_path + '\*.dcm')
        slices = [pdcm.read_file(dcm) for dcm in slices]
        slices.sort(key=lambda x: float(x.ImagePositionPatient[2]))
        pet_array = np.stack([s.pixel_array for s in slices], axis=0)
        return slices, pet_array

    def collect_metadata(self):
        metadata = dict()
        # Patient information
        metadata['patient_id'] = self.slices[0].PatientID
        metadata['patient_weight'] = self.slices[0].PatientWeight*1000
        if (metadata['patient_weight'] < 20000 or
                metadata['patient_weight'] == None):
            logger.warning('%s weight is anomalous', metadata['patient_id'])

        # Radiopharmacutical information
        rad_pharm_seq = self.slices[0].RadiopharmaceuticalInformationSequence[0]

        metadata['radiopharmaceutical'] = rad_pharm_seq.Radiopharmaceutical
        metadata['injected_dose'] = float(rad_pharm_seq.RadionuclideTotalDose)
        metadata['units'] = self.slices[0].Units
        metadata['half_life'] = float(rad_pharm_seq.RadionuclideHalfLife)
        metadata['decayconstant'] = np.log(2)/metadata['half_life']

        rad_pharm_start = (self.slices[0].SeriesDate +
                           rad_pharm_seq.RadiopharmaceuticalStartTime)
        rad_pharm_start = datetime.strptime(rad_pharm_start, '%Y%m%d%H%M%S.%f')
        metadata['radiopharmaceutical_start'] = rad_pharm_start
        # Here also a check could be made for radpharmstart < seriesstart

        metadata['rescale_slope'] = float(self.slices[0].RescaleSlope)
        metadata['rescale_intercept'] = float(self.slices[0].RescaleIntercept)

        # Acquisition information
        metadata['acquisition_datetime'] = datetime.strptime(
            self.slices[0].AcquisitionDate + self.slices[0].AcquisitionTime,
            '%Y%m%d%H%M%S.%f')
        metadata['series_datetime'] = datetime.strptime(
            self.slices[0].SeriesDate + self.slices[0].SeriesTime,
            '%Y%m%d%H%M%S')

        return metadata

    def calculate_suv(self):
        suv_array = (self.array * self.metadata['rescale_slope']
                     + self.metadata['rescale_intercept'])

        decay_time = (self.metadata['radiopharmaceutical_start']
                      - self.metadata['acquisition_datetime']).total_seconds()

        decayed_dose = self.metadata['injected_dose']*(
            2**(-decay_time/self.metadata['half_life']))

        suv_bwscalefactor = (self.metadata['patient_weight']/decayed_dose)

        suv_bw_array = suv_array*suv_bwscalefactor

        return suv_bw_array

    # the pet loader did not get fully finished because it became out of scope 
    # for the project. I'm leving it in here so someone might finish it in the
    # future or at the very least, take inspiration from it.

if __name__ == '__main__':
    pass
