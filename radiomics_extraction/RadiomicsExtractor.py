"""
RadiomicsExtractor.py

This script contains the RadiomicsExtractor class, which is used to extract 
radiomics features from medical images.

Author: Jacob Menzinga
Date: 25-01-2024
Version: 2.0
"""
# Standard Imports 
import os
from os.path import join
import time
import yaml
import logging

# External Imports
import SimpleITK as sitk
import numpy as np
from radiomics import featureextractor
import pandas as pd

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

class RadiomicsExtractor:
    """
    A class for extracting radiomics features from medical images and dose data.
    """

    def __init__(self, patient_id:str) -> None:
        """
        Initializes a RadiomicsExtractor object.
        The image and mask must have the same shape so both must be cropped.

        Parameters:
        - patient_id (str): The ID of the patient.

        Returns:
        None
        """
        self.patient_id = patient_id
        self.features_df = pd.DataFrame(index=[patient_id])
        self.extractor = featureextractor.RadiomicsFeatureExtractor()
        self.extractor.binWidth = 0.25
        
    def set_bin_width(self, bin_width):
        """
        Sets the bin width for feature extraction.

        Parameters:
        - bin_width (float): The bin width value.

        Returns:
        None
        """
        self.extractor.binWidth = bin_width
        
    def extract_features(self, masked_ct:sitk.Image, masked_dose:dict, 
                         mask:sitk.Image, return_df:bool = False) -> pd.DataFrame:
        """
        Extracts radiomics features from the masked CT image and dose data.

        Parameters:
        - masked_ct (sitk.Image): The masked CT image.
        - masked_dose (dict): The masked dose data.
        - mask (sitk.Image): The mask image.
        - return_df (bool): Whether to return the features as a DataFrame. Default is False.

        Returns:
        - pd.DataFrame: The extracted features as a DataFrame, if return_df is True.
        None, otherwise.
        """
        assert masked_ct['name'] == mask['name'] == masked_dose['name']
        
        struct_features = self.extractor.execute(masked_ct['image'], 
                                                 mask['image'])
        struct_features = pd.DataFrame.from_dict(struct_features, 
                                                 orient='index').T
        struct_features.index = [self.patient_id]
        struct_features.columns = [f'{masked_ct["name"]}_{col}'
                                      for col in struct_features.columns]
        
        dose_features = self.extract_dose_features(masked_dose)
        dose_features.index = [self.patient_id]
        dose_features.columns = [f'{masked_dose["name"]}_{col}'
                                for col in dose_features.columns]
        
        features_df = pd.concat([self.features_df, struct_features, dose_features], axis=1)   
        self.features_df = features_df
        
        if return_df:
            return features_df
    
    def extract_dose_features(self, dose:dict) -> pd.DataFrame:
        """
        Extracts dose features from the dose data.

        Parameters:
        - dose (dict): The dose data.

        Returns:
        - pd.DataFrame: The extracted dose features as a DataFrame.
        """
        array = dose['array'][~np.isnan(dose['array'])]
        dose_features = {}
        dose_features['max_dose'] = array.max()
        dose_features['mean_dose'] = array.mean()
        dose_features['std_dose'] = array.std()
        dose_features = pd.DataFrame.from_dict(dose_features, orient='index').T
        return dose_features
    
if __name__ == '__main__':
    pass
    
    
    
