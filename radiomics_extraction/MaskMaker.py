"""
This module provides a class for generating lung masks, creating masks from 
contours, manipulating masks, and save them as SITK images.

Author: Jacob Menzinga
Version: 2.0
"""
# Standard imports
import os
from os.path import join

# External imports
import numpy as np
import SimpleITK as sitk
from skimage.draw import polygon
from lungmask import LMInferer

class MaskMaker:
    @classmethod
    def generate_lungmask(self, ct):
        """
        Generates lung masks based on the input CT image.

        Args:
            ct (CTImage): The input CT image.

        Returns:
            tuple: A tuple containing two dictionaries representing the lung masks.
                The first dictionary represents the right lung mask, 
                the second dictionary represents the left lung mask
                Each dictionary has two keys: 'name' (the name of the mask) and
                'array' (the mask array).
        """
        lm_inferer = LMInferer()
        lung_mask = lm_inferer.apply(ct.image)

        l_re = np.where(lung_mask == 1, 1, np.nan)
        l_li = np.where(lung_mask == 2, 1, np.nan)

        l_re[l_re == 0] = np.nan
        l_li[l_li == 0] = np.nan

        l_re = {'name': 'lung_right', 'array': l_re}
        l_li = {'name': 'lung_left', 'array': l_li}

        return l_re, l_li       
    
    @classmethod  
    def make_mask(self, contour, ct, fill=True):
        """
        Create a mask from a contour.

        Args:
            contour (dict): The contour information.
            ct (CTLoader): The ct object.
            fill (bool, optional): Whether to fill the mask. Defaults to True.

        Returns:
            dict: The mask information.

        """
        if fill:
            mask = np.zeros(ct.array.shape)
            for i in range(len(contour['contour_sequence'])):
                x, y, z = zip(*contour['contour_sequence'][i][:, :])
                # Spacing, maybe.
                rr, cc = polygon(y, x)
                mask[z[0], rr, cc] = 1
                
        else:
            mask = np.zeros(ct.array.shape)
            for i in range(len(contour['contour_sequence'])):
                x, y, z = zip(*contour['contour_sequence'][i][:, :])

                mask[z[0], y, x] = 1
                
        mask = np.roll(mask, -1, axis=0)
        mask[mask == 0] = np.nan
        
        return {'name': contour['name'], 'array': mask}
        
    @classmethod    
    def manipulate_mask(self, mask1, mask2, operation='add'):
        """ 
        Adds two masks together.

        Args:
            - mask1 (numpy array): The first mask.
            - mask2 (numpy array): The second mask.
            - operation (str): The operation to be performed on the masks. Must be
                either 'add' or 'subtract'.
                
        Returns:
            - The sum of the two masks.
        """
        # Cant add or subtract nan's, so convert them to 0.
        m1 = np.nan_to_num(mask1['array'])
        m2 = np.nan_to_num(mask2['array'])
        
        if operation == 'add':
            combined_mask = np.clip(m1 + m2, 0, 1)
            combined_mask[combined_mask == 0] = np.nan
            
            return {'name': f'{mask1["name"]}+{mask2["name"]}',
                    'array': combined_mask}
            
        elif operation == 'subtract':
            combined_mask = np.clip(m1 - m2*2, 0, 1)
            # pretty sure the *2 is unnecessary.
            combined_mask[combined_mask == 0] = np.nan
            
            return {'name': f'{mask1["name"]}-{mask2["name"]}',
                    'array': combined_mask}
            
        else:
            raise ValueError('Invalid operation, must be "add" or "subtract"')
    @classmethod    
    def create_masked_array(self, array, mask):
        """
        Creates a masked ct array from a ct and a mask.

        Args:
            - array (ndarray): The ct / dose array.
            - mask (dict): The mask to be applied to the ct.

        Returns:
            - masked_array (dict): A dict containing the structure name and the
                                   corresponding masked array.
        """
        try:
            masked_array = array * mask['array']
            masked_array = {'name': mask['name'], 'array': masked_array}
            return masked_array  
    
        except Exception as e:
            print(f'Could not create masked ct array for {mask["name"]} due to {e}')
        
    @classmethod
    def create_image(self, array_dict):
        """
        Creates images from the given array
        In the pipeline, this is used to create cropped images.
        
        Args:
            - array_dict (dict): A dictionary containing the structure name and 
                the corresponding mask array.
            
        Returns:
            - image (dict): A dict containing the structure name and the 
                corresponding sitk image.
        """
        
        name = array_dict['name']
        image = sitk.GetImageFromArray(array_dict['array'])
        
        return {'name': name, 'image': image}
        
    @classmethod
    def save_image(self, image, type, patient_id, folder = 'cropped_structures'):
        """
        Saves the given image to the specified folder.
        
        Args:
            - image (dict): A dict containing the structure name and the 
                corresponding sitk image.
            - type (str): The type of image, either 'mask' or 'ct'.
            - folder (str): The path to the folder where the images will be 
                saved.
            
        """                
        try:
            if not os.path.exists(folder):
                os.makedirs(folder)
                
            sitk.WriteImage(
                image['image'], 
                join(folder, 
                     f"{str(patient_id)}_{type}_{image['name']}.nii.gz"))
        
        except Exception as e:
            print(f'Could not save {type} {image["name"]} for patient {patient_id} due to {e}')
            
if __name__ == "__main__":
    pass