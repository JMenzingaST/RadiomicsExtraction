"""
This script contains the cropper class, used for cropping 3D arrays
"""

# External Import
import numpy as np

class Cropper:
    @classmethod
    def get_cropping_bounds(self, array) -> tuple:
        """
        Get the cropping bounds for a given mask.

        Parameters:
        mask (ndarray): The mask array.

        Returns:
        tuple: A tuple containing the minimum and maximum values for each axis 
        (x_min, x_max, y_min, y_max, z_min, z_max).
        """
        array = np.nan_to_num(array)
        # Get the non-zero slices along each axis to define the cropping region
        x_non_zero = np.where(~np.all(np.all(array == 0, axis=2), axis=1))[0]
        y_non_zero = np.where(~np.all(np.all(array == 0, axis=2), axis=0))[0]
        z_non_zero = np.where(~np.all(np.all(array == 0, axis=1), axis=0))[0]

        # Expand the cropping region by 1 pixel in each direction
        x_min, x_max = x_non_zero.min() - 1, x_non_zero.max() + 2
        y_min, y_max = y_non_zero.min() - 1, y_non_zero.max() + 2
        z_min, z_max = z_non_zero.min() - 1, z_non_zero.max() + 2

        return x_min, x_max, y_min, y_max, z_min, z_max
    
    @classmethod
    def crop_3d_array(self, mask, masked_array=None):
        """
        Crop a 3D array based on the provided mask and return the cropped array.

        Args:
            mask (dict): The mask dictionary containing the 'name' and 'array' keys.
            masked_array (dict, optional): The masked array dictionary 
            containing the 'name' and 'array' keys. Defaults to None.

        Returns:
            dict: The cropped array dictionary containing the 'name' and 'array' keys.
        """
        if masked_array:
            cropping_bounds = self.get_cropping_bounds(mask['array'])
            array = masked_array['array'][cropping_bounds[0]:cropping_bounds[1],
                                          cropping_bounds[2]:cropping_bounds[3],
                                          cropping_bounds[4]:cropping_bounds[5]]
            
            return {'name': masked_array['name'], 'array': array}
        
        if not masked_array:    
            cropping_bounds = self.get_cropping_bounds(mask['array'])
            array = mask['array'][cropping_bounds[0]:cropping_bounds[1],
                                  cropping_bounds[2]:cropping_bounds[3],
                                  cropping_bounds[4]:cropping_bounds[5]]
        
            return {'name': mask['name'], 'array': array}
        
        
    
