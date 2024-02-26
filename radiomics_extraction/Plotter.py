"""
This script contains functions for plotting CT, dose, and mask images for a given slice number and patient ID.

Author: Jacob Menzinga
Version: 1.1
"""

# Standard Imports
import os
from os.path import join

# External Imports
from Loaders import CTLoader, RTSTRUCTLoader, RTDOSELoader
from MaskMaker import MaskMaker
import matplotlib.pyplot as plt
import numpy as np

def plot(slice_nr:int, patient_id:int, struct_name:str, ct:CTLoader=False, 
         dose:RTDOSELoader=False, mask=False, masked_array=False,
         mask_fill:bool=True, save:bool=False, save_path=None):
    """
    Plot the CT, dose, and mask images for a given slice number and patient ID.
    
    Parameters:d
        slice_nr (int): The slice number to plot.
        patient_id (int): The ID of the patient.
        ct (CTLoader, optional): The CTLoader object containing the CT image. Defaults to False.
        dose (RTDOSELoader, optional): The RTDOSELoader object containing the dose image. Defaults to False.
        mask (mask dict from MM, optional): Whether to plot the mask image. Defaults to False.
        mask_fill (bool, optional): Whether the supplied mask is filled. Defaults to False.
        save (bool, optional): Whether to save the plot as an image. Defaults to False.
        save_path (str, optional): The path to save the image. Defaults to None.
    """
    
    figure = plt.figure(figsize=(10, 10), dpi=300)
    
    if ct and not masked_array:
        plt.imshow(ct.array[slice_nr, :, :], cmap='gray',
                   vmin= -1000, vmax=1000)
        
    if ct and masked_array:
        plt.imshow(ct['array'][slice_nr, :, :], cmap='gray',
                   vmin= -1000, vmax=1000)

    if dose and not masked_array:
        dose_alphas = np.where(dose.resampled_array[slice_nr, :, :] < 20, 
                               0, 0.3)
        plt.imshow(dose.resampled_array[slice_nr, :, :], cmap='jet', 
                   alpha=dose_alphas, vmin=20)
        
        cbar = plt.colorbar(shrink=0.5)
        cbar.set_label('Dose (Gy)')
        
    if dose and masked_array:
        dose_alphas = np.where(dose['array'][slice_nr, :, :] < 20, 
                               0, 0.3)
        plt.imshow(dose['array'][slice_nr, :, :], cmap='jet', 
                   alpha=dose_alphas, vmin=20)
        
        cbar = plt.colorbar(shrink=0.5)
        cbar.set_label('Dose (Gy)')

    if mask is not False and mask_fill is True:
        mask_alphas = np.where(mask['array'][slice_nr, :, :] == 0, 0, 0.3)
        plt.imshow(mask['array'][slice_nr, :, :], cmap='Greens',
                   alpha=mask_alphas)

    if mask is not False and mask_fill is False:
        mask_alphas = np.where(mask['array'][slice_nr, :, :] == 0, 0, 0.99)
        plt.imshow(mask['array'][slice_nr, :, :], cmap='Reds',
                   alpha=mask_alphas)

    title = 'Patient: {}, Slice: {}'.format(patient_id, slice_nr)
    plt.title(title)
    plt.axis('off')
    
    if save:
        try:
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            

            filename = f'{patient_id}_{struct_name}_slice_{slice_nr}.png'
            plt.savefig(join(save_path, filename), bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f'Could not save {struct_name} for patient {patient_id} due to {e}')
    
    if not save:
        plt.show()
        plt.close('all')
    
def get_slice_nr(array, position='middle'):
    """
    Returns the slice number based on the position parameter.

    Parameters:
    array (ndarray): The input array.
    position (str): The position of the slice. Can be 'middle', 
    'cranial', or 'caudal'. Defaults to 'middle'.

    Returns:
    int: The slice number.

    """
    if position == 'middle':
        return int(array.shape[0] / 2)
    
    elif position == 'cranial':
        return int(array.shape[0] * 0.75)
    
    elif position == 'caudal':
        return int(array.shape[0] * 0.25)
    
    else:
        raise ValueError('Position must be either "middle", "cranial", or "caudal".')

if __name__ == '__main__':
    pass