"""
This script contains the implementation of the Registrator class
The Registrator class is used to perform deformable registration of a tumor mask
Thanks to Luis de la Arevalo for providing scripts and notebooks that served as
the basis and inspiration for this class.

Author: Jacob Menzinga
Version 2.0
"""

# Standard imports 
from os.path import join
import logging

# External imports
import SimpleITK as sitk
from lungmask import LMInferer
import numpy as np
import matplotlib.pyplot as plt
from Loaders import CTLoader, RTSTRUCTLoader, RTDOSELoader

# setting up a logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

class Registrator:
    """
    Class for performing registration of CT images and masks.
    
    Args:
        - fixed_ct (CTLoader): The fixed CT image.
        - moving_ct (CTLoader): The moving CT image.
        - tumor_mask (sitk.Image): The tumor mask.
        - rtdose (dict): The RT dose.
    """
    def __init__(self, fixed_ct:CTLoader, moving_ct:CTLoader,
                 tumor_mask:sitk.Image, rtdose:dict) -> None:
        
        self.fixed_ct = fixed_ct
        self.moving_ct = moving_ct
        self.tumor_mask = tumor_mask
        self.rtdose = rtdose
        
        self.fixed_lung_mask = self.gen_lung_mask(fixed_ct.image)
        self.moving_lung_mask = self.gen_lung_mask(moving_ct.image)
        
        self.min_indices, self.max_indices = \
            self.get_bounds(self.fixed_lung_mask)
        
        self.c_fixed_ct_arr = self.crop_array(self.min_indices, 
                                              self.max_indices, 
                                              fixed_ct.array)
        
        self.c_fixed_lungmask_arr = self.crop_array(self.min_indices, 
                                                    self.max_indices,
                                                    self.fixed_lung_mask)
        
        self.c_moving_ct_arr = self.crop_array(self.min_indices,
                                               self.max_indices,
                                               self.moving_ct.array)
        
        self.c_moving_lungmask_arr = self.crop_array(self.min_indices,
                                                     self.max_indices,
                                                     self.moving_lung_mask)
        
        self.c_tumor_mask_arr = self.crop_array(self.min_indices,
                                                self.max_indices,
                                                tumor_mask['array'])
        
        self.c_rtdose_arr = self.crop_array(self.min_indices,
                                            self.max_indices,
                                            rtdose.resampled_array)
        
        self.c_reg_tumor_mask_img = self.deformably_register(
            self.create_image(self.c_fixed_ct_arr, moving_ct.image),
            self.create_image(self.c_moving_ct_arr, moving_ct.image),
            self.create_image(self.c_tumor_mask_arr, moving_ct.image))
        
        self.c_reg_tumor_mask_arr = self.reg_img_to_array(self.c_reg_tumor_mask_img)
         

    def gen_lung_mask(self, ct_image):
        """
        Generates lungmask from a CT image. using the lungmask package.
        
        Args:
            - ct_image (sitk.Image): The CT image.
            
        Returns:
            - lungmask (numpy.ndarray): The lungmask.
        """
        inferer = LMInferer()
        lungmask = inferer.apply(ct_image)
        lungmask = lungmask.clip(0,1)
        return lungmask

    def get_bounds(self, lungmask, padding = 0):
        """
        Get the cropping bounds for a given mask.
        
        Args:
            - mask (ndarray): The mask array.
            - padding (int): The padding to add to the cropping bounds.
            
        Returns:
            - min_indices (list): A list containing the minimum values for each 
              axis (x_min, y_min, z_min).
            - max_indices (list): A list containing the maximum values for each 
              axis (x_max, y_max, z_max).
        """
        mask_bb = np.where(lungmask>0)
        min_indices = [np.min(bb) for bb in mask_bb]
        max_indices = [np.max(bb) for bb in mask_bb]
        
        min_indices = [max(0, idx - padding) for idx in min_indices]
        max_indices = [min(dim_size - 1, idx + padding) for idx, dim_size in 
                       zip(max_indices, lungmask.shape)]
        
        return min_indices, max_indices

    def crop_array(self, min_indices, max_indices, array):
        """
        Crop a 3D array based on the provided mask and return the cropped array.
        
        Args:
            - min_indices (list): A list containing the minimum values for each
                axis (x_min, y_min, z_min).
            - max_indices (list): A list containing the maximum values for each
                axis (x_max, y_max, z_max).
            - array (numpy.ndarray): The array to crop.
            
        Returns:
            - cropped_array (numpy.ndarray): The cropped array.
        """
        
        cropped_array = array[min_indices[0]:max_indices[0]+1, 
                            min_indices[1]:max_indices[1]+1, 
                            min_indices[2]:max_indices[2]+1]

        return cropped_array

    def create_image(self, array, ref_ct_img):
        """
        Create a SimpleITK image from a numpy array and set its metadata based on a reference CT image.

        Args:
            - array (numpy.ndarray): The input array representing the image.
            - ref_ct_img (SimpleITK.Image): The reference CT image to set the metadata.

        Returns:
            SimpleITK.Image: The created image with metadata set.
        """
        image = sitk.GetImageFromArray(array)
        image.SetOrigin(ref_ct_img.GetOrigin())
        image.SetDirection(ref_ct_img.GetDirection())
        image.SetSpacing(ref_ct_img.GetSpacing())

        return image

    def b_spline(self, fixed_image, moving_image):
        """
        Performs B-spline registration between a fixed image and a moving image.
        
        Args:
        - fixed_image: The fixed image to register.
        - moving_image: The moving image to register.
        
        Returns:
        - final_transform: The final transformation obtained from the 
        """
                     
        control_point = 3
        order = 1
        
        registration_method = sitk.ImageRegistrationMethod()
        #Metrics
        registration_method.SetMetricAsCorrelation()
        #registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
        registration_method.SetMetricSamplingPercentage(0.3)
        #Optimizer
        registration_method.SetOptimizerAsRegularStepGradientDescent(
            learningRate=0.1, minStep=1e-4, numberOfIterations=200)
        
        bspline_transform = sitk.BSplineTransformInitializer(
            fixed_image, 
            transformDomainMeshSize=[control_point]*fixed_image.GetDimension(), 
            order=order)
        
        registration_method.SetInitialTransform(bspline_transform, 
                                                inPlace=False)
        
        registration_method.SetShrinkFactorsPerLevel([4, 2, 1])
        registration_method.SetSmoothingSigmasPerLevel([2, 1, 0])

        # Perform the registration
        final_transform = registration_method.Execute(fixed_image, 
                                                        moving_image)
        
        evaluation_metric = registration_method.GetMetricValue()
        initial_metric_value = registration_method.MetricEvaluate(fixed_image, 
                                                                    moving_image)
        
        print(f"Initial metric value: {initial_metric_value}")
        print(f"Final metric value: {evaluation_metric}")

        return final_transform

    def deformably_register(self, fixed_lung_mask_img, 
                            moving_lung_mask_img, 
                            mask_to_deform_img):
        
        """
        Deformably registers the mask_to_deform_img to the fixed_lung_mask_img using the moving_lung_mask_img as a reference.
        
        Parameters:
            - fixed_lung_mask_img (SimpleITK.Image): The fixed lung mask image.
            - moving_lung_mask_img (SimpleITK.Image): The moving lung mask image.
            - mask_to_deform_img (SimpleITK.Image): The image to be deformed and registered.

        Returns:
            SimpleITK.Image: The registered image after deformation.
        """
        
        logger.info('Deformably registering the tumor mask to the fixed CT phase')
        
        fixed_lung_mask_img = sitk.Cast(fixed_lung_mask_img, sitk.sitkFloat32)
        moving_lung_mask_img = sitk.Cast(moving_lung_mask_img, sitk.sitkFloat32)
        mask_to_deform_img = sitk.Cast(mask_to_deform_img, sitk.sitkInt16)
        
        # Perform the registration
        # fixed image and moving image here are the masks
        final_transform = self.b_spline(fixed_lung_mask_img, moving_lung_mask_img)

        registered_image = sitk.Resample(mask_to_deform_img, fixed_lung_mask_img, 
                                            final_transform,
                                            sitk.sitkNearestNeighbor, 
                                            0.0, mask_to_deform_img.GetPixelID())

        return registered_image
    
    def reg_img_to_array(self, reg_img):
        """
        turns the registered image into an array and sets the 0's to np.nan.

        Args:
            reg_img (sitk.Image): the registered image.

        Returns:
            reg_arr (np.ndarray): the registered array.
        """
        reg_arr = sitk.GetArrayFromImage(reg_img)
        reg_arr = reg_arr.astype(float)
        reg_arr[reg_arr == 0] = np.nan
        return reg_arr
    
    def plot_registration(self, slice_nr):
        """
        Plots the registration resykts of the tumor mask to the fixed CT phase.
        
        Args:
            - slice_nr (int): The slice number to plot.
        """
        fig, (ax1, ax2) = plt.subplots(1,2)
        fig.set_size_inches(8, 4)
        fig.set_dpi(300)
        
        ax1.axis('off')
        ax1.imshow(self.c_moving_ct_arr[slice_nr,:,:], cmap='gray')
        ax1.contour(self.c_moving_lungmask_arr[slice_nr,:,:], colors='r')
        ax1.imshow(self.c_rtdose_arr[slice_nr,:,:], alpha=0.2, cmap='jet')
        ax1.imshow(self.c_tumor_mask_arr[slice_nr,:,:], alpha=0.5)

        ax2.axis('off')
        ax2.imshow(self.c_fixed_ct_arr[slice_nr,:,:], cmap='gray')
        ax2.contour(self.c_moving_lungmask_arr[slice_nr,:,:], colors='r')
        ax2.contour(self.c_fixed_lungmask_arr[slice_nr,:,:], colors='b')
        ax2.imshow(self.c_rtdose_arr[slice_nr,:,:], alpha=0.2, cmap='jet')
        ax2.imshow(self.c_tumor_mask_arr[slice_nr,:,:], alpha=0.3)
        ax2.imshow(self.c_reg_tumor_mask_arr[slice_nr,:,:], alpha=0.3, cmap='autumn')
        
        return fig

if __name__ == '__main__':
    pass