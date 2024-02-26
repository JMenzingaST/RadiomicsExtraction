"""
This script performs deformable registration of a tumor mask and radiomics
extraction from the lungs.

Author: Jacob Menzinga
Version: 1.1
"""

# Standard imports
import yaml
import traceback
import os
from os.path import join
import logging
import time

# External imports
from Loaders import CTLoader, RTSTRUCTLoader, RTDOSELoader
from MaskMaker import MaskMaker as mm
from Cropper import Cropper as cr
from Plotter import plot, get_slice_nr
from RadiomicsExtractor import RadiomicsExtractor
from Registrator import Registrator
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import SimpleITK as sitk

### Setting up the nececarry parts ###
with open('config.yaml', 'r') as stream:
    
    config = yaml.safe_load(stream)
    log_folder = config['log_folder']
    output_folder = config['radiomics_0%BP_output']
    scans = pd.read_excel(config['scan_availability']
                          ).drop(['MRI', 'MRI_path'], axis=1
                                 ).rename(columns={'Unnamed: 0': 'p_id'}
                                          ).set_index('p_id')
                                 
logging.basicConfig(filename=join(log_folder, 'radiomics_extraction_log_'+ 
                                  f'{time.strftime("%Y%m%d_%H%M",time.localtime())}.log'),
                    format='%(levelname)s - %(name)s - %(asctime)s - %(message)s',
                    level=logging.INFO)

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

### Selecting patients with relevant scans available ###
df = scans[['CT', 'CT_path', 'RTSTRUCT', 'RTSTRUCT_path',
            'RTDOSE', 'RTDOSE_path','60%', '60%_path', '0%', '0%_path']
           ][(scans['CT'] == True) &
             (scans['RTSTRUCT'] == True) &
             (scans['RTDOSE'] == True) &
             (scans['60%'] == True) &
             (scans['0%'] == True)]


### main ###
def main(pat):
    start = time.time()
    try:
        ### Load the CT, RTSTRUCT and RTDOSE ###
        logger.info(f'Starting on patient {pat}')

        bp0 = CTLoader(df.loc[pat, '0%_path'])
        bp60 = CTLoader(df.loc[pat, '60%_path'])
        rtstruct = RTSTRUCTLoader(df.loc[pat, 'RTSTRUCT_path'], bp60)
        rtdose = RTDOSELoader(df.loc[pat, 'RTDOSE_path'], bp60)

        pat = str(pat)
        ### making masks ###
        logger.info(f'{pat} - Making tumor mask image')
        tumor_found = False
        try:
            tumor = mm.make_mask(rtstruct.contours['gtv'], bp60)
            tumor_found = True
        except KeyError:
            logger.warning(f'{pat} - No GTV contour found, using ITV instead')
            tumor = mm.make_mask(rtstruct.contours['itv'], bp60)
            tumor_found = True 
            
        if not tumor_found:
            logger.error(f'{pat} - No tumor contour found, skipping patient | early exit to avoid waiting for registration')
            return None
        
        ### Deformable registration of the tumor mask to different CT Phases. ###
        logger.info(f'{pat} - starting the registration, this will take ~15 minutes')
        rs = Registrator(bp0, bp60, tumor, rtdose)
        f = rs.plot_registration(len(bp0.array)//2)
        f.savefig(join(output_folder, 'registration_visualisations', f'{pat}_registration.png'))
        
        # Manipulating the masks into images.        
        # First we need to get the arrays, all the same size. 
        # Admittedly this is pretty conveluted but I dont have the time, 
        # nor mental space to rewrite the classes to faccilitate this easier
        
        c_bp0_lung_mask = {'name':'bp0_lung_mask', 
                           'array':rs.c_fixed_lungmask_arr}
        c_m_bp0 = mm.create_masked_array(rs.c_fixed_ct_arr, c_bp0_lung_mask)
        c_m_res_dose = mm.create_masked_array(rs.c_rtdose_arr, c_bp0_lung_mask)
        
        # !!!
        # !!! Somewhere, the mask array is not all nan's so the average dose is 
        # !!! Not correct. This gets fixed in 'non_nan_mask_fix.py'
        # !!!
        
        ### Visualizing the masks ###
        slcs = [get_slice_nr(c_bp0_lung_mask['array'], position=pos) 
                for pos in ['cranial', 'middle', 'caudal']]
        
        for slc in slcs:
            plot(
                slc, bp0.metadata['patient_id'], 
                struct_name=c_m_bp0['name'],
                ct=c_m_bp0, dose=c_m_res_dose, masked_array=True,
                save=True, save_path=join(output_folder, 'mask_visualisations') #!
                )
            
        c_mask_img = mm.create_image(c_bp0_lung_mask)
        c_m_bp0_img = mm.create_image(c_m_bp0)
        c_m_res_dose_img = mm.create_image(c_m_res_dose)
        
        mm.save_image(c_mask_img, 'mask', pat, 
                      folder = join(output_folder, 'mask_images'))
        mm.save_image(c_m_bp0_img, 'ct', pat, 
                      folder = join(output_folder, 'mask_images'))
        mm.save_image(c_m_res_dose_img, 'dose', pat, 
                      folder = join(output_folder, 'mask_images'))
        
        ### Radiomics Extraction ###
        logger.info(f'{pat} - Extracting radiomics')
        
        rx = RadiomicsExtractor(pat)
        rx.extract_features(c_m_bp0_img,
                            c_m_res_dose,
                            c_mask_img)
        
        if not os.path.exists(join(output_folder, 'radiomics')):
            os.mkdir(join(output_folder, 'radiomics'))
        
        rx.features_df.to_csv(join(output_folder, 'radiomics', 
                                   f'{pat}_bp0_radiomics.csv'))
        
        ### Ending ###
        time_elapsed = time.time() - start
        logger.info(f'{pat} - Finished in {time_elapsed:.1f} seconds')
        
    except Exception as e:
        logger.error(f'Error in patient {pat}:')
        logger.error(traceback.format_exc())
    
if __name__ == '__main__':
    # You could multiprocess this if you want, same to the 'main.py' script
    # The PC I developed this on did not have the memory to do so.
    for pat in df.index:
        main(pat)
    
    
    

