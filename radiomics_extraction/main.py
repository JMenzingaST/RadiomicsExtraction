"""
In this script, the radiomics features are extracted from average ct images, 
using only the deliniations found in RTSTRUCT files.

Author: Jacob Menzinga
Version: 2
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
import pandas as pd

with open('config.yaml', 'r') as stream:
    
    config = yaml.safe_load(stream)
    log_folder = config['log_folder']
    output_folder = config['radiomics_output']
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

# Only select patients with relevant scans available    
df = scans[['CT', 'CT_path', 'RTSTRUCT', 'RTSTRUCT_path',
            'RTDOSE', 'RTDOSE_path','60%', '60%_path']
           ][(scans['CT'] == True) &
             (scans['RTSTRUCT'] == True) &
             (scans['RTDOSE'] == True) &
             (scans['60%'] == True)]

# real work.
def main(pat):
    start = time.time()
    try:
        # Load the CT, RTSTRUCT and RTDOSE
        logger.info(f'Starting on patient {pat}')
        logger.info(f'{pat} - Loading CT')
        ct = CTLoader(df.loc[pat, 'CT_path'])
        logger.info(f'{pat} - Loading RTSTRUCT')
        rtstruct = RTSTRUCTLoader(df.loc[pat, 'RTSTRUCT_path'], ct)
        logger.info(f'{pat} - Loading RTDOSE')
        rtdose = RTDOSELoader(df.loc[pat, 'RTDOSE_path'], ct)

        pat = str(pat)
        #making masks
        logger.info(f'{pat} - Making masks')
        
        lungs_found = False
        tumor_found = False
        
        try:
            lung_li = mm.make_mask(rtstruct.contours['lung_left'], ct)
            lung_re = mm.make_mask(rtstruct.contours['lung_right'], ct)
            lungs = mm.manipulate_mask(lung_li, lung_re, 'add')
            lungs_found = True
            
        except KeyError:
            logger.warning(f'{pat} - No individual lung contours found, using combined lung contour instead')
            lungs = mm.make_mask(rtstruct.contours['lungs'], ct)
            lungs_found = True
            
        try:
            tumor = mm.make_mask(rtstruct.contours['gtv'], ct)
            tumor_found = True
        except KeyError:
            logger.warning(f'{pat} - No GTV contour found, using ITV instead')
            tumor = mm.make_mask(rtstruct.contours['itv'], ct)
            tumor_found = True
        
            
        if not lungs_found and not tumor_found:
            try: 
                lungs_minus_tumor = mm.make_mask(rtstruct.contours['gtv_lung'], ct)

            except KeyError:
                pass
                
        else:
            lungs_minus_tumor = mm.manipulate_mask(lungs, tumor, 'subtract')
        
        # Deformable registration of the (lung) mask(s) to different CT Phases.
        logger.info('Deformable registration')
        
        # Manipulating the masks into images.
        # Potentially, I could make a list of masks (lungs, tumor, etc) and then
        # loop over them to create the images.
        
        logger.info('Creating images')
        masked_ct = mm.create_masked_array(ct.array, lungs_minus_tumor)
        masked_dose = mm.create_masked_array(rtdose.resampled_array, 
                                             lungs_minus_tumor)
        
        c_mask = cr.crop_3d_array(lungs_minus_tumor)
        c_masked_ct = cr.crop_3d_array(lungs_minus_tumor, masked_ct)
        c_masked_dose = cr.crop_3d_array(lungs_minus_tumor, masked_dose)
       
        # Visualizing the masks
        slcs = [get_slice_nr(c_mask['array'], position=pos) 
                for pos in ['cranial', 'middle', 'caudal']]
        
        for slc in slcs:
            plot(
                slc, ct.metadata['patient_id'], 
                struct_name=c_masked_ct['name'],
                ct=c_masked_ct, dose=c_masked_dose, masked_array=True,
                save=False, save_path=join(output_folder, 'mask_visualisations') #!
                )
            
        c_mask_img = mm.create_image(c_mask)
        c_masked_ct_img = mm.create_image(c_masked_ct)
        c_masked_dose_img = mm.create_image(c_masked_dose)
        
        mm.save_image(c_mask_img, 'mask', pat, 
                      folder = join(output_folder, 'mask_images'))
        mm.save_image(c_masked_ct_img, 'ct', pat, 
                      folder = join(output_folder, 'mask_images'))
        mm.save_image(c_masked_dose_img, 'dose', pat, 
                      folder = join(output_folder, 'mask_images'))
        
    except Exception as e:
        logger.error(f'Error in patient {pat}:')
        logger.error(traceback.format_exc())
    
if __name__ == '__main__':
    import multiprocessing as mp
    with mp.Pool(6) as pool:
        pool.map(main, df.index)
        pool.close()

    
    
    

