"""
This script performs radiomics extraction on lung cancer patients.
It loads CT, RTSTRUCT, and RTDOSE data, creates masks, crops the data, visualizes the masks,
and extracts radiomics features. The extracted features are saved as CSV files.

Author: Jacob Menzinga
Date: 26-01-2024
Version: 1.2
"""

# Standard imports
import os
from os.path import join
import logging
import time
import yaml
import traceback

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
    output_folder = config['radiomics_LMInf_output']
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
        logger.info(f'{str(pat)} - Loading CT')
        ct = CTLoader(df.loc[pat, 'CT_path'])
        logger.info(f'{str(pat)} - Loading RTSTRUCT')
        rtstruct = RTSTRUCTLoader(df.loc[pat, 'RTSTRUCT_path'], ct)
        logger.info(f'{str(pat)} - Loading RTDOSE')
        rtdose = RTDOSELoader(df.loc[pat, 'RTDOSE_path'], ct)

        pat = str(pat)
        #making masks
        logger.info(f'{pat} - Making masks')
        
        try:
            tumor = mm.make_mask(rtstruct.contours['gtv'], ct)
        except KeyError:
            logger.warning(f'{pat} - No ITV contour found, using GTV instead')
            tumor = mm.make_mask(rtstruct.contours['itv'], ct)
            
        lung_li, lung_re = mm.generate_lungmask(ct)
        lungs = mm.manipulate_mask(lung_li, lung_re, 'add')
        lungs_minus_tumor = mm.manipulate_mask(lungs, tumor, 'subtract') 
        
        logger.info(f'{pat} - Creating images')
        masked_ct = mm.create_masked_array(ct.array, lungs_minus_tumor)
        masked_dose = mm.create_masked_array(rtdose.resampled_array, lungs_minus_tumor)
        
        c_mask = cr.crop_3d_array(lungs_minus_tumor)
        c_masked_ct = cr.crop_3d_array(lungs_minus_tumor, masked_ct)
        c_masked_dose = cr.crop_3d_array(lungs_minus_tumor, masked_dose)
       
        # Visualizing the masks
        slcs = [get_slice_nr(c_mask['array'], position=pos) for pos in ['cranial', 'middle', 'caudal']]
        for slc in slcs:
            plot(slc, ct.metadata['patient_id'], 
                 struct_name=c_masked_ct['name'],
                 ct=c_masked_ct, dose=c_masked_dose, masked_array=True,
                 save=True, save_path=join(output_folder, 'mask_visualisations')
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
        
        # Radiomics Extraction
        logger.info(f'{pat} - Extracting radiomics')
        
        radiomics_extractor = RadiomicsExtractor(pat)
        radiomics_extractor.extract_features(c_masked_ct_img, 
                                             c_masked_dose,
                                             c_mask_img)
        
        if not os.path.exists(join(output_folder, 'radiomics')):
            os.mkdir(join(output_folder, 'radiomics'))
        
        radiomics_extractor.features_df.to_csv(join(output_folder, 
                                                    'radiomics', 
                                                    f'{pat}_radiomics.csv'))
        
        # Ending
        time_elapsed = time.time() - start
        logger.info(f'Finished patient {pat} in {time_elapsed:.1f} seconds')
        
    except Exception as e:
        logger.error(f'Error in patient {pat}:')
        logger.error(traceback.format_exc())
    
if __name__ == '__main__':
    # You could multiprocess this if you want, same to the 'main.py' script
    # The PC I developed this on did not have the memory to do so.
    for pat in df.index:
        main(pat)
    
    
    

