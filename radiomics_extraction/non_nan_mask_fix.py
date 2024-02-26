"""
This script fixes an issue I had with the deformable registered mask arrays
having 0 values instead of NaNs. 

Author: Jacob Menzinga
Date: 23-01-2024
Version: 1.0
"""

# Standard Imports
import os
from os.path import join

# External Imports
import pandas as pd
import numpy as np
import SimpleITK as sitk
from natsort import natsorted
import matplotlib.pyplot as plt
from RadiomicsExtractor import RadiomicsExtractor

root = r'.\radiomics_0%BP_output'
maskpath = join(root, 'mask_images')
imgs = os.listdir(maskpath)
pat_nums = [i.split('_')[0] for i in imgs]
pat_nums = natsorted(list(set(pat_nums)))
i = 1
df = pd.DataFrame()

for num in pat_nums[:1]:
    print(num)
    print(f'{i} / {len(pat_nums)}')
    masks = [join(maskpath, img) for img in imgs if num in img]
    
    if len(masks) != 3:
        print('error')
        
    mask = [img for img in masks if '_mask_' in img][0]
    dose = [img for img in masks if '_dose_' in img][0]
    ct = [img for img in masks if '_ct_' in img][0]
    
    mask = sitk.ReadImage(mask)
    dose = sitk.ReadImage(dose)
    ct = sitk.ReadImage(ct)

    mask = sitk.GetArrayFromImage(mask)
    dose = sitk.GetArrayFromImage(dose)
    ct = sitk.GetArrayFromImage(ct)

    ct = ct.astype(float)
    ct[ct == 0] = np.nan
    dose[dose == 0] = np.nan
    
    ct = {'name': 'bp0_lung_minus_tumor', 'image':sitk.GetImageFromArray(ct)}
    mask = {'name': 'bp0_lung_minus_tumor', 'image':sitk.GetImageFromArray(mask)}
    dose = {'name': 'bp0_lung_minus_tumor', 'array': dose}
    try:
        rx = RadiomicsExtractor(num)
        rxdf = rx.extract_features(ct, dose, mask, True)
        rxdf.to_csv(join(root, 'radiomics', f'{num}_bp0_radiomics.csv'))
        i += 1

    except:
        print(f'error with {num}')
        i += 1
        continue
