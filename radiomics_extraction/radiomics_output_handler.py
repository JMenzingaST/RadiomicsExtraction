"""
This script has two functions to handle the radiomcis output. This could be extended on a 'as needed' basis.
Author: Jacob Menzinga
Version 2.0
"""

# Standard imports
import re
import os
from os.path import join

# External imports
import pandas as pd
def merge_output(folder, exclusion_dose = None):
    """
    Merges all radiomics output files into one dataframe and saves it as an 
    excel file. Excludes patients that have a MLD < 1.
    folder: radiomics, radiomics_LMInf, radiomics_60%bp, radiomics_0%bp"""
    
    pattern = ".*-.{3}_"

    if folder == 'radiomics':
        root = 'radiomics_output'
    elif folder == 'radiomics_LMInf':
        root = 'radiomics_LMInf_output'
    elif folder == 'radiomics_60%bp':
        root = 'radiomics_60%BP_output'
    elif folder == 'radiomics_0%bp':
        root = 'radiomics_0%BP_output'
    else:
        raise ValueError('folder should be radiomics or radiomics_LMInf')
    
    exclusion_file = join(root, 'excluded_patients_in_merging.txt')
    with open(exclusion_file, 'w') as xf:
        xf.write('Patient ID \t Reason for exclusion\n')
    
    radiomics_folder = join(root, 'radiomics')

    p_ids = [file.split('_')[0] for file in os.listdir(radiomics_folder)]
    radiomics_files = [os.path.join(radiomics_folder, file) 
                       for file in os.listdir(radiomics_folder)]
    df = pd.DataFrame()
    
    for file in radiomics_files:
        f = pd.read_csv(file, index_col=0)
        f.columns = [re.sub(pattern, '', col) for col in f.columns]
        
        if exclusion_dose is not None and\
            f[f.columns[-2]].values[0] < exclusion_dose:
            with open(exclusion_file, 'a') as xf:
                xf.write(f'{f.index.values[0]} \t MLD < 1\n')
            continue
        
        else:
            df = pd.concat([df, f], axis=0)
        
    radiomics_diagnostics = [col for col in df.columns if 'diagnostics' in col]
    df.drop(radiomics_diagnostics, axis=1, inplace=True)
        
    df.to_excel(join(root, f'cleaned_{folder}.xlsx'))

def delta_radiomics(radiomics_bp1_path, radiomics_bp2_path, output_path, 
                    absolute_values=False):
    
    radiomics_bp1 = pd.read_excel(radiomics_bp1_path, index_col=0)
    radiomics_bp2 = pd.read_excel(radiomics_bp2_path, index_col=0)
    
    shared_patients = \
        list(set(radiomics_bp1.index).intersection(set(radiomics_bp2.index)))
        
    vals1 = radiomics_bp1.loc[shared_patients].values
    vals2 = radiomics_bp2.loc[shared_patients].values
    delta_vals = vals1 - vals2
    
    if absolute_values:
        delta_vals = abs(delta_vals)
        
    cols = [f"delta_{'_'.join(col.split('_')[-2:])}" 
            for col in radiomics_bp1.columns]
    
    delta_df = pd.DataFrame(delta_vals, index=shared_patients, columns=cols)
    delta_df.to_excel(output_path)
    
if __name__ == '__main__':
    merge_output('radiomics')
    merge_output('radiomics_LMInf')
    merge_output('radiomics_0%bp')
    merge_output('radiomics_60%bp')
    
    bp0_path = 'radiomics_0%BP_output\cleaned_radiomics_0%bp.xlsx'
    bp60_path = 'radiomics_60%BP_output\cleaned_radiomics_60%bp.xlsx'
    
    delta_radiomics(bp0_path, bp60_path, 
                    'data_output\delta_radiomics_0%bp_60%bp.xlsx')
    
    print('Done!') 