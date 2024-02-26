"""
This script performs the cleaning and processing of a scan database.
It merges structured and unstructured data, removes excluded patients,
converts study dates to datetime format, adds RT start dates to the database,
reorders columns, removes scans performed after the RT start date,
creates a scan availability frame, adds 4D phases to the scan availability frame,
and saves the cleaned scan database and scan availability frame to files.

The script requires a configuration file (config.yaml) that specifies the paths
to the input data files and the output directories.

By: Jacob Menzinga
Date: 14-11-2023

"""
# Standard imports
import os
from os.path import join
import logging
import time
import re

# External imports
import yaml
import pandas as pd

with open(r'\\zkh\appdata\RTDicom\Projectline - Modelling lung cancer outcomes [panama code]\Users\Jacob Menzinga\Scripts\Data pipeline\0_retrieving_and_cleaning_scandata\config.yaml', 'r') as stream:
    config = yaml.safe_load(stream)
    unstructured = pd.read_csv(config['unstructured'])
    structured_pet = pd.read_csv(config['structured_pet'])
    structured_ct = pd.read_csv(config['structured_ct'])
    endpoints = pd.read_excel(config['endpoints'])
    log_output = config['log_output']
    data_output = config['data_output']

logging.basicConfig(filename=join(log_output, 'cleaning_scan_db_' +
    time.strftime("%Y%m%d_%H%M", time.localtime()) + '.log'),
    format='%(levelname)s - %(name)s - %(asctime)s - %(message)s',
    level=logging.INFO)

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

#    ________________________ Merging datasources ________________________    #
# First we need to merge the structured and unstructured data.

logger.info('Merging structured and unstructured data')
df = pd.concat([unstructured, structured_pet, structured_ct],
               axis=0, ignore_index=True)
df.set_index('p_id', inplace=True)

if df.isna().sum().sum() > 0:
    logger.info('Missing values in the database:\n%s', df.isna().sum())
    logger.info('Filling them up with: "name not registered"')
    df.fillna('name not registered', inplace=True)

logger.info('Removing excluded patients')
clean_endpoints = endpoints[endpoints['Flowchart'] == 'inlc']
logger.info('Removed %s patients', len(endpoints) - len(clean_endpoints))
rt_startdate = clean_endpoints.set_index('patient_id')['RT_Start_Date']

#   ______________________ Study date to date time ______________________    #

nodate = df[df['study_date'] == '--']
nodate = df.loc[nodate.index].groupby('p_id')['study_date'].unique()
nodate = nodate.apply(lambda x: x[1])
check_id = nodate.index
nodate = zip(nodate.index, nodate.values)

for i, date in nodate:
    df.loc[i, 'study_date'] = date

df['study_date'] = pd.to_datetime(df['study_date'], format='mixed')

# Patients that are excluded from the endpoints
excluded_patients = \
    df.index.unique()[(~df.index.unique().isin(clean_endpoints['patient_id']))
                      & (df.index.unique().isin(endpoints['patient_id']))]

# Patients that are in the scan database but not in the endpoints
no_endpoint = \
    df.index.unique()[(~df.index.unique().isin(clean_endpoints['patient_id']))
                      & (~df.index.unique().isin(endpoints['patient_id']))]

# Lets safe the patients we're going to exclude for later reference.
with open(join(data_output, 'exclusion_on_endpoint.txt'), 'w+', encoding='utf-8') as doc:
    doc.write(
        f'excluded patients: {len(excluded_patients.tolist())}, {excluded_patients.to_list()} \n')
    doc.write(
        f'patients with no endpoint:{len(no_endpoint.tolist())}, {no_endpoint.tolist()}')

# Now we can remove the excluded patients from the scan database
clean_df = df.drop(excluded_patients, axis=0)
print(len(clean_df))
clean_df = clean_df.drop(no_endpoint, axis=0)
logger.info('Removed patients without endpoints from the scan database')
logger.info('The database now contains %s patients',
            len(clean_df.index.unique()))

# adding RT_startdate to the scan database
clean_df['RT_startdate'] = rt_startdate.loc[clean_df.index]

# reordering the columns for better overview
clean_df = clean_df[['modality', 'study_date', 'RT_startdate',
                     'study_description', 'series_description',
                     'folder_path']]

# Removing scans that have a later study date than the RT start date
clean_df.reset_index(inplace=True)
scan_after_rtstart = clean_df[clean_df['study_date']
                              > clean_df['RT_startdate']]
scan_after_rtstart.to_excel(join(data_output, 'scans_after_rtstart.xlsx'))
logger.info('Saved scans that were performed after the RT start date to %s',
            join(data_output, 'scans_after_rtstart.xlsx'))

clean_df.drop(scan_after_rtstart.index, inplace=True)
# clean_df.set_index('p_id', inplace=True)
clean_df.to_csv(
    join(data_output, 'cleaned_collection_of_scans_per_patient.csv'))

logger.info('Saved cleaned scan database to %s',
            join(data_output, 'cleaned_collection_of_scans_per_patient.csv'))

#  ______________________ Making availability sheet ______________________    #

scantype_dict = {'CT': 'Thorax|Ave|MIP', '4D-phases': '\d{1,3}%',
                 'PET': 'PET|FDG', 'MRI': "MR|MRI",
                 'RTDOSE': 'DOSE', 'RTSTRUCT': 'STRUCT'}

scan_available = pd.DataFrame(index=clean_df['p_id'].unique())

logger.info('Creating scan availability frame')
for _, row in clean_df.iterrows():
    for key, value in scantype_dict.items():
        if key == '4D-phases':
            if re.search(value, row['series_description']) and \
                    not re.search(scantype_dict['CT'], row['series_description']):
                scan_available.loc[row['p_id'], key] = True
                scan_available.loc[row['p_id'], key +
                                   '_path'] = row['folder_path']

        if key == 'RTDOSE':
            if re.search(value, row['modality']) \
                and not re.search('VMAT', row['series_description']):
                
                scan_available.loc[row['p_id'], key] = True
                scan_available.loc[row['p_id'], key +
                                   '_path'] = row['folder_path']
        
        else:
            if re.search(value, row['series_description']):
                scan_available.loc[row['p_id'], key] = True
                scan_available.loc[row['p_id'], key +
                                   '_path'] = row['folder_path']
            elif re.search(value, row['modality']):
                scan_available.loc[row['p_id'], key] = True
                scan_available.loc[row['p_id'], key +
                                   '_path'] = row['folder_path']

# ______________________ Adding 4D phases ______________________    #

logger.info('Adding 4D phases to the scan availability frame')
df_4d = clean_df[(clean_df['modality'] == 'CT') &
                 (clean_df['series_description'].str.contains('\d{1,2}%')) &
                 ~(clean_df['series_description'].str.contains('Ave|MIP'))]

# Rounding phases to 10% intervals
# First some support funcitons


def get_base_phase(row):
    row = row.split('  ')
    base = ' '.join(row[:-1])
    phase = ' '.join(row[-1:])

    return base, phase


def round_to_reference_values(value, reference_values):
    # Calculate the absolute differences between the value and each reference value
    differences = [abs(value - ref) for ref in reference_values]
    # Find the index of the reference value with the smallest difference
    min_difference_index = differences.index(min(differences))

    # Return the reference value at the index with the smallest difference
    return reference_values[min_difference_index]

# Proplem: The phases are a number followed by string that both are important to keep.
# I'll have to create a function that splits the numvers from the phase indicator
# Then rounds the number to the nearest of the top 8 phases.
# Finally it should merge the number and phase indicator back together.


def convert_phases(row, reference_values):
    num, phase = row.split('%')
    num = round_to_reference_values(int(num), reference_values)
    return str(num) + '% ' + phase


phase_binned = df_4d[df_4d['series_description'].str.contains(
    'RespLow|dal|top')]
phase_binned['split'] = phase_binned['series_description'].apply(
    get_base_phase)

phase_binned['base'], phase_binned['phase'] = zip(*phase_binned['split'])
phase_binned.drop(phase_binned[phase_binned['phase'].str.contains('iMAR')
                               ].index, inplace=True)

# Now to round the phase numbers to the nearest decimal.
reference_vals = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
phase_binned['converted_phase'] = phase_binned['phase'].apply(
    convert_phases, reference_values=reference_vals)

# And now lets update the series description to the new rounded phases
phase_binned['series_description'] = phase_binned['base'] + \
    ' ' + phase_binned['converted_phase']
phase_binned.drop(
    ['split', 'base', 'phase', 'converted_phase'], axis=1, inplace=True)

# Now they need to get merged back into the full 4D-DF
df_4d.loc[phase_binned.index, 'series_description'] = \
    phase_binned['series_description']

df_4d.set_index('p_id', inplace=True)

# A dict with filters to map the 4D phases to the right phase percentage
phase_fltrs = {'0%': ['T=0%', 'TRIGGER_DELAY 0%', r'100%\s+Ex', 'Qr40 S3 0%'],
               '10%': ['T=10%', 'TRIGGER_DELAY 10%', r'90%\s+Ex', 'Qr40 S3 10%'],
               '20%': ['T=20%', 'TRIGGER_DELAY 20%', r'60%\s+Ex', 'Qr40 S3 20%'],
               # No matching RespLow phase
               '30%': ['T=30%', 'TRIGGER_DELAY 30%', 'Qr40 S3 30%'],
               '40%': ['T=40%', 'TRIGGER_DELAY 40%', r'30%\s+Ex', 'Qr40 S3 40%'],
               '50%': ['T=50%', 'TRIGGER_DELAY 50%', r'5%\s+Ex', 'Qr40 S3 50%'],
               '60%': ['T=60%', 'TRIGGER_DELAY 60%', r'0%\s+Ex|0% In', 'Qr40 S3 60%'],
               '70%': ['T=70%', 'TRIGGER_DELAY 70%', r'40%\s+In', 'Qr40 S3 70%'],
               '80%': ['T=80%', 'TRIGGER_DELAY 80%', r'60%\s+In', 'Qr40 S3 80%'],
               '90%': ['T=90%', 'TRIGGER_DELAY 90%', r'90%\s+In', 'Qr40 S3 90%'],
               '100%': ['T=100%', 'TRIGGER_DELAY 100%', r'100%\s+In', 'Qr40 S3 100%']}

phases = pd.DataFrame(index=df_4d.index.unique())
for index, row in df_4d.iterrows():
    for key, value in phase_fltrs.items():
        for fltr in value:
            if re.search(fltr, row['series_description']):
                phases.loc[index, key] = True
                # adding path, might want to dissable this
                phases.loc[index, key+'_path'] = row['folder_path']

scans_and_phases_available = pd.merge(scan_available, phases, how='left',
                                      left_index=True, right_index=True)

scans_and_phases_available.fillna(False, inplace=True)

scans_and_phases_available.to_excel(join(data_output,
                                         'scan_and_phase_availability.xlsx'))

logger.info('Saved scan and phase availability to %s',
            join(data_output, 'scan_and_phase_availability.xlsx'))