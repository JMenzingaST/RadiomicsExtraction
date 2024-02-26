"""
A script that looks through a path and creates an overview of all the patients
and the scans that were performed on them.

By: Jacob Menzigna
Version 2
"""

# standard imports
import os
import logging
from datetime import datetime

# external imports
import yaml
from ScanInfoCrawler import ScanCrawler


# Load config file
with open('config.yaml', 'r', encoding='utf-8') as stream:
    config = yaml.safe_load(stream)
    log_output = config['log_output']
    if log_output is not None:
        log_output = log_output + datetime.now().strftime('%Y-%m-%d_%H%M') + '.log'
    data_dir = config['data_dir']
    csv_output = config['csv_output']
    patient_excel_output = config['patient_excel_output']


# Configure logger
logging.basicConfig(filename=log_output,
                    format='%(levelname)s - %(name)s - %(asctime)s - %(message)s',
                    level=logging.INFO)

# Add logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

def main():
    # Get all the patient numbers in the unstructured data folder
    crawler = ScanCrawler(data_dir, csv_output)
    
    all_patients = crawler.get_patient_list()
    logger.info('Found %s patients in %s', len(all_patients),
                data_dir)

    logger.info('writing header to csv file')
    crawler.csv_writer(crawler.header, 'w+')

    logger.info('starting collection of scan info')
    for patient in all_patients:
        crawler.get_scans_info(patient)

    logger.info('creating excel file')
    ScanCrawler.create_excel(csv_output, patient_excel_output)
    logger.info('all done. Excel file saved to %s', patient_excel_output)


if __name__ == '__main__':
    main()
