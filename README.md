# RadiomicsExtraction
Lung cancer is the second most common cancer worldwide, and radiotherapy (RT) is a common treatment for non operable lung cancer. However, radiotherapy can cause side effects which include radiation pneumonitis. Normal tissue complication probability (NTCP) models can be used to estimate the risk of RP and in turn to create a personalised treatment plan. Current NTCP models however do leave room for improvement. One potential approach to improve these models is by including radiomics features extracted from a pre-RT CT-scan.

In broad strokes, the code in the repo does three things: create a database of available medical images and related files, extract radiomics and create and evaluate NTCP models.

* The scripts to create a database of available medical imaging files are in the 'creating_a_database_of_available_scans' folder.
*  A pipeline for radiomics extraction from these scans is found in the 'radiomics_extraction' folder.
*  Lastly, the R scripts for making the NTCP models are found in the 'NTCP_modeling' folder.

See the below image for a visual walk-through the radiomics extraction process

![Visualisation of the radiomics extrection process](https://github.com/JMenzingaST/RadiomicsExtraction/blob/main/RE_visualisation.png?raw=true "Visualisation of the Radiomics Extraction process")

Each of the folders contais it's own readme doc providing relevant information.

## Set up virtual environment

To run the scripts provided in this repo, you need python version 3.7.10 installed, together with the libraries proveided in 'requirements.txt'.

Assuming you have a conda version of python installed on your PC, you can create a virtual enviornment using:

    conda create -n myenv python=3.7.10

After you have created your environment and activated it, please install the requirements in this environment with the command:

    pip install -r /path/to/requirements.txt

