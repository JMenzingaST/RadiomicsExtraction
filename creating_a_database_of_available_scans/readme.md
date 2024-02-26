In this folder there's a few scripts that work together to create a Excel sheet that shows the availabillity of various imaging modalities and their path.

E.g.:
| CT    | CT_path   | PET   | PET_path  |
|------ |-----------|-------|---------- |
| True  | C:\Files  | False |   None    |

First edit the config.yaml file to refer to the relevant directories, then run 'retrieve_scan_info_main' script. Lastly, run the 'cleaning_scan_db' script.
