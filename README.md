# Identifying and analyzing shifted wobble
This pipeline will take the search output from the RCSB protein databank to identify and characterize different orientations of wobbles. 
## Setup instructions
- Install Python >= 3.8
- Install required libraries:
```sh
pip install -r requirements.txt
```
- For characterizing the 3D structure Dissecting the Spatial Structure of RNA (DSSR, version: v2.2.1-2021jan12) was used
## Folder structure

    identify_and_characterize_shifted_wobble
      ├── data                                                 #Search output from RCSB protein data bank
      ├── results                                              #Output files from each python script
      ├── scripts                                              #All python scripts and jupyter notebooks
      ├── requirements.txt                                     #Required python packages and libraries
      ├── README.md                                            #This file
      ├── LICENSE
## Instructions for each scripts

