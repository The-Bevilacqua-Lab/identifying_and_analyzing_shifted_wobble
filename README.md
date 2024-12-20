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
### 1. step_1_data_preparation.py
This script will read the CSV files (custom search results from RCSB Protein Data Bank (https://www.rcsb.org/)) stored in the 'data' folder  and output the cleaner version of the result within user-defined resolution cut-off. 
```sh
python step_1_data_preparation.py '3.2'
```
### 2. step_2_download_and_characterize.py
This script downloads the structures and characterizes them using DSSR (version v2.2.1-2021jan12 is required for this pipeline). You need to specify the directory where the CSV file (which will be stored in results by default) generated from step_1_data_preparation.py is stored and the directory where the downloaded structure files and JSON output files from DSSR characterization will be saved.
```sh
python step_2_download_and_characterize.py 'results/data.csv' '/directory_for_structures/'
```
### 3. step_3_filter_out_bp.py
This script filters out all the base pairs of interest within the structures downloaded and characterized in step 2. For this, it will be all G-U base pairs. To run the script you need the same CSV file used in step 2 and the directory of the structure and JSON output files. The output will be stored in the 'result' folder. 
```sh
python step_3_filter_out_bp.py 'results/data.csv' '/directory_for_structures/'
```
### 4. step_4_filter_out_registers.py
This script will identify the different registers for G-U wobbles. However, if the user wants to identify registers from other base pairs, the required hydrogen bonds can be specified in lines #186 and 187 of the current version of the script. To run this script you need the CSV file generated in the previous step (step_3_filter_out_bp.py) and the directory of the structure and JSON files. 
```sh
python step_4_filter_out_registers.py 'results/data-from-step-3.csv' '/directory_for_structures/'
```
### 5. step_5_data_extraction_RCSB.py
This script will web-scrape metadata (for instance, Experimental_Method, Resolution, Molecule, Source_organism) for the structures from the RCSB Protein Data Bank dTo run this script, you need to specify the directory containing the CSV file from step 4 and indicate whether you want data for the standard or shifted wobble. Enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_5_data_extraction_RCSB.py 'results/data-from-step-4.csv' '1' or '2'
```
### 6. step_6_data_extraction_structure.py
This script calculates the dihedral angles for hydrogen bonds associated with the specified register of interest, as well as the distances between three heteroatoms on the Watson-Crick-Franklin face of G and three heteroatoms on the Watson-Crick-Franklin face of U. It will also extract the temperature factor for the nucleobase G and U forming the shifted wobble. The script requires the CSV file generated in step 5 and an integer input to specify the register: enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_6_data_extraction_structure.py 'results/data-from-step-5.csv' '1' or '2'
```
### 7. step_7_filter_and_quality_check.py
This script will check the hydrogen bond quality and filter out those with poor quality. The script requires the CSV file generated in step 6 and an integer input to specify the register: enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_7_filter_and_quality_check.py 'results/data-from-step-6.csv' '1' or '2'
```
### 8. step_8_adjust_res_index.py
This script 
### 9. step_9_redundancy_check.py
This script will identify representative structures for redundant wobble examples, where multiple wobble instances from the same location in the same RNA of the same organism are available. The script requires the CSV file generated in step 8 and an integer input to specify the register: enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_9_redundancy_check.py 'results/data-from-step-8.csv' '1' or '2'
```
### 10. step_10_analyze_map_model_cc.py

### 11. step_11_identifying_structural_clusters.py

### 12. step_12_stem_check.py

### 13. step_13_non_WCF_check.py

### 14. step_14_assigning_consensus_1D_2D.ipynb

