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
1. step_1_data_preparation.py
2. step_2_download_and_characterize.py
3. step_3_filter_out_bp.py
4. step_4_filter_out_registers.py
5. step_5_data_extraction_RCSB.py
6. step_6_data_extraction_structure.py
7. step_7_filter_and_quality_check.py
8. step_8_adjust_res_index.py
9. step_9_redundancy_check.py
10. step_10_analyze_map_model_cc.py
11. step_11_identifying_structural_clusters.py
12. step_12_stem_check.py
13. step_13_non_WCF_check.py
14. step_14_assigning_consensus_1D_2D.ipynb

