# Identifying and analyzing shifted wobble
This pipeline will take the search output from the RCSB protein databank to identify and characterize different orientations of wobbles. 
## Setup instructions
- Install Python >= 3.8
- Install required libraries:
```sh
pip install -r requirements.txt
```
- For characterizing the 3D structure, Dissecting the Spatial Structure of RNA (DSSR, version: v2.2.1-2021jan12) was used
- The Phenix software package (version: 1.21.1-5286) was used to compare electron density maps and modeled structures.  
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
This script will read the CSV files (custom search results from [RCSB Protein Data Bank](https://www.rcsb.org/)) stored in the 'data' folder  and output the cleaner version of the result within user-defined resolution cut-off. 
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
This script will web-scrape metadata (for instance, Experimental_Method, Resolution, Molecule, Source_organism) for the structures from the RCSB Protein Data Bank. To run this script, you need to specify the directory containing the CSV file from step 4 and indicate whether you want data for the standard or shifted wobble. Enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_5_data_extraction_RCSB.py 'results/data-from-step-4.csv' '1' or '2'
```
### 6. step_6_data_extraction_structure.py
This script calculates the dihedral angles for hydrogen bonds associated with the specified register of interest, as well as the distances between three heteroatoms on the Watson-Crick-Franklin face of G and three heteroatoms on the Watson-Crick-Franklin face of U. It will also extract the average temperature factor for the nucleobase G and U forming the shifted wobble. The script requires the CSV file generated in step 5 and an integer input to specify the register: enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_6_data_extraction_structure.py 'results/data-from-step-5.csv' '1' or '2'
```
### 7. step_7_filter_and_quality_check.py
This script will check the hydrogen bond quality and filter out those with poor quality. The script requires the CSV file generated in step 6 and an integer input to specify the register: enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_7_filter_and_quality_check.py 'results/data-from-step-6.csv' '1' or '2'
```
### 8. step_8_get_fastas.py
This script takes the data from step 7 and pulls fasta formatted sequences from pymol and the pdb server, as well as makes a reference fasta with one representative sequence per RNA type and organism.

```sh
python step_8_get_fastas.py -i 'results/data-from-step-7.csv' -ft pymol -o pymol.fasta
python step_8_get_fastas.py -i 'results/data-from-step-7.csv' -ft pdb -o pdb.fasta
python step_8_get_fastas.py -i 'results/data-from-step-7.csv' -ft ref -pdb pdb.fasta -o reference.fasta
```

### 9. step_9_adjust_res_index.py
This script takes the csv from step 7 along with three fasta files derived from step_8_get_fastas.py. The ouput is the same csv from step 7 with additional columns for: ref_org, pdb_adj_res1, pdb_adj_res2, flank1, flank2, ref_pdb_chain, adj_res1, and adj_res2

```sh
python step_9_adjust_res_index.py -i 'results/data-from-step-7.csv' -py pymol.fasta -pdb pdb.fasta -ref reference.fasta -o 'results/data-from-step-9.csv
```

### 10. step_10_redundancy_check.py
This script will identify representative structures from redundant wobble examples, where multiple wobble instances from the same location in the same RNA of the same organism are available. The script requires the CSV file generated in step 9 and an integer input to specify the register: enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_10_redundancy_check.py 'results/data-from-step-9.csv' '1' or '2'
```
### 11. step_11_analyze_map_model_cc.py
This script will extract the correlation coefficient between the electron density map and modeled structure (map-model cc) for the nucleobases forming shifted wobbles, as well as the mean and median map-model cc for all residues within the chain containing the corresponding shifted wobble. This script will take the CSV file generated in step 10 and the directory containing the calculated map-model cc files (as .txt or .csv format) as input. 
```sh
python step_11_analyze_map_model_cc.py 'results/data-from-step-10.csv' '/directory_of_map_model_cc_files/'
```
### 12. step_12_identifying_structural_clusters.py
This script will identify structural clusters within the non-redundant dataset of the shifted wobbles. For that first, it will generate clipped motifs with the user-defined flanking sequence length which is currently one. These motifs will be stored in a directory which must be specified as an input. The generated structures will be aligned with each other to generate a distance matrix. This distance matrix will then be used to perform a hierarchical clustering. Currently, the distance cut-off is 1.23 Å and requires at least 4 members within a group to be identified as a cluster. 
```sh
python step_12_identifying_structural_clusters.py 'results/data-from-step-11.csv' '/directory_to_store_clipped_structures/'
```
### 13. step_13_stem_check.py
This script will take the output from step 12 and identify the location of the wobbles in the secondary structure motifs. The other inputs will be the directory where all the structures and DSSR output (generated in step 2) are stored, and an integer specifying either standard ('1') or shifted wobble ('2'). 
```sh
python step_13_stem_check.py 'results/data-from-step-12.csv' '/directory_with_all_DSSR_output/' '1' or '2'
```
### 14. step_14_non_WCF_check.py
This script will identify non-Watson-Crick-Franklin interactions within 3.4 Å  of G(N7), G(O6), G(N1), G(N2), G(N3), G(O4'), G(O2'), U(O4), U(N3), U(O2), U(O4'), and U(O2') atoms. The script requires the CSV file generated in step 13 and an integer input to specify the register: enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_14_non_WCF_check.py 'results/data-from-step-13.csv' '1' or '2'
```


