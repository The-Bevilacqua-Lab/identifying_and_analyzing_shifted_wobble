# Identifying and analyzing shifted wobble
This pipeline will take the search output from the RCSB protein databank to identify and characterize different orientations of wobbles. 
## Setup instructions
- Install [Miniconda](https://docs.anaconda.com/miniconda/install/), and create a new environment
- Install Python >= 3.11
- Install PyMOL by running the following line in terminal
```sh
conda install conda-forge::pymol-open-source=3.0.0
```
- Install required libraries:
```sh
pip install -r requirements.txt
```
- You might need to install EMBOSS Needle manually. First, use the following command to check if it is already installed:
```sh
embossversion
```
Use the following command to install EMBOSS Needle if it is not already installed:
```sh
conda install -c bioconda emboss=6.6.0
```

- Detailed instructions for creating the environments can be found in the env_instructions.rtf file
- For characterizing the 3D structure, Dissecting the Spatial Structure of RNA (DSSR, version: v2.2.1-2021jan12) was used
- The Phenix software package (version: 1.21.1-5286) was used to compare electron density maps and modeled structures.  
## Folder structure

    identify_and_characterize_shifted_wobble
      ├── data                                   #Search output from RCSB protein data bank
      ├── results                                #Output files from each Python script
      ├── scripts                                              
          ├── analysis_scripts                   #All Python scripts to extract and analyze structural features
          ├── plot_scripts                       #All Jupyter notebooks to visualize the outputs 
      ├── requirements.txt                       #Required Python packages and libraries
      ├── env_instructions.rtf                   #Instructions for creating a new environment for this workflow
      ├── README.md                              #This file
      ├── LICENSE
## Instructions for each script
### 1. step_1_data_preparation.py
Before running the steps below, you need to obtain the custom search files in .csv format collected from the RCSB protein data bank (https://www.rcsb.org/). To begin this, start by choosing to make an advance search on the website's homepage. Next, select the 'Polymer Molecular Features' dropdown and select 'Polyer Entity Type'. From there, select the Entity Type to be RNA and then click to run a search. Then, select 'Create Custom Report' from the Tabular Report dropdown located at the top of the screen. Select: PDB ID, Experimental Method, and all of the resolution options. Finally, click the Run Report button to display your custom report. Select to download the files in CSV format. After this, you will see the downloadable link for the splited CSVs for your query. This script will read these CSV files stored in the directory specified by the user and output the cleaner version of the result within user-defined resolution cut-off. Please make sure that there are no extra CSV files with a filename starting with 'rcsb_pdb_custom_report_' except the relevant CSV files. 
```sh
python step_1_data_preparation.py '/directory_for_custom_RCSB_search_results' '3.2'
```
### 2. step_2_download_and_characterize.py
This script downloads the structures and characterizes them using DSSR (version v2.2.1-2021jan12 is required for this pipeline). You need to specify the directory where the CSV file (which will be stored in results by default) generated from step_1_data_preparation.py is stored and the directory where the downloaded structure files and JSON output files from DSSR characterization will be saved.
```sh
python step_2_download_and_characterize.py 'results/data-from-step-1.csv' '/directory_for_structures/'
```
### 3. step_3_filter_out_bp.py
This script filters out all the base pairs of interest within the structures downloaded and characterized in step 2. To run the script you need the same CSV file used in step 2, the directory of the structure and JSON output files, and the base pair of interest (order and case do not matter, 'AU', 'Ua', 'au', and 'uA' all indicates the same base pair). The output will be stored in the 'result' folder. 
```sh
python step_3_filter_out_bp.py 'results/data-from-step-1.csv' '/directory_for_structures/' 'GU'
```
### 4. step_4_filter_out_registers.py
This script will identify the different registers for a base of interest based on the hydrogen bond required to form that pair. To run this script you need the CSV file generated in the previous step (step_3_filter_out_bp.py) and the directory of the structure and JSON files. The script will ask for the hydrogen bond. For standard and shifted wobble, you can write 'G.O6-U.N3,G.N1-U.O2' and 'G.N1-U.O4,G.N2-U.N3', respectively. Do not add any quotation mark to specify the hydrogen bonds. There should not be any space in between the hydrogen bonds, and multiple hydrogen bonds should only be separated by ','.
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
This script will check the hydrogen bond quality and filter out those with poor quality. In this step the source organism and RNA type (column name: 'Molecule') items will also be generalized, therefore, examples such as '16S rRNA', '16S ribosomal RNA', or '16S RIBOSOMAL RNA' will be considered the same RNA in the later steps. However, after executing this step, you might need to manually check for any inconsistency in organism or RNA type information, specifically, the RNA type examples with 'RNA (X-MER)' where X is the length of RNA, if X is comparable with the lengths other RNA types from the same organism, you can replace the 'RNA (X-MER)' with the name of other RNA types. The script requires the CSV file generated in step 6 and an integer input to specify the register: enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_7_filter_and_quality_check.py 'results/data-from-step-6.csv' '1' or '2'
```
### 8. step_8_get_fastas.py
This script takes the data from step 7 and pulls fasta formatted sequences from pymol and the pdb server, as well as makes a reference fasta with one representative sequence per RNA type and organism. The output files are placed in the folder where the python script was run.

```sh
python step_8_get_fastas.py -i 'results/data-from-step-7.csv' -ft pymol -o pymol
python step_8_get_fastas.py -i 'results/data-from-step-7.csv' -ft pdb -o pdb
python step_8_get_fastas.py -i 'results/data-from-step-7.csv' -ft ref -pdb pdb.fasta -o reference
```

### 9. step_9_adjust_res_index.py
This script takes the csv from step 7 along with three fasta files derived from step_8_get_fastas.py. The ouput is the same csv from step 7 with additional columns for: ref_org, pdb_adj_res1, pdb_adj_res2, flank1, flank2, ref_pdb_chain, adj_res1, and adj_res2

```sh
python step_9_adjust_res_index.py -i 'results/data-from-step-7.csv' -py pymol.fasta -pdb pdb.fasta -ref reference.fasta -o 'results/data-from-step-9.csv'
```

### 10. step_10_redundancy_check.py
This script will identify representative structures from redundant wobble examples, where multiple wobble instances from the same location in the same RNA of the same organism are available. The script requires the CSV file generated in step 9 and an integer input to specify the register: enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_10_redundancy_check.py 'results/data-from-step-9.csv' '1' or '2'
```
### 11. step_11_prepare_structure_files.py
This script prepares unique structural files for the non-redundant wobble dataset, where each structure may contain multiple standard or shifted wobbles. As part of the preparation process, all modified residues are removed to ensure compatibility with Phenix. This allows Phenix to utilize the structure and the raw electron density map file to calculate map-model correlation coefficients. The inputs for this script are the CSV file generated in Step 10 and the directory designated for storing the clipped structure files.
```sh
python step_11_prepare_structure_files.py 'results/data-from-step-10.csv' '/directory_for_clipped_structures/'
```
### 12. step_12_analyze_map_model_cc.py
Before executing this script, files containing the correlation coefficient between the electron density maps and modeled structures (map-model cc) must be generated. This involves comparing the clipped structures created in Step 11 with the corresponding electron density files obtained from the [RCSB Protein Data Bank](https://www.rcsb.org/) using the Phenix software package (version: 1.21.1-5286). This script will extract the map-model cc for the nucleobases forming shifted wobbles, as well as the mean and median map-model cc for all residues within the chain containing the corresponding shifted wobble. This script will take the CSV file generated in step 10 and the directory containing the calculated map-model cc files (as .txt or .csv format) as input. 
```sh
python step_12_analyze_map_model_cc.py 'results/data-from-step-10.csv' '/directory_of_map_model_cc_files/'
```
### 13. step_13_identifying_structural_clusters.py
This script will identify structural clusters within the non-redundant dataset of the shifted wobbles. For that first, it will generate clipped motifs with the user-defined flanking sequence length which is currently one. These motifs will be stored in a directory which must be specified as an input. The generated structures will be aligned with each other to generate a distance matrix. This distance matrix will then be used to perform a hierarchical clustering. Currently, the distance cut-off is 1.23 Å and requires at least 4 members within a group to be identified as a cluster. 
```sh
python step_13_identifying_structural_clusters.py 'results/data-from-step-12.csv' '/directory_for_motifs/' 'length of flanking sequences above and below of two residues forming wobble of interest'
```
### 14. step_14_stem_check.py
This script will take the output from step 13 and identify the location of the wobbles in the secondary structure motifs. The other inputs will be the directory where all the structures and DSSR output (generated in step 2) are stored, and an integer specifying either standard ('1') or shifted wobble ('2'). 
```sh
python step_14_stem_check.py 'results/data-from-step-13.csv' '/directory_with_all_DSSR_output/' '1' or '2'
```
### 15. step_15_non_WCF_check.py
This script will identify non-Watson-Crick-Franklin interactions within 3.4 Å  of G(N7), G(O6), G(N1), G(N2), G(N3), G(O4'), G(O2'), U(O4), U(N3), U(O2), U(O4'), and U(O2') atoms. The script requires the CSV file generated in step 14 and an integer input to specify the register: enter '1' for the standard G•U wobble or '2' for the shifted G•U wobble.
```sh
python step_15_non_WCF_check.py 'results/data-from-step-14.csv' '1' or '2'
```


