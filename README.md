# scRNAseq_pipelines
These scripts were written to run key scRNAseq analysis steps as "blocks" that can be called from a master script. 
This facilitates the analysis of scRNAseq datasets for researcher with little coding experience. 
In the script 0_Master_script, reserachers can call different functions and specify inputs such as path to data or parameters. 
The functions' output differ depending on the function and can include csv files or other outputs. 
These functions were written to accumplish specific functions, as identified by the project needs and resources at that time.

Scripts 2_1_runENfcsFuncs.R and 2_2_ENfcsFuncs were adapted from https://www.nature.com/articles/nbt.4314 and written by Evan W Newell as indicated in the respective scripts.
