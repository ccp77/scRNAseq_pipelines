# scRNAseq_pipelines
These scripts were written in 2018 to run key scRNAseq analysis steps based mostly on Seurat Version 3 https://satijalab.org/seurat/archive/v3.0/. These scripts can be run as blocks from a master script, by chanigng the inputs and output names, which could greatly facilitate the analysis of scRNAseq datasets for researchers with little coding experience. 

In the script 0_Master_script, reserachers can call different functions that are written in the other scripts included, and specify inputs such as paths to data or parameters. 
The functions' outputs differ depending on the function and can include csv files or other outputs. 
These functions were written to accumplish specific functions, as identified by the project needs and resources at that time.

Most scrips were adapted from Seurat libraries https://satijalab.org/seurat/. 
Scripts 2_1_runENfcsFuncs.R and 2_2_ENfcsFuncs were adapted from https://www.nature.com/articles/nbt.4314 and written by Evan W Newell as indicated in the respective scripts.
