############################################################################################################
########################## Master Script scRNAseq - Seurat V3  ##########################
############################################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~scRNAseq analysis one dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("~/path/1_SrtV3_singledataset_analysis.R")
dir.create("Srt")
Matrix <- read.table("TPM_file.txt",sep="\t",header=T,row.names=1,check.names = F)
seurat_analysis(TPM_file=Matrix,                                  # TPM_file.txt
                base_name="Srt",                                       # Name given to all output files
                numpc=50,                                                     # numpc = number of PCs calculated by Jackstraw 
                resolution = 0.85)                                            # KNN clustering resolutio   

source("~/path/2_1_runENfcsFuncs.R")
dir.create("DimRed")
DimRed(outputSuffix = "DimRed/Srt_minD0.05_k20",
       DotSNE = F,
       tSNEperpelxity = 30 ,     #default = 30; increased perplexity => increased spread
       DoUMAP = T,
       min_dist_UMAP = 0.05,          #default = 0.2
       n_neighbors_UMAP = 15,         #default = 15
       Dophenograph = T,             #clusters cells using Rphenograpy. 
       kValue = 20, 
       names2_file = "Srt/Srt_names2.csv",
       embedding_file_folder = "Srt/Srt_fcs/")



dir.create("GLOBAL")
source("~/path/4_Merge_Metadata.R")
Merge_metadata(Input_folders = c("Metadata/"),# list all the folders that contain Metadata for the cells you're intersted in
               output_file_name= "GLOBAL/name.txt") 






