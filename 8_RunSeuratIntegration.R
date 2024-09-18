RunSeuratIntegration <- function(Merged_matrix_file ="../Data/Merged_matrix.txt",
                     sample_batch ="../Data/CellID_batch.txt",
                     base_name="SrtV3_intgration",
                     k_filter = 200){
  
  
  library(Seurat)  
  library(dplyr)
  library(Matrix)
  
  
  Merged_matrix <- read.table(Merged_matrix_file,sep="\t",header=T,row.names=1,check.names = F)
  
  metadata <- read.table(sample_batch,sep="\t",header=T,row.names=1,check.names = F)
  
  batch_element <- levels(metadata$batch)
  
  
  bone_marrow <- CreateSeuratObject(counts=Merged_matrix, meta.data = metadata, min.cells = 3, min.features = 200)
  
  bone_marrow.list <- SplitObject(object = bone_marrow, split.by = "batch")
  
  for (i in 1:length(x = bone_marrow.list)) {
    bone_marrow.list[[i]] <- NormalizeData(object = bone_marrow.list[[i]], verbose = FALSE)
    bone_marrow.list[[i]] <- FindVariableFeatures(object = bone_marrow.list[[i]], 
                                                  selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  
  reference.list <- bone_marrow.list[batch_element]
  bone_marrow.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, k.filter= k_filter) # Might need to reduce k.filter
  bone_marrow.integrated <- IntegrateData(anchorset = bone_marrow.anchors, dims = 1:30)
  
  
  DefaultAssay(object = bone_marrow.integrated) <- "integrated"
  
  save(bone_marrow.integrated, file = paste0(base_name,"Seurat_objet_after_integration.Robj"))
  
  bone_marrow.integrated <- ScaleData(object = bone_marrow.integrated, verbose = FALSE)
  
  bone_marrow.integrated <- RunPCA(object = bone_marrow.integrated, npcs = 50, verbose = FALSE, do.print = TRUE, pcs.print = 1:5, 
                                   genes.print = 5)
  
  write.table(Embeddings(object = bone_marrow.integrated, reduction = "pca"), 
              paste(base_name,"PCA_embedding.txt"),sep="\t",row.names = T,col.names =NA )
  
  numpc <- 50
  ifpng <- FALSE
  perplexity=10
  
  
  chunk <- function(x,n){
    vect <- c(1:x)
    num <- ceiling(x/n)
    split(vect,rep(1:num,each=n,len=x))
  }
  
  
  pdf(paste(base_name,"_PCA.pdf",sep=""),width=9)
  DimPlot(object = bone_marrow.integrated, reduction = "pca")
  dev.off()
  
  if(ifpng){
    png(paste(base_name,"_PCA.png",sep=""),width=2*480/7*9, height = 2*480, res = 2*72)
    PCAPlot(object = bone_marrow.integrated, dim.1 = 1, dim.2 = 2)
    dev.off()
  }
  
  lapply(chunk(numpc,9),function(l){
    pdf(paste(base_name,"_PC",l[1],"-",l[length(l)],"_heatmap.pdf",sep=""),width=9,height=10)
    print(DimHeatmap(object = bone_marrow.integrated, dims = l, cells = 100, balanced = TRUE))
    dev.off()
    
    if(ifpng){
      png(paste(base_name,"_PC",l[1],"-",l[length(l)],"_heatmap.png",sep=""),width=2*480/7*9,height=2*480/7*10, res = 2*72)
      print(DimHeatmap(object = bone_marrow.integrated, dims = l, cells = 100, balanced = TRUE))
      dev.off()
    }
  })
  
  pdf(paste(base_name,"_ElbowPlot.pdf",sep=""))
  print(ElbowPlot(object = bone_marrow.integrated))
  dev.off()
  
  bone_marrow.integrated <- JackStraw(bone_marrow.integrated, reduction ="pca",num.replicate = 100,dims = 50)
  bone_marrow.integrated<- ScoreJackStraw(bone_marrow.integrated,dims = 1:50, do.plot = T)
  
  pdf(paste(base_name,"_JackStraw.pdf"),width=9,height=10)
  print(JackStrawPlot(object = bone_marrow.integrated, dims=1:50))
  dev.off()
  
  pc_pval <- as.data.frame(bone_marrow.integrated@reductions$pca@jackstraw@overall.p.values)
  sig_pc <- which(pc_pval$Score<0.05)
  
  get_1st <- function(x){
    temp <- strsplit(x, "_")[[1]]
    temp_len <- length(temp)
    strsplit(x, "_")[[1]][temp_len]
  }
  
  pca_res <- Embeddings(object = bone_marrow.integrated, reduction = "pca")
  row.names(pca_res) <- sapply(row.names(pca_res),get_1st)
  write.table(pca_res,paste(base_name,"_seurat_pca.txt",sep=""),sep="\t",col.names=NA)
  dir.create(paste(base_name,"_fcs",sep=''))
  write.csv(pca_res,paste(base_name,"_fcs","/seurat_pca_for_DimRed.csv",sep=""),sep=',',row.names = F)
  
  names2 <- as.data.frame(colnames(pca_res))
  colnames(names2) <- "x"
  names2$all <- ""
  names2[sig_pc,"all"] <- "y"
  write.csv(names2,paste(base_name,"_names2.csv",sep=""),sep=",",row.names = F, na = "")
  

  
}


