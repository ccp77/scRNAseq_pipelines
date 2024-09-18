seurat_analysis <- function(TPM_file=TPM_sub,
                            # TPM_file="TPM.txt",
                            sample_file= NULL,
                            group_col="GroupID",
                            base_name="GC", 
                            resolution = 0.85,
                            mt="^Mt-", 
                            scale=F,
                            scale.factor = 1e2,
                            perplexity=10,
                            deg_method="bimod",
                            col.low = "#FF00FF",col.mid = "#000000", col.high = "#FFFF00",
                            ifpng = F,
                            logfc.threshold = 0,
                            x.low.cutoff = 0.0125, y.cutoff = 0.5,
                            ifdeg=F,
                            numpc = 20, 
                            JackStraw= T,
                            Louvain =F,
                            RunTSNE=F, 
                            RegressOutCellCycle_Mouse =F,
                            RegressOutCellCycle_Human =F){            
  # library(xlsx)
  library(Seurat)
  library(dplyr)
  library(Matrix)
  # pbmc.data <- read.table(TPM_file,sep="\t",header=T,row.names=1,check.names = F)
  pbmc.data <- TPM_file
  
  if(!is.null(sample_file)){
    sample <- read.table(sample_file,sep="\t",header=T,row.names=1)
    colnames(pbmc.data) <- paste(sample[colnames(pbmc.data),group_col],colnames(pbmc.data),sep="_")
  }
  
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = base_name, min.cells = 3, min.features = 200)
  
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = mt)
  
  
  pdf(paste(base_name,"_QC.pdf",sep=""))
  print(VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
  
  if(ifpng){
    
    png(paste(base_name,"_QC.png",sep=""),height = 2*480, width=2*480,res = 2*72)
    print(VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
    dev.off()
  }
  

    pdf(paste(base_name,"_QC_nUMI_mito_nGene.pdf",sep=""))
    plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    CombinePlots(plots = list(plot1, plot2))
    dev.off()

  # pbmc <- subset(pbmc, subset = nFeature_RNA >200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize")
  pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(pbmc),10)
  
  pdf(paste(base_name,"_variableGenes.pdf",sep=""))
  plot1 <- VariableFeaturePlot(pbmc)
  plot2 <- LabelPoints(plot = plot1, points = top10, repe1=T)
  print(plot2)
  dev.off()
  
  if(ifpng){
    png(paste(base_name,"_variableGenes.png",sep=""),height = 480, width = 480, res = 72)
    plot1 <- VariableFeaturePlot(pbmc)
    plot2 <- LabelPoints(plot = plot1, points = top10, repe1=T)
    print(plot2)
    dev.off()
  }
  
  write.table(VariableFeatures(pbmc),paste(base_name,"_most_variable_genes.txt",sep=""), sep="\t", row.names  = F)
  
  if(RegressOutCellCycle_Mouse==T) {
    print("Mouse")
    
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }

    s.genes <- cc.genes$s.genes
    s.genes_lower <- tolower(s.genes)
    s.genes_final <- firstup(s.genes_lower)

    g2m.genes <- cc.genes$g2m.genes
    g2m.genes_lower <- tolower(g2m.genes)
    g2m.genes_final <- firstup(g2m.genes_lower)

    pbmc <- CellCycleScoring(pbmc, s.features = s.genes_final, g2m.features = g2m.genes_final, set.ident = TRUE)
    pbmc <- ScaleData(pbmc, verbose = T,vars.to.regress = c("S.Score", "G2M.Score"))

  } else if(RegressOutCellCycle_Human==T){
    print("Human")
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    file.integrated <- CellCycleScoring(file.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    file.integrated <- ScaleData(file.integrated, verbose = T,vars.to.regress = c("S.Score", "G2M.Score"))

  } else {
    print("none")
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)
    
  }
  
  
  
  pbmc <- RunPCA(object = pbmc, 
                 features = VariableFeatures(object = pbmc), 
                 do.print = TRUE, 
                 npcs = numpc)
  
  print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
  
  chunk <- function(x,n){
    vect <- c(1:x)
    num <- ceiling(x/n)
    split(vect,rep(1:num,each=n,len=x))
  }
  
  pdf(paste(base_name,"_PCA.pdf",sep=""),width=9)
  DimPlot(object = pbmc, reduction= "pca")
  dev.off()
  
  if(ifpng){
    png(paste(base_name,"_PCA.png",sep=""),width=2*480/7*9, height = 2*480, res = 2*72)
    DimPlot(object = pbmc, reduction= "pca")
    dev.off()
  }
  
  lapply(chunk(numpc,9),function(l){
    pdf(paste(base_name,"_PC",l[1],"-",l[length(l)],"_heatmap.pdf",sep=""),width=9,height=10)
    DimHeatmap(object = pbmc, dims = l, cells = 500, balanced = TRUE)
    dev.off()
    
    if(ifpng){
      png(paste(base_name,"_PC",l[1],"-",l[length(l)],"_heatmap.png",sep=""),width=2*480/7*9,height=2*480/7*10, res = 2*72)
      DimHeatmap(object = pbmc, dims = l, cells = 500, balanced = TRUE)
      dev.off()
    }
  })
  
  if(JackStraw == T) {
    print("Do JackStraw")
    pbmc <- JackStraw(object = pbmc, num.replicate = 100, dims = numpc)
    pbmc <- ScoreJackStraw(pbmc, dims = 1:numpc)

    pdf(paste(base_name,"_JackStraw.pdf",sep=""),width=9,height=10)
    print(JackStrawPlot(object = pbmc,dims = 1:numpc))
    dev.off()

    if(ifpng){
      png(paste(base_name,"_JackStraw.png",sep=""),width=2*480/7*9,height=2*480/7*10, res = 2*72)
      print(JackStrawPlot(object = pbmc,dims = 1:numpc))
      dev.off()
    }

    pc_pval <- as.data.frame(pbmc@reductions$pca@jackstraw@overall.p.values)
    sig_pc <- which(pc_pval$Score<0.05)
  } else {
    print("Don't do JackStraw")
    sig_pc <- 1:50
  }

  save(pbmc, file = paste(base_name,"_seurat_beforeClustering.Robj",sep=""))
  
  pdf(paste(base_name,"_PC_SD.pdf",sep=""))
  print(ElbowPlot(object = pbmc))
  dev.off()
  
  if(ifpng){
    png(paste(base_name,"_PC_SD.png",sep=""), width = 2*480, height = 2*480, res = 2*72)
    print(ElbowPlot(object = pbmc))
    dev.off()
  }
  
  get_1st <- function(x){
    temp <- strsplit(x, "_")[[1]]
    temp_len <- length(temp)
    strsplit(x, "_")[[1]][temp_len]
  }
  
  
  if(RunTSNE==T) {
    pbmc <- RunTSNE(object = pbmc, dims.use = sig_pc, do.fast = TRUE,perplexity=perplexity)
    pdf(paste(base_name,"_tSNE.pdf",sep=""),width=9)
    TSNEPlot(object = pbmc)
    dev.off()
    
    if(ifpng) {
      png(paste(base_name,"_tSNE.png",sep=""),width=2*480/7*9,height = 2*480, res = 2*72)
      TSNEPlot(object = pbmc)
      dev.off()
    }
    tsne_res <- Embeddings(object = pbmc, reduction = "tsne")
    row.names(tsne_res) <- sapply(row.names(tsne_res),get_1st)
    write.table(tsne_res,paste(base_name,"_seurat_tsne.txt",sep=""),sep="\t",col.names=NA)
    
  }
  
  if(Louvain == T){
    pbmc <- FindNeighbors(pbmc, dims = 1:last(sig_pc))
    pbmc <- FindClusters(object = pbmc, 
                         resolution = resolution)
    
    
    pbmc <- RunTSNE(object = pbmc, dims.use = sig_pc, do.fast = TRUE,perplexity=perplexity)
    
    pdf(paste(base_name,"_tSNE_cluster.pdf",sep=""),width=9)
    TSNEPlot(object = pbmc)
    dev.off()
    
    if(ifpng){
      png(paste(base_name,"_tSNE_cluster.png",sep=""),width=2*480/7*9, height = 2*480, res = 2*72)
      TSNEPlot(object = pbmc)
      dev.off()
    }
    
    pdf(paste(base_name,"_PCA_cluster.pdf",sep=""),width=9)
    PCAPlot(object = pbmc)
    dev.off()
    
    if(ifpng){
      png(paste(base_name,"_PCA_cluster.png",sep=""),width=2*480/7*9,height = 2*480, res = 2*72)
      PCAPlot(object = pbmc)
      dev.off()
    }
    cluster_res <- Idents(object = pbmc)
    names(cluster_res) <- sapply(names(cluster_res),get_1st)
    write.table(data.frame(cell=names(cluster_res),cluster=cluster_res),paste(base_name,"_seurat_cluster.txt",sep=""),sep="\t",row.names=F,quote=F)
    
  }
  

  pca_res <- Embeddings(object = pbmc, reduction = "pca")
  #row.names(pca_res) <- sapply(row.names(pca_res),get_1st)
  write.table(pca_res,paste(base_name,"_seurat_pca.txt",sep=""),sep="\t",col.names=NA)
  dir.create(paste(base_name,"_fcs",sep=''))
  write.csv(pca_res,paste(base_name,"_fcs","/seurat_pca_for_DimRed.csv",sep=""),sep=',',row.names = F)
  
  names2 <- as.data.frame(colnames(pca_res))
  colnames(names2) <- "x"
  names2$all <- ""
  names2[sig_pc,"all"] <- "y"
  write.csv(names2,paste(base_name,"_names2.csv",sep=""),sep=",",row.names = F, na = "")
  
  
  
  save(pbmc, file = paste(base_name,"_seurat.Robj",sep=""))
  
  if(ifdeg){
    pbmc.markers <- FindAllMarkers(object = pbmc, test.use = deg_method, logfc.threshold = logfc.threshold)
    pbmc.markers$p_val_adj_BH <- p.adjust(pbmc.markers$p_val, method = "BH")
    
    deg <- pbmc.markers[pbmc.markers$p_val_adj<0.05,]
    deg <- deg[order(deg$cluster,deg$avg_logFC,decreasing = T),]
    
    write.table(deg,paste(base_name,"_seurat_bimod_DEG.txt",sep=""),sep="\t",col.names=NA)
    #  write.xlsx(deg,paste(base_name,"_seurat_bimod_DEG.xls",sep=""))
    
    write.table(pbmc.markers,paste(base_name,"_seurat_bimod.txt",sep=""),sep="\t",col.names=NA)
    # write.xlsx(pbmc.markers,paste(base_name,"_seurat_bimod.xls",sep=""))
    
    
    # pbmc.markers %>% group_by(cluster) %>% top_n(5, avg_logFC) -> temp
    # 
    # for(s in unique(temp$cluster)){
    #   
    #   pdf(paste(base_name,"_cluster",s,"_DEG_violin.pdf",sep=""))
    #   VlnPlot(pbmc, temp$gene[temp$cluster==s])
    #   dev.off()
    #   
    #   png(paste(base_name,"_cluster",s,"_DEG_violin.png",sep=""),width = 2*480, height = 2*480, res = 2*72)
    #   VlnPlot(pbmc, temp$gene[temp$cluster==s])
    #   dev.off()
    #   
    # }
    
    # sapply(unique(temp$cluster),function(s){
    #   
    #   pdf(paste(base_name,"_cluster",s,"_DEG_tSNE.pdf",sep=""))
    #   FeaturePlot(pbmc, temp$gene[temp$cluster==s],cols.use = c("grey","blue"))
    #   dev.off()
    #   
    #   if(ifpng){
    #     png(paste(base_name,"_cluster",s,"_DEG_tSNE.png",sep=""),width = 2*480, height = 2*480, res = 2*72)
    #     FeaturePlot(pbmc, temp$gene[temp$cluster==s],cols.use = c("grey","blue"))
    #     dev.off()
    #   }
    #   
    # })
    
    #pbmc.markers %>% group_by(cluster) %>% top_n(100, avg_logFC) -> top10
    
    # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
    deg <- pbmc.markers[pbmc.markers$p_val_adj<0.05 & pbmc.markers$avg_logFC>0,]
    deg <- deg[order(deg$cluster,deg$avg_logFC,decreasing = T),]
    pdf(paste(base_name,"_cluster_DEG_hm.pdf",sep=""))
    print(DoHeatmap(object = pbmc, genes.use = row.names(deg), slim.col.label = TRUE, remove.key = TRUE),group.label.rot =T)
    dev.off()
    
    if(ifpng){
      png(paste(base_name,"_cluster_DEG_hm.png",sep=""),width = 2*480, height = 2*480, res = 2*72)
      print(DoHeatmap(object = pbmc, genes.use = row.names(deg), slim.col.label = TRUE, remove.key = TRUE),group.label.rot =T)
      dev.off()
    }
  }
  ####################
}


JackStraw_pval <- function (object, PCs = 1:5, score.thresh = 1e-05) {
  library(reshape2)
  
  pAll <- GetDimReduction(object, reduction.type = "pca", slot = "jackstraw")@emperical.p.value
  pAll <- pAll[, PCs, drop = FALSE]
  pAll <- as.data.frame(pAll)
  pAll$Contig <- rownames(x = pAll)
  pAll.l <- melt(data = pAll, id.vars = "Contig")
  colnames(x = pAll.l) <- c("Contig", "PC", "Value")
  qq.df <- NULL
  score.df <- NULL
  for (i in PCs) {
    q <- qqplot(x = pAll[, i], y = runif(n = 1000), plot.it = FALSE)
    pc.score <- suppressWarnings(prop.test(x = c(length(x = which(x = pAll[, 
                                                                           i] <= score.thresh)), floor(x = nrow(x = pAll) * 
                                                                                                         score.thresh)), n = c(nrow(pAll), nrow(pAll)))$p.val)
    if (length(x = which(x = pAll[, i] <= score.thresh)) == 
        0) {
      pc.score <- 1
    }
    if (is.null(x = score.df)) {
      score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
    }
    else {
      score.df <- rbind(score.df, data.frame(PC = paste0("PC", 
                                                         i), Score = pc.score))
    }
    if (is.null(x = qq.df)) {
      qq.df <- data.frame(x = q$x, y = q$y, PC = paste0("PC", 
                                                        i))
    }
    else {
      qq.df <- rbind(qq.df, data.frame(x = q$x, y = q$y, 
                                       PC = paste0("PC", i)))
    }
  }
  return(score.df)
}

#seurat_analysis(TPM_file="TPM.txt",sample_file="sample.txt",group_col="PatientID",base_name="GC_PatientID", resolution = 0.85)
#seurat_analysis(TPM_file="TPM.txt",sample_file="sample.txt",group_col="Phenotype",base_name="GC_Phenotype", resolution = 0.85)
#seurat_analysis(TPM_file="TPM.txt",sample_file="sample.txt",group_col="Site",base_name="GC_Site", resolution = 0.85)
#seurat_analysis(TPM_file="TPM.txt",sample_file="sample.txt",group_col="Batch",base_name="GC_Batch", resolution = 0.85)


#seurat_analysis(TPM_file="TPM_geneSymbol_mean_selectedMarker.txt",sample_file=NULL,group_col=NULL,base_name="Bat_Bcell_selectedMarker", resolution = 0.85)


