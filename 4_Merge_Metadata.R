Merge_metadata <-function(Input_folders = c("../path1/", "../path2/"),
                          output_file_name= "../path/name.txt") {
  
  list_folders <- as.list(Input_folders)
  
  Merged_tables <- lapply(list_folders,function(l){
    files <- list.files(path = l)
    files_txt <- grep(pattern = ".txt", x= files, value = T)
    
    tables <- lapply(files_txt,function(l2){
      path <- paste0(l,l2)
      table <- read.table(path, header = T,sep = "\t", check.names = F)
      colnames(table)[1] <- "CellID"
      return(table)
    })
    
    merged_table <- Reduce(function(x, y){merge(x, y, by ="CellID")},tables)
    return(merged_table)
  })
  
  Merged_table <- Reduce(function(x, y){merge(x, y, by ="CellID")},Merged_tables)
  colnames(Merged_table)
  Merged_table[,(ncol(Merged_table)+1)] <- 1:nrow(Merged_table)
  splitted <-unlist(strsplit(output_file_name,"/"))
  splitted <- splitted[length(splitted)]
  splitted <- gsub(pattern = ".txt","",x=splitted)
  colnames(Merged_table)[ncol(Merged_table)] <- paste0("equ_nb_",splitted)
  Merged_table_noPCs <- Merged_table[,-grep("PC", colnames(Merged_table))]
  write.table(Merged_table,output_file_name,sep="\t",row.names = F)
  write.table(Merged_table_noPCs,paste0(output_file_name,"noPCs.txt"),sep="\t",row.names = F)
  write.csv(Merged_table[,-1],paste0(output_file_name,"_NUM.csv"),row.names = F)
  write.csv(Merged_table_noPCs[,-1],paste0(output_file_name,"_NUM_noPCs.csv"),sep="\t",row.names = F)
  
  
}

