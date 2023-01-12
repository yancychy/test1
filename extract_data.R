## extract DNA methylation data (beta values) and SNP data (allele) for overlapping SNPs and CpG probe

extract_methyl_data <-function(INPUT_D, OUTPUT_R, OUTPUT_T, window_nt){
    DNAmethylData <- local(get(load(paste(OUTPUT_R,"combined_all_DNAmethyl_data.RData",sep=""))))
    DATA_Table <- paste(INPUT_D,"DNA_meth_450/FILE_SAMPLE_MAP.txt",sep="")
    metadata_T <- read.table(DATA_Table,sep= "\t",header=T,stringsAsFactors = FALSE)
    colnames(DNAmethylData) <- substr(metadata_T$barcode.s.,1,12)
    out_CpG_SNP <- paste(OUTPUT_T, "out_CpG_SNP_", window_nt,"nt_filter_var.bed",sep="")
    out_CpG_SNP_DATA <- read.table(out_CpG_SNP, sep= "\t", header=T)
    matrix_CpG_overlap <- DNAmethylData[rownames(DNAmethylData) %in% out_CpG_SNP_DATA$CpG_ID,] 
    save(matrix_CpG_overlap, file = paste(OUTPUT_R,"matrix_CpG_",window_nt,"nt_overlap.RData", sep=""))
    write.table(matrix_CpG_overlap, file = paste(OUTPUT_T,"matrix_CpG_",window_nt,"nt_overlap.txt", sep=""), sep = "\t", quote=F)
}

extract_methyl_data_sliding_window <-function(INPUT_D, OUTPUT_R, OUTPUT_T, max_boundary, window_nt, i){
    DNAmethylData <- local(get(load(paste(OUTPUT_R,"combined_all_DNAmethyl_data.RData",sep=""))))
    DATA_Table <- paste(INPUT_D,"DNA_meth_450/FILE_SAMPLE_MAP.txt",sep="")
    metadata_T <- read.table(DATA_Table,sep= "\t",header=T,stringsAsFactors = FALSE)
    colnames(DNAmethylData) <- substr(metadata_T$barcode.s.,1,12)
    out_CpG_SNP <- paste(OUTPUT_T,"out_CpG_SNP_",max_boundary*2,"nt_",window_nt,"nt_",i,".bed" ,sep="")
    out_CpG_SNP_DATA <- read.table(out_CpG_SNP, sep= "\t", header=T)
    matrix_CpG_overlap <- DNAmethylData[rownames(DNAmethylData) %in% out_CpG_SNP_DATA$CpG_ID,] 
    save(matrix_CpG_overlap, file = paste(OUTPUT_R,"matrix_CpG_",max_boundary*2,"nt_",window_nt,"nt_",i,"_overlap.RData" ,sep=""))
    write.table(matrix_CpG_overlap, file = paste(OUTPUT_T,"matrix_CpG_",max_boundary*2,"nt_",window_nt,"nt_",i,"_overlap.txt" ,sep=""), sep = "\t", quote=F)
}

extract_SNP_data <- function(INPUT_D, OUTPUT_R, OUTPUT_T, window_nt){
    SNP_Data <- local(get(load(paste(OUTPUT_R,"combined_all_SNP_data.RData",sep=""))))
    SNP_Data<-SNP_Data[!(duplicated(SNP_Data$dbsnp_rs_id)),]
    rownames(SNP_Data) <- SNP_Data$dbsnp_rs_id
    SNP_Data$Row.names <- NULL
    SNP_Data<-SNP_Data[,grep(".CEL",colnames(SNP_Data))]
    DATA_Table <- paste(INPUT_D, "SNP_array/FILE_SAMPLE_MAP.txt", sep="")
    metadata_T <- read.table(DATA_Table,sep= "\t",header=T,stringsAsFactors = FALSE)
    colnames(SNP_Data) <- substr(metadata_T$barcode.s.,1,12)
    
    out_CpG_SNP <- paste(OUTPUT_T, "out_CpG_SNP_", window_nt,"nt_filter_var.bed",sep="")
    out_CpG_SNP_DATA <- read.table(out_CpG_SNP, sep= "\t", header=T)
    matrix_SNP_overlap <- SNP_Data[rownames(SNP_Data) %in% out_CpG_SNP_DATA$SNP_ID,] 
    save(matrix_SNP_overlap, file = paste(OUTPUT_R,"matrix_SNP_",window_nt,"nt_overlap.RData", sep=""))
    write.table(matrix_SNP_overlap, file = paste(OUTPUT_T,"matrix_SNP_",window_nt,"nt_overlap.txt", sep=""), sep = "\t", quote=F)
}

extract_SNP_data_sliding_window <- function(INPUT_D, OUTPUT_R, OUTPUT_T, max_boundary, window_nt, i){
    SNP_Data <- local(get(load(paste(OUTPUT_R,"combined_all_SNP_data.RData",sep=""))))
    SNP_Data<-SNP_Data[!(duplicated(SNP_Data$dbsnp_rs_id)),]
    rownames(SNP_Data) <- SNP_Data$dbsnp_rs_id
    SNP_Data$Row.names <- NULL
    SNP_Data<-SNP_Data[,grep(".CEL",colnames(SNP_Data))]
    DATA_Table <- paste(INPUT_D, "CNV_SNP_T/FILE_SAMPLE_MAP.txt", sep="")
    metadata_T <- read.table(DATA_Table,sep= "\t",header=T,stringsAsFactors = FALSE)
    colnames(SNP_Data) <- substr(metadata_T$barcode.s.,1,12)
    out_CpG_SNP <- paste(OUTPUT_T,"out_CpG_SNP_",max_boundary*2,"nt_",window_nt,"nt_",i,".bed" ,sep="")
    out_CpG_SNP_DATA <- read.table(out_CpG_SNP, sep= "\t", header=T)
    matrix_SNP_overlap <- SNP_Data[rownames(SNP_Data) %in% out_CpG_SNP_DATA$SNP_ID,] 
    save(matrix_SNP_overlap, file = paste(OUTPUT_R,"matrix_SNP_",max_boundary*2,"nt_",window_nt,"nt_",i,"_overlap.RData" ,sep=""))
    write.table(matrix_SNP_overlap, file = paste(OUTPUT_T,"matrix_SNP_",max_boundary*2,"nt_",window_nt,"nt_",i,"_overlap.txt", sep=""), sep = "\t", quote=F)
}