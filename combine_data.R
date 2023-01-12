## combine all SNP data and DNA methylation data

combine_all_DNA_methyl_data <- function(INPUT_D, OUTPUT_R, OUTPUT_T) {
    ## combine all beta values and then calculate mean and variance
    
    DNAmethyl_DATA <- paste(INPUT_D,"DNA_meth_450/Level_3/",sep="")
    listxt<-dir(path=DNAmethyl_DATA, pattern="*.txt")
    
    # process one file and prepare a matrix
    temp <- read.delim(paste(DNAmethyl_DATA,listxt[1],sep=""), sep="", header=T)
    File_size<-length(listxt); S_ncol<-ncol(temp); S_nrow<-nrow(temp)
    DNAmethyl_loc<-1; rawData<-temp[2:S_nrow,DNAmethyl_loc]
    
    x <- matrix(data=NA, nrow=S_nrow-1, ncol=File_size)  
    rownames(x)<-rownames(temp)[2:S_nrow];  x[,1]<-as.character(rawData)
    
    for (i in 2:length(listxt)){
        temp<-read.delim(paste(DNAmethyl_DATA,listxt[i],sep=""), sep="", header=T)
        rawData<-temp[2:S_nrow,1]
        x[,i]<-as.character(rawData); print(i)
    }
    
    y<-apply(x,2,as.numeric); rownames(y)<-rownames(x); y <- data.frame(y)
    TCGA_head<-46; TCGA_tail<-57
    for (i in 1:length(listxt)) { colnames(y)[i]<-substring(listxt[i],TCGA_head,TCGA_tail) }
    y<-y[-grep("rs",rownames(y)),]; y<-na.omit(y)  # remove SNPs and NAs
    
    # compute mean and variance
    options(scipen=500)
    y_mean<-apply(y,1,mean); y_var<-apply(y,1,var);
    y_mean_var<-NULL; y_mean_var$CpG_ID<-rownames(y);  y_mean_var$mean<-y_mean; y_mean_var$var<-y_var
    y_mean_var<-data.frame(y_mean_var)
    y_mean_var<-y_mean_var[-grep("rs", y_mean_var[,1]),];  y_mean_var<-na.omit(y_mean_var)   # remove SNPs and NAs
    
    save(y, file= paste(OUTPUT_R, "combined_all_DNAmethyl_data.RData", sep=""))
    write.table(y, file= paste(OUTPUT_T,"combined_all_DNAmethyl_data.txt", sep=""),sep = "\t")
    write.table(y_mean_var, file= paste(OUTPUT_T,"combined_DNAmethyl_mean_var.txt", sep=""), sep= "\t", row.names=FALSE)
}

combine_all_SNP_data <- function(HOME, INPUT_D, OUTPUT_R, OUTPUT_T){
    # combine all SNP data
    library(ff); library(crlmm); library(genomewidesnp6Crlmm)
    path2Cels <- paste(INPUT_D,"SNP_array/Level_1",sep="")
    celFiles <- list.celfiles(path2Cels, full.names=TRUE, pattern=".CEL")
    numberOfFiles<-length(celFiles)
    tmp_dir<- paste(HOME,"tmp_dir/",sep="");dir.create(tmp_dir); ldPath(tmp_dir)
    cnSet <- genotype(filenames=celFiles,cdfName="genomewidesnp6",batch=rep("B",numberOfFiles),SNRMin=1,gender=rep(1,numberOfFiles))
    #     print(cnSet$SNR)
    
    pmatrix<-confs(cnSet[1:906600,]) ## confidence score
    temp <- calls(cnSet[1:906600])   ## genotype call
    
    library(pd.genomewidesnp.6)
    con <- pd.genomewidesnp.6@getdb(); dbListTables(con)
    fsnps<-dbGetQuery(con,"select * from featureSet")  ## snp refernce
    rownames(fsnps) <-fsnps[,2]
    SNP_all<-merge(temp,fsnps,by="row.names")
    
    save(SNP_all, file=paste(OUTPUT_R,"combined_all_SNP_data.RData",sep=""))
    write.table(SNP_all, file=paste(OUTPUT_T,"combined_all_SNP_data.txt", sep=""), sep= "\t")
}
