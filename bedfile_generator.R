SNP_bedFile_Gen<-function(OUTPUT_T, window_nt){ 
    options(scipen=500)
    ## input.F: combined SNP data, OUTPUT_T: SNP bed file

    input.F<-read.table(paste(OUTPUT_T, "combined_all_SNP_data.txt", sep=""),sep="\t", header=T, stringsAsFactors=FALSE)
    fourth_clm<- data.frame(paste(input.F$dbsnp_rs_id,input.F$allele_a, input.F$allele_b, sep=" "))
    bedFormat<-cbind(input.F$chrom, as.numeric(input.F$physical_pos), as.numeric(input.F$physical_pos)+1, fourth_clm)
    rm(input.F,fourth_clm)
    bedFormat$score<-0; bedFormat$strand<-"+"
    bedFormat<-na.omit(bedFormat)
    #     bedFormat<-format(bedFormat,scientific=FALSE)
    write.table(bedFormat,file=paste(OUTPUT_T,"SNP_",window_nt,"nt.bed",sep=""),col.names=F, row.names=F,sep="\t", quote=F)
}

DNAmethyl_bedFile_Gen <- function(INPUT_D, OUTPUT_T, window_nt) { 
    options(scipen=500)
    ## prepare DNAmethyl_data.bed file for overlap selection
    ## window_nt is the covered up- and down-stream length (e.g., 25, 50, 100)
    
    window_nt<-as.numeric(window_nt); 
    DNAmethyl_DATA <- paste(INPUT_D,"DNA_meth_450/Level_3/",sep="")
    listxt<-dir(path=DNAmethyl_DATA, pattern="*.txt")    
    input.F <- read.delim(paste(DNAmethyl_DATA,listxt[1],sep=""), sep="\t", header=T, stringsAsFactors=FALSE)
    for (i in 1:ncol(input.F)){ colnames(input.F)[i]<-as.character(input.F[1,i])}
    
    input.mean_var<-read.table(paste(OUTPUT_T,"combined_DNAmethyl_mean_var.txt", sep=""), header=T)
    rownames(input.mean_var)<-input.mean_var$CpG_ID
    input.F<-input.F[-1,]; input.F<-na.omit(input.F)   # remove header row
    rownames(input.F)<-input.F[,1]
    comb.data<-merge(input.F,input.mean_var,by="row.names")
    rownames(comb.data)<-comb.data[,1]
    comb.data<-comb.data[,c(-1,-2)]; comb.data[comb.data==""]<-0
    rm(input.F, input.mean_var)
    header<-c("Chr", "Genomic_start", "Genomic_end", "CpGID","score","strand")
    fourth_clm<- data.frame(paste(rownames(comb.data),comb.data$Gene_Symbol, comb.data$mean, comb.data$var, sep=" "))
    
    bedFormat<-cbind(as.numeric(comb.data$Genomic_Coordinate), as.numeric(comb.data$Genomic_Coordinate)+window_nt)
    bedFormat<-zapsmall(bedFormat)
    bedFormat<-cbind(comb.data$Chromosome, bedFormat, fourth_clm)
    bedFormat$score<-0; bedFormat$strand<-"+"
    colnames(bedFormat)<-header    
    write.table(bedFormat,file=paste(OUTPUT_T,"DNAmeth","_",window_nt,"nt.bed",sep=""),col.names=F, row.names=F, sep="\t", quote=F)
}

DNAmethyl_bedFile_sliding_window_Gen <- function(INPUT_D, OUTPUT_T, window_nt, overlap_nt, max_boundary) { 
    options(scipen=500)
    ## prepare DNAmethyl_data.bed file for overlap selection
    ## window_nt is the covered up- and down-stream length (e.g., 25, 50, 100)
    
    window_nt<-as.numeric(window_nt); overlap_nt<-as.numeric(overlap_nt); max_boundary<-as.numeric(max_boundary)
    DNAmethyl_DATA <- paste(INPUT_D,"DNA_meth_450/Level_3/",sep="")
    listxt<-dir(path=DNAmethyl_DATA, pattern="*.txt")    
    input.F <- read.delim(paste(DNAmethyl_DATA,listxt[1],sep=""), sep="\t", header=T, stringsAsFactors=FALSE)
    for (i in 1:ncol(input.F)){ colnames(input.F)[i]<-as.character(input.F[1,i])}
    
    input.mean_var<-read.table(paste(OUTPUT_T,"combined_DNAmethyl_mean_var.txt", sep=""), header=T)
    rownames(input.mean_var)<-input.mean_var$CpG_ID
    input.F<-input.F[-1,]; input.F<-na.omit(input.F)   # remove header row
    rownames(input.F)<-input.F[,1]
    comb.data<-merge(input.F,input.mean_var,by="row.names")
    rownames(comb.data)<-comb.data[,1]
    comb.data<-comb.data[,c(-1,-2)]; comb.data[comb.data==""]<-0
    rm(input.F, input.mean_var)
    header<-c("Chr", "Genomic_start", "Genomic_end", "CpGID","score","strand")
    fourth_clm<- data.frame(paste(rownames(comb.data),comb.data$Gene_Symbol, comb.data$mean, comb.data$var, sep=" "))
    comb.data<-comb.data[as.numeric(comb.data$Genomic_Coordinate) >= max_boundary,] 
    
    for (i in 1:length(seq(0,max_boundary*2,overlap_nt))){      
        bedFormat<-cbind(as.numeric(comb.data$Genomic_Coordinate)-max_boundary+overlap_nt*(i-1), as.numeric(comb.data$Genomic_Coordinate)-max_boundary+window_nt+overlap_nt*(i-1))
        bedFormat<-zapsmall(bedFormat)
        bedFormat<-cbind(comb.data$Chromosome, bedFormat, fourth_clm)
        bedFormat$score<-0; bedFormat$strand<-"+"
        colnames(bedFormat)<-header    
        write.table(bedFormat,file=paste(OUTPUT_T,"DNAmeth_",max_boundary*2, "nt_",window_nt,"nt_",i,".bed",sep=""),col.names=F, row.names=F, sep="\t", quote=F)
    }
}