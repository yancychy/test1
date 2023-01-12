filter_var_after_overlapSelect <- function(INPUT_F, OUTPUT_T, window_nt){
    CpG_SNP<-read.table(paste(OUTPUT_T,INPUT_F,sep=""))
    header<-c("chromosome", "C_Start", "C_END", "SNP_ID","Alley1","Alley2","Zero","Starnd", "chromosome2", "before_CpG", "after_CpG", "CpG_ID", "Gene", "Mean", "Var", "Zero2", "Strand2")
    colnames(CpG_SNP)<-header
    CpG_SNP$chromosome<-paste("chr",CpG_SNP$chromosome,sep="")
    CpG_SNP<-CpG_SNP[CpG_SNP$Var!=0,]
    CpG_SNP_NewMatrix<-cbind(CpG_SNP$chromosome,CpG_SNP$C_Start,CpG_SNP$C_END, as.character(CpG_SNP$SNP_ID), CpG_SNP$Zero, as.character(CpG_SNP$Starnd), CpG_SNP$before_CpG, CpG_SNP$after_CpG, as.character(CpG_SNP$CpG_ID), as.character(CpG_SNP$Gene), CpG_SNP$Mean, CpG_SNP$Var)
    header2<-c("chromosome","C_Start","C_END", "SNP_ID", "Zero", "Starnd", "before_CpG", "after_CpG", "CpG_ID", "Gene", "Mean", "Var")
    colnames(CpG_SNP_NewMatrix)<-header2
    write.table(CpG_SNP_NewMatrix, file=paste(OUTPUT_T,"out_CpG_SNP_",window_nt,"nt_filter_var.bed",sep=""), row.names=F, col.names=T, sep="\t", quote=FALSE)
}

filter_var_after_overlapSelect_sliding_window <- function(INPUT_F, OUTPUT_T){
    CpG_SNP<-read.table(paste(OUTPUT_T,INPUT_F,sep=""))
    header<-c("chromosome", "C_Start", "C_END", "SNP_ID","Alley1","Alley2","Zero","Starnd", "chromosome2", "before_CpG", "after_CpG", "CpG_ID", "Gene", "Mean", "Var", "Zero2", "Strand2")
    colnames(CpG_SNP)<-header
    CpG_SNP$chromosome<-paste("chr",CpG_SNP$chromosome,sep="")
    CpG_SNP<-CpG_SNP[CpG_SNP$Var!=0,]
    CpG_SNP_NewMatrix<-cbind(CpG_SNP$chromosome,CpG_SNP$C_Start,CpG_SNP$C_END, as.character(CpG_SNP$SNP_ID), CpG_SNP$Zero, as.character(CpG_SNP$Starnd), CpG_SNP$before_CpG, CpG_SNP$after_CpG, as.character(CpG_SNP$CpG_ID), as.character(CpG_SNP$Gene), CpG_SNP$Mean, CpG_SNP$Var)
    header2<-c("chromosome","C_Start","C_END", "SNP_ID", "Zero", "Starnd", "before_CpG", "after_CpG", "CpG_ID", "Gene", "Mean", "Var")
    CpG_SNP_NewMatrix<-data.frame(CpG_SNP_NewMatrix)
    colnames(CpG_SNP_NewMatrix)<-header2
    write.table(CpG_SNP_NewMatrix, file=paste(OUTPUT_T,"out",substring(INPUT_F,8,),sep=""), row.names=F, col.names=T, sep="\t", quote=FALSE)
}