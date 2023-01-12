annova_surv <- function(INPUT_D, OUTPUT_R, OUTPUT_T){
    
    setwd("/Users/asl061/Documents/results/BRCA/BRCA_patients/output_R/")
    
    beta_overlap <- local(get(load("matrix_CpG_T_overlap.RData")))
    SNP_overlap <- local(get(load("SNP_overlap.RData")))
    
    setwd("/Users/asl061/Documents/results/BRCA/BRCA_patients/script_R")
    source("kmplot.R")
    library(survival)
    source("/Users/asl061/Desktop/01-Projects/AMIA/tcga_integ_fun.R")
    cancer_type<-"BRCA"
    surv<-read.csv("/Users/asl061/Documents/results/BRCA/BRCA_patients/Clinical/Biotab/nationwidechildrens.org_clinical_patient_prad.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)  # important to set "stringsAsFactors=FALUSE" to modify matrix entities
    surv1<-extract_survival_info(surv)
    rownames(surv)<-surv$bcr_patient_barcode
    # substring(rownames(surv),1,4)<-cancer_type
    # substring(rownames(surv1),1,4)<-cancer_type
    
    tmp_clinic<-surv[rownames(surv1),]
    
    path_stage_intermediate<-c(5,6)
    path_stage_high<-c(8,9,10)
    
    for (i in 1:nrow(tmp_clinic)){
        if (tmp_clinic$score[i] %in% path_stage_intermediate | (tmp_clinic$score[i]=="7" & tmp_clinic$pattern_primary[i]=="3")) {
            surv1$pathologic_stage[i]<-"medium_risk"
        }    else    surv1$pathologic_stage[i]<-"high_risk"
        #     print(i)
        #     print(surv1$pathologic_stage[i])    
    }
    
    rm(tmp_clinic, surv)
    
    DATA <- "/Users/asl061/Documents/results/BRCA/BRCA_patients/overlap_file/"
    CpG_SNP_Data <- read.table(paste(DATA,"out_CpG_SNP_filter_var.bed",sep=""), header= T, sep = "\t")
    CpG_SNP_Data <- data.frame(CpG_SNP_Data)
    CpG_SNP_Data$SNP_ID <- as.character(CpG_SNP_Data$SNP_ID)
    CpG_SNP_Data$CpG_ID <- as.character(CpG_SNP_Data$CpG_ID)
    CpG_SNP_Data$Gene <- as.character(CpG_SNP_Data$Gene)
    
    matrix_CpG_T <- beta_overlap[,unique(colnames(beta_overlap)[(colnames(beta_overlap) %in% colnames(SNP_overlap))])]
    
    Tumour_SNPs <- NULL;snpids1 <- NULL; CpG_id1 <- NULL; geneid1 <- NULL; freq1 <- NULL; freq2 <- NULL; freq3 <- NULL; out3 <-NULL
    out4 <- NULL; out5 <- NULL
    
    # plot related
    group_label<-c('1/1','1/0','0/0')
    
    alternative<-1; aeterozygous<-2; reference<-3
    for (i in 1:nrow(matrix_CpG_T)){ 
        snpids<-CpG_SNP_Data$SNP_ID[grep(rownames(matrix_CpG_T)[i],CpG_SNP_Data$CpG_ID)]
        
        for(j in 1:length(snpids)){
            # SNP_overlap: genotype; matrix_CpG_T: beta value
            m1<-rbind(SNP_overlap[snpids[j],colnames(matrix_CpG_T)], matrix_CpG_T[i,]) 
            
            Tumour_SNPs <- rbind(m1, Tumour_SNPs)
            m2 <- t(m1)
            colnames(m2)<-c("genotype_group","beta")
            m2 <- data.frame(m2)
            
            if (length(table(m2$genotype_group)) > 2){  # for anova, it requires three genotype_groups (1,2,3); alternative, heterogeneous, observed)
                m2$genotype_group <- as.factor(m2$genotype_group)
                rs <- anova(aov(beta ~ genotype_group,data=m2 ))[1,"Pr(>F)"] # aov: fit an analysis of variance model
                surv2 <- surv1[rownames(surv1) %in% rownames(m2),]   
                surv2$genotype_group <- m2[rownames(surv2),1]
                fit <- survfit(Surv(days_to_death,vital_status)~genotype_group, data=surv2)
                
                sdiff <- survdiff(eval(fit$call$formula), data = eval(fit$call$data))
                pval <- pchisq(sdiff$chisq,length(sdiff$n) - 1,lower.tail = FALSE) # chi-squared distribution
                pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))
                
                print(i); print(rs); print(pval)
                ## STOP HERE
                geneid<-CpG_SNP_Data$Gene[grep(rownames(matrix_CpG_T)[i],CpG_SNP_Data$CpG_ID)]
                CpG_id <- rownames(matrix_CpG_T)[i]
                # genotype frequency
                genotype_freq1 <- rbind(sum(m1[1,]==alternative, genotype_freq1); genotype_freq2 <- rbind(sum(m1[1,]==aeterozygous), genotype_freq2); genotype_freq3 <- rbind(sum(m1[1,]==reference), genotype_freq3)
                                        
                                        CpG_id1 <- rbind(CpG_id, CpG_id1); snpids1 <- rbind(snpids[j], snpids1); geneid1 <- rbind(geneid[j], geneid1)
                                        out3 <- rbind(i,out3); out4 <- rbind(rs, out4); out5 <- rbind(pval, out5)
                                        
                                        output <- cbind( out3, out4, out5, CpG_id1, snpids1, geneid1, genotype_freq1, genotype_freq2, genotype_freq3)
                                        colnames(output) <- c("S.No", "beta-p", "surv_p", "CpG_ID", "SNP_ID", "GENE", "AA", "AB","BB" ) 
            }
            
            #         filename = paste(Output,"/","survival_slhuet_",cancer_type,".png",sep="")
            if (out4 <0.05 & out5 <0.05){
                filename = paste(snpids1,"_survival",".png",sep="")
                png(file=filename,width = 1024, height = 768)
                kmplot(fit, mark='',
                       #        xaxis.at=c(0,0.5,1:9)*365, xaxis.lab=c(0,0.5,1:9), # n.risk.at
                       xaxis.at=c(0,1:9)*365, xaxis.lab=c(0,1:9), # n.risk.at
                       #        lty.surv=c(1,2), lwd.surv=1, col.surv=c(1,2), # survival.curves
                       lty.surv=c(1), lwd.surv=3, col.surv=c(1,2,3,4), # survival.curves
                       lty.ci=0,        lwd.ci=1,   col.ci=1, # confidence intervals not plotted 
                       #        group.names=c('C1','C2','C3'),
                       group.names=group_label,
                       group.order=seq(1:3), # order of appearance in the n.risk.at table and legend.
                       #        group.order=c(1,2,3), # order of appearance in the n.risk.at table and legend.
                       extra.left.margin=6, label.n.at.risk=FALSE, draw.lines=TRUE, 
                       cex.axis=1.5, xlab='Years', ylab='Survival Probability', # labels
                       grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
                       legend=TRUE, loc.legend='bottomleft',
                       cex.lab=1.5, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
                )
                mtext(snpids1, side= 3,  outer = TRUE)
                dev.off()       
            }
        }
    }
    output<-output[,-3]; output<-data.frame(output)
    sig_output<-output[as.numeric(as.character(output$beta.p)) <0.05 & as.numeric(as.character(output$surv_p)) <0.05,]
    write.table(sig_output,file="sig_output.txt",sep="\t")
}