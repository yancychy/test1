## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  Version: 0.1                                            ##
##  Date: 2018-05-09                                        ##
##  Title: variaton/methylation/risk analysis               ##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## -------- user-defined HOME directory: modify below to your home directory --------
HOME <- "/Users/ASH/Documents/NGS/BRCA/SNP_CpG_Surv_Analysis_SlidingWindow_5k/"

## -------- directory definition ----------------------------------------------------
#  SCRIPT: source code;     OUTPUT_T: text based output (e.g., txt, bed)
#  OUTPUT_R: RData output;  OUTPUT_F: figure output

#  creat folder if doesn't exsit and set the path; transfer data to input_D folder (SNP, DNA methylation, etc)
# dir.create(paste(HOME,"output_R",sep="")); dir.create(paste(HOME,"output_T",sep="")) 
# dir.create(paste(HOME,"output_F",sep="")); dir.create(paste(HOME,"input_D",sep="")) 
SCRIPT<-paste(HOME,"script_R/",sep=""); OUTPUT_T<-paste(HOME,"output_T/",sep=""); 
OUTPUT_R<-paste(HOME,"output_R/",sep=""); OUTPUT_F<-paste(HOME,"output_F/",sep=""); 
INPUT_D<-paste("/Users/ASH/Documents/NGS/BRCA/SNP_CpG_input_D/",sep="")

source(paste(SCRIPT,"combine_data.R", sep="")) 
source(paste(SCRIPT,"bedfile_generator.R", sep=""))
source(paste(SCRIPT,"filter_var_after_overlapSelect.R", sep=""))
source(paste(SCRIPT,"extract_data.R", sep=""))
source(paste(SCRIPT,"annova_surv.R", sep=""))


# extract all beta values for DNA methylation
# combine_all_DNA_methyl_data(INPUT_D, OUTPUT_R, OUTPUT_T)
# combine_all_SNP_data(HOME, INPUT_D, OUTPUT_R, OUTPUT_T)

# generate bed file
window_nt<-100; overlap_nt<-50; max_boundary<-5000
# SNP_bedFile_Gen(OUTPUT_T, window_nt)
# DNAmethyl_bedFile_sliding_window_Gen(INPUT_D, OUTPUT_T, window_nt, overlap_nt, max_boundary)

overlapSelect<-paste(SCRIPT,"overlapSelect_exec.sh", sep="")
for (i in 1:length(seq(0,max_boundary*2,overlap_nt))){ 
    DNAmethy_bedfile<-paste(OUTPUT_T,"DNAmeth_",max_boundary*2,"nt_",window_nt,"nt_",i,".bed" ,sep="")
    SNP_bedfile<-paste(OUTPUT_T,"SNP_",window_nt,"nt.bed",sep="")
    output_bedfile<-paste(OUTPUT_T,"overlap_CpG_SNP_",max_boundary*2,"nt_",window_nt,"nt_",i,".bed" ,sep="")
    command<-paste(overlapSelect,DNAmethy_bedfile, SNP_bedfile, output_bedfile)
    system(command)
    INPUT_F<-paste("overlap_CpG_SNP_",max_boundary*2,"nt_",window_nt,"nt_",i,".bed" ,sep="")
    filter_var_after_overlapSelect_sliding_window(INPUT_F, OUTPUT_T)  ## filter here
    extract_methyl_data_sliding_window(INPUT_D, OUTPUT_R, OUTPUT_T, max_boundary, window_nt, i)
    extract_SNP_data_sliding_window(INPUT_D, OUTPUT_R, OUTPUT_T, max_boundary, window_nt, i)
}



