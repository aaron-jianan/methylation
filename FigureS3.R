##Import original file
load("Figure2_Raw.Rdata")
library(ggDCA)
library(rms)
library(ggplot2)
library(tidyverse)

##Plotting the DCA curve of CCP binding to CpGs
colnames(data_raw_Roc_RA)
data_raw_Roc_RA$group=as.numeric(data_raw_Roc_RA$group)
data_raw_Roc_RA$group=ifelse(str_detect(data_raw_Roc_RA$group,"2"),1,0)
class(data_raw_Roc_RA$group)

colnames(data_raw_Roc_RA)[1:7]=paste0("CpG_",colnames(data_raw_Roc_RA)[1:7])

CCP_cg15692052 <- lrm(group ~ cg15692052 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898116 <- lrm(group ~ CpG_46898116 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898066<- lrm(group ~ CpG_46898066 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898048 <- lrm(group ~ CpG_46898048 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898042 <- lrm(group ~ CpG_46898042 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898024 <- lrm(group ~ CpG_46898024 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898006 <- lrm(group ~ CpG_46898006 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898004 <- lrm(group ~ CpG_46898004 + CCP,data_raw_Roc_RA,maxit=1000)

fit_DCA_all<- dca(CCP_46898116,CCP_46898066,CCP_46898048,
                  CCP_46898042,CCP_46898024,CCP_46898006,
                  CCP_46898004,CCP_cg15692052)

P107=ggplot(fit_DCA_all)

##Plotting the DCA curve of RF binding to CpGs
RF_cg15692052 <- lrm(group ~ cg15692052 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898116 <- lrm(group ~ CpG_46898116 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898066<- lrm(group ~ CpG_46898066 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898048 <- lrm(group ~ CpG_46898048 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898042 <- lrm(group ~ CpG_46898042 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898024 <- lrm(group ~ CpG_46898024 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898006 <- lrm(group ~ CpG_46898006 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898004 <- lrm(group ~ CpG_46898004 + RF,data_raw_Roc_RA,maxit=1000)
fit_DCA_all_RF<- dca(RF_46898116,RF_46898066,RF_46898048,
                     RF_46898042,RF_46898024,RF_46898006,
                     RF_46898004,RF_cg15692052)

P109=ggplot(fit_DCA_all_RF)


##Plotting the DCA curve of CCP binding to haplotype
load("Figure3_FigureS2_Raw.Rdata")
Haplotype_data_RA$group=as.numeric(Haplotype_data_RA$group)
Haplotype_data_RA$group=ifelse(str_detect(Haplotype_data_RA$group,"2"),1,0)

CCCCCCC_CCP <- lrm(group ~ CCCCCCC + CCP,Haplotype_data_RA,maxit=1000)
TCCCCCC_CCP<- lrm(group ~ TCCCCCC + CCP ,Haplotype_data_RA,maxit=1000)
TTTTTTT_CCP<- lrm(group ~ TTTTTTT + CCP ,Haplotype_data_RA,maxit=1000)
TTCCCCC_CCP<- lrm(group ~ TTCCCCC + CCP ,Haplotype_data_RA,maxit=1000)
CCCTCCC_CCP<- lrm(group ~ CCCTCCC + CCP ,Haplotype_data_RA,maxit=1000)
TTTTTCC_CCP<- lrm(group ~ TTTTTCC + CCP ,Haplotype_data_RA,maxit=1000)
CTCCCCC_CCP<- lrm(group ~ CTCCCCC + CCP ,Haplotype_data_RA,maxit=1000)
CCTCCCC_CCP<- lrm(group ~ CCTCCCC + CCP ,Haplotype_data_RA,maxit=1000)
TTTTCCC_CCP<- lrm(group ~ TTTTCCC + CCP ,Haplotype_data_RA,maxit=1000)
TCCTCCC_CCP<- lrm(group ~ TCCTCCC + CCP ,Haplotype_data_RA,maxit=1000)
TCTCCCC_CCP<- lrm(group ~ TCTCCCC + CCP ,Haplotype_data_RA,maxit=1000)
TTTTTTC_CCP<- lrm(group ~ TTTTTTC + CCP ,Haplotype_data_RA,maxit=1000)
CCCCCCT_CCP<- lrm(group ~ CCCCCCT + CCP ,Haplotype_data_RA,maxit=1000)
CCCCCTC_CCP<- lrm(group ~ CCCCCTC + CCP ,Haplotype_data_RA,maxit=1000)
TTTCTCC_CCP<- lrm(group ~ TTTCTCC + CCP ,Haplotype_data_RA,maxit=1000)


fit_DCA_all_1<- dca(CCCCCCC_CCP,
                    TCCCCCC_CCP,
                    TTTTTTT_CCP,
                    TTCCCCC_CCP,
                    CCCTCCC_CCP,
                    TTTTTCC_CCP,
                    CTCCCCC_CCP,
                    CCTCCCC_CCP)


fit_DCA_all_2=dca( TTTTCCC_CCP,
                   TCCTCCC_CCP,
                   TCTCCCC_CCP,
                   TTTTTTC_CCP,
                   CCCCCCT_CCP,
                   CCCCCTC_CCP,
                   TTTCTCC_CCP)

P110=ggplot(fit_DCA_all_1)

P111=ggplot(fit_DCA_all_2)


##Plotting the DCA curve of RF binding to haplotype
CCCCCCC_RF <- lrm(group ~ CCCCCCC + RF,Haplotype_data_RA,maxit=1000)
TCCCCCC_RF<- lrm(group ~ TCCCCCC + RF ,Haplotype_data_RA,maxit=1000)
TTTTTTT_RF<- lrm(group ~ TTTTTTT +RF ,Haplotype_data_RA,maxit=1000)
TTCCCCC_RF<- lrm(group ~ TTCCCCC + RF ,Haplotype_data_RA,maxit=1000)
CCCTCCC_RF<- lrm(group ~ CCCTCCC + RF ,Haplotype_data_RA,maxit=1000)
TTTTTCC_RF<- lrm(group ~ TTTTTCC + RF ,Haplotype_data_RA,maxit=1000)
CTCCCCC_RF<- lrm(group ~ CTCCCCC + RF ,Haplotype_data_RA,maxit=1000)
CCTCCCC_RF<- lrm(group ~ CCTCCCC + RF ,Haplotype_data_RA,maxit=1000)
TTTTCCC_RF<- lrm(group ~ TTTTCCC + RF ,Haplotype_data_RA,maxit=1000)
TCCTCCC_RF<- lrm(group ~ TCCTCCC + RF ,Haplotype_data_RA,maxit=1000)
TCTCCCC_RF<- lrm(group ~ TCTCCCC + RF ,Haplotype_data_RA,maxit=1000)
TTTTTTC_RF<- lrm(group ~ TTTTTTC + RF ,Haplotype_data_RA,maxit=1000)
CCCCCCT_RF<- lrm(group ~ CCCCCCT + RF ,Haplotype_data_RA,maxit=1000)
CCCCCTC_RF<- lrm(group ~ CCCCCTC + RF ,Haplotype_data_RA,maxit=1000)
TTTCTCC_RF<- lrm(group ~ TTTCTCC + RF ,Haplotype_data_RA,maxit=1000)

fit_DCA_all_3<- dca(CCCCCCC_RF,
                    TCCCCCC_RF,
                    TTTTTTT_RF,
                    TTCCCCC_RF,
                    CCCTCCC_RF,
                    TTTTTCC_RF,
                    CTCCCCC_RF,
                    CCTCCCC_RF)


fit_DCA_all_4=dca( TTTTCCC_RF,
                   TCCTCCC_RF,
                   TCTCCCC_RF,
                   TTTTTTC_RF,
                   CCCCCCT_RF,
                   CCCCCTC_RF,
                   TTTCTCC_RF)

P112=ggplot(fit_DCA_all_3)

P113=ggplot(fit_DCA_all_4)


library(patchwork)
P107+P109+P110+P111+P112+P113
ggsave(filename = "FigureS3.pdf",width = 15,height = 8)

save.image(file = "FigureS3.Rdata")

