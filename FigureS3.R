##Import original file
load("FigureS3_raw_data1.Rdata")
library(ggDCA)
library(rms)
library(ggplot2)
library(tidyverse)
class(data_raw_Roc_RA$group)
data_raw_Roc_RA$group=as.numeric(data_raw_Roc_RA$group)
data_raw_Roc_RA$group=ifelse(str_detect(data_raw_Roc_RA$group,"2"),1,0)
colnames(data_raw_Roc_RA)
data_raw_Roc_RA=data_raw_Roc_RA[,c(1:8,21:22)]
data_raw_Roc_RA=cbind(data_raw_Roc_RA,mean_data)

CCP_cg15692052 <- lrm(group ~ cg15692052_35 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898116 <- lrm(group ~ cg15692052_35_75 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898066<- lrm(group ~ cg15692052_35_125 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898048 <- lrm(group ~ cg15692052_35_143 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898042 <- lrm(group ~ cg15692052_35_149 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898024 <- lrm(group ~ cg15692052_35_167 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898006 <- lrm(group ~ cg15692052_35_185 + CCP,data_raw_Roc_RA,maxit=1000)
CCP_46898004 <- lrm(group ~ cg15692052_35_187 + CCP,data_raw_Roc_RA,maxit=1000)

fit_DCA_075 <- dca(cg15692052_35_75_CCP)
fit_DCA_0125 <- dca(cg15692052_35_125_CCP)
fit_DCA_0143 <- dca(cg15692052_35_143_CCP)
fit_DCA_0149 <- dca(cg15692052_35_149_CCP)
fit_DCA_0167 <- dca(cg15692052_35_167_CCP)
fit_DCA_0185 <- dca(cg15692052_35_185_CCP)
fit_DCA_0187 <- dca(cg15692052_35_187_CCP)

write.csv(fit_DCA_075,file ="075_DCA_CCP_results.csv" )
write.csv(fit_DCA_0125,file ="0125_DCA_CCP_results.csv" )
write.csv(fit_DCA_0143,file ="0143_DCA_CCP_results.csv" )
write.csv(fit_DCA_0149,file ="0149_DCA_CCP_results.csv" )
write.csv(fit_DCA_0167,file ="0167_DCA_CCP_results.csv" )
write.csv(fit_DCA_0185,file ="0185_DCA_CCP_results.csv" )
write.csv(fit_DCA_0187,file ="0187_DCA_CCP_results.csv" )
fit_DCA_all<- dca(CCP_46898116,CCP_46898066,CCP_46898048,CCP_46898042,
                  CCP_46898024,CCP_46898006,CCP_46898004,CCP_cg15692052)

P1=ggplot(fit_DCA_all)

##RF
RF_cg15692052 <- lrm(group ~ cg15692052_35 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898116 <- lrm(group ~ cg15692052_35_75 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898066<- lrm(group ~ cg15692052_35_125 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898048 <- lrm(group ~ cg15692052_35_143 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898042 <- lrm(group ~ cg15692052_35_149 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898024 <- lrm(group ~ cg15692052_35_167 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898006 <- lrm(group ~ cg15692052_35_185 + RF,data_raw_Roc_RA,maxit=1000)
RF_46898004 <- lrm(group ~ cg15692052_35_187 + RF,data_raw_Roc_RA,maxit=1000)

##draw all
fit_DCA_all_RF<- dca(RF_46898116,RF_46898066,RF_46898048,RF_46898042,
                     RF_46898024,RF_46898006,RF_46898004,RF_cg15692052)

P2=ggplot(fit_DCA_all_RF)
P2
rm(BMI_data,combine_data_raw_BMI,combine_data_raw_merge_data,cor_matr_BMI_Hmisc,
   cor_matr_merge_Hmisc)
save.image(file = "20220808_CCP_RA_DCA_RAW.Rdata")
load("20220808_CCP_RA_DCA_RAW.Rdata")

##haplotype_CCP
load("Haplotype_DCA_RAW.Rdata")
data_raw_Roc_RA=data_raw_Roc_RA[,c(1:16,29:30)]
colnames(data_raw_Roc_RA)[1:15]=toupper(colnames(data_raw_Roc_RA)[1:15])
data_raw_Roc_RA$group=ifelse(str_detect(data_raw_Roc_RA$group,"0"),1,0)
class(data_raw_Roc_RA$cg15692052_35)
data_raw_Roc_RA=cbind(data_raw_Roc_RA,mean_data)


CCCCCCC_CCP <- lrm(group ~ CCCCCCC + CCP,data_raw_Roc_RA,maxit=1000)
TCCCCCC_CCP<- lrm(group ~ TCCCCCC + CCP ,data_raw_Roc_RA,maxit=1000)
TTTTTTT_CCP<- lrm(group ~ TTTTTTT + CCP ,data_raw_Roc_RA,maxit=1000)
TTCCCCC_CCP<- lrm(group ~ TTCCCCC + CCP ,data_raw_Roc_RA,maxit=1000)
CCCTCCC_CCP<- lrm(group ~ CCCTCCC + CCP ,data_raw_Roc_RA,maxit=1000)
TTTTTCC_CCP<- lrm(group ~ TTTTTCC + CCP ,data_raw_Roc_RA,maxit=1000)
CTCCCCC_CCP<- lrm(group ~ CTCCCCC + CCP ,data_raw_Roc_RA,maxit=1000)
CCTCCCC_CCP<- lrm(group ~ CCTCCCC + CCP ,data_raw_Roc_RA,maxit=1000)
TTTTCCC_CCP<- lrm(group ~ TTTTCCC + CCP ,data_raw_Roc_RA,maxit=1000)
TCCTCCC_CCP<- lrm(group ~ TCCTCCC + CCP ,data_raw_Roc_RA,maxit=1000)
TCTCCCC_CCP<- lrm(group ~ TCTCCCC + CCP ,data_raw_Roc_RA,maxit=1000)
TTTTTTC_CCP<- lrm(group ~ TTTTTTC + CCP ,data_raw_Roc_RA,maxit=1000)
CCCCCCT_CCP<- lrm(group ~ CCCCCCT + CCP ,data_raw_Roc_RA,maxit=1000)
CCCCCTC_CCP<- lrm(group ~ CCCCCTC + CCP ,data_raw_Roc_RA,maxit=1000)
TTTCTCC_CCP<- lrm(group ~ TTTCTCC + CCP ,data_raw_Roc_RA,maxit=1000)

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

P3=ggplot(fit_DCA_all_1)
P3
P4=ggplot(fit_DCA_all_2)
P4

##haplotype_RF
CCCCCCC_RF <- lrm(group ~ CCCCCCC + RF,data_raw_Roc_RA,maxit=1000)
TCCCCCC_RF<- lrm(group ~ TCCCCCC + RF ,data_raw_Roc_RA,maxit=1000)
TTTTTTT_RF<- lrm(group ~ TTTTTTT +RF ,data_raw_Roc_RA,maxit=1000)
TTCCCCC_RF<- lrm(group ~ TTCCCCC + RF ,data_raw_Roc_RA,maxit=1000)
CCCTCCC_RF<- lrm(group ~ CCCTCCC + RF ,data_raw_Roc_RA,maxit=1000)
TTTTTCC_RF<- lrm(group ~ TTTTTCC + RF ,data_raw_Roc_RA,maxit=1000)
CTCCCCC_RF<- lrm(group ~ CTCCCCC + RF ,data_raw_Roc_RA,maxit=1000)
CCTCCCC_RF<- lrm(group ~ CCTCCCC + RF ,data_raw_Roc_RA,maxit=1000)
TTTTCCC_RF<- lrm(group ~ TTTTCCC + RF ,data_raw_Roc_RA,maxit=1000)
TCCTCCC_RF<- lrm(group ~ TCCTCCC + RF ,data_raw_Roc_RA,maxit=1000)
TCTCCCC_RF<- lrm(group ~ TCTCCCC + RF ,data_raw_Roc_RA,maxit=1000)
TTTTTTC_RF<- lrm(group ~ TTTTTTC + RF ,data_raw_Roc_RA,maxit=1000)
CCCCCCT_RF<- lrm(group ~ CCCCCCT + RF ,data_raw_Roc_RA,maxit=1000)
CCCCCTC_RF<- lrm(group ~ CCCCCTC + RF ,data_raw_Roc_RA,maxit=1000)
TTTCTCC_RF<- lrm(group ~ TTTCTCC + RF ,data_raw_Roc_RA,maxit=1000)

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

P5=ggplot(fit_DCA_all_3)
P5
P6=ggplot(fit_DCA_all_4)
P6

library(patchwork)
P1+P2+P3+P4+P5+P6
ggsave(filename = "ccp_rf_DCA.pdf",width = 15,height = 8)

save.image(file = "FigureS3_raw_data2.Rdata")

