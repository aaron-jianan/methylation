####Read the raw data and process it, it need to use the raw data data_raw_Roc_RA of Figure 2
load("Figure2_Raw.Rdata")
library(pROC)
library(ggplot2)
library(tidyverse)

##1. Different CpG diagnoses the ROC calculation of multiple RA types, after outputting the calculation results, manually merge the 5 RA types

##Calculation of AUC for RA diagnosis alone for different CpGs
{ROCStatFunc <- function(dat, group, var,retype = c("threshold", "specificity", "sensitivity"),
                        auc = T,youden = T, digit = 3){
  subgroup <- levels(as.factor(dat[[group]]))
  subgroup1 <- paste0(subgroup[2], " vs ", subgroup[1])
  rocmodel <- roc(dat[[group]], dat[[var]])
  other <- coords(rocmodel, "b", ret = retype)
  other <- round(other, digit)
  if(auc == T){
    auc <- round(ci.auc(rocmodel),digit)
    auc <- paste0(auc[2],"(",auc[1],"-",auc[3],")")
    if(youden == T){
      abc <- coords(rocmodel, "b", ret = c("specificity", "sensitivity"))
      youdenres <- abc[1] + abc[2] - 1
      youdenres <- round(youdenres, digit)
      result <- c(group, subgroup1, auc, other, youdenres)
      names(result) <- c("group", "subgroup","auc(95%CI)", retype, "youden")
    }else{
      result <- c(group, subgroup1, auc, other)
      names(result) <- c("group", "subgroup", "auc(95%CI)", retype)
    }
  }else{
    if(youden == T){
      abc <- coords(rocmodel, "b", ret = c("specificity", "sensitivity"))
      youdenres <- abc[1] + abc[2] - 1
      youdenres <- round(youdenres, digit)
      result <- c(group, subgroup1, other, youdenres)
      names(result) <- c("group","subgroup", retype, "youden")
    }else{
      result <- c(group, subgroup1,other)
      names(result) <- c("group", "subgroup",retype)
    }
  }
  return(result)
}

quiteROCFunc <- quietly(ROCStatFunc)
multigroup <- colnames(data_raw_Roc_RA)[1:8]
rocRes <- lapply(multigroup, function(x) quiteROCFunc(data_raw_Roc_RA, "group", x)$result)
rocResDat <- do.call(rbind, rocRes)
rocResDat=as.data.frame(rocResDat)
rocResDat$variable=multigroup 
rocResDat<- apply(rocResDat,2,as.character)
write.csv(rocResDat,file ="Figure3_raw_data1_RA.csv")}

##AUC of CpGs alone to diagnose 4 other types of RA
RF_CCP_seprate_data=data_raw_Roc_RA[1:12]
RF_CCP_seprate_data=RF_CCP_seprate_data[RF_CCP_seprate_data$group==1,]

##RF-/CCP-
{RF_CCP_seprate_data$group=ifelse(RF_CCP_seprate_data$RF<20 & RF_CCP_seprate_data$CCP<25,1,0)
RF_CCP_seprate_data$group=as.factor(RF_CCP_seprate_data$group)
rocRes <- lapply(multigroup, function(x) quiteROCFunc(RF_CCP_seprate_data, "group", x)$result)
rocResDat <- do.call(rbind, rocRes)
rocResDat
rocResDat=as.data.frame(rocResDat)
rocResDat$variable=multigroup 
rocResDat<- apply(rocResDat,2,as.character)
write.csv(rocResDat,file ="Figure3_raw_data2_RF_neg_CCP_neg_ROC.csv")}

##RF-/CCP+
{
  RF_negtive_seprate_data=RF_CCP_seprate_data
  RF_negtive_seprate_data$group=ifelse(RF_negtive_seprate_data$RF<20&RF_negtive_seprate_data$CCP>25,1,0)
  RF_negtive_seprate_data$group=as.factor(RF_negtive_seprate_data$group)
  rocRes1 <- lapply(multigroup, function(x) quiteROCFunc(RF_negtive_seprate_data, "group", x)$result)
  rocResDat1 <- do.call(rbind, rocRes1)
  rocResDat1
  rocResDat1=as.data.frame(rocResDat1)
  rocResDat1$variable=multigroup 
  rocResDat1<- apply(rocResDat1,2,as.character)
  write.csv(rocResDat1,file ="Figure3_raw_data3_RF_neg_CCP_pos_ROC.csv")
}

##RF+/CCP-
{
  CCP_negtive_seprate_data=RF_CCP_seprate_data
  CCP_negtive_seprate_data$group=ifelse( CCP_negtive_seprate_data$RF>20& CCP_negtive_seprate_data$CCP<25,1,0)
  CCP_negtive_seprate_data$group=as.factor(CCP_negtive_seprate_data$group)
  rocRes2 <- lapply(multigroup, function(x) quiteROCFunc(CCP_negtive_seprate_data, "group", x)$result)
  rocResDat2 <- do.call(rbind, rocRes2)
  rocResDat2
  rocResDat2=as.data.frame(rocResDat2)
  rocResDat2$variable=multigroup 
  rocResDat2<- apply(rocResDat2,2,as.character)
  write.csv(rocResDat2,file ="Figure3_raw_data4_RF_pos_CCP_neg_ROC.csv")
}

##RF+/CCP+
{ 
  RF_CCP_postive_seprate_data=RF_CCP_seprate_data
  RF_CCP_postive_seprate_data$group=ifelse( RF_CCP_postive_seprate_data$RF>20& RF_CCP_postive_seprate_data$CCP>25,1,0)
  RF_CCP_postive_seprate_data$group=as.factor(RF_CCP_postive_seprate_data$group)
  rocRes7 <- lapply(multigroup, function(x) quiteROCFunc(RF_CCP_postive_seprate_data, "group", x)$result)
  rocResDat7 <- do.call(rbind, rocRes7)
  rocResDat7
  rocResDat7=as.data.frame(rocResDat7)
  rocResDat7$variable=multigroup 
  rocResDat7<- apply(rocResDat7,2,as.character)
  write.csv(rocResDat7,file ="Figure3_raw_data5_RF_pos_CCP_pos_ROC.csv")}


##2. Different haplotypes diagnose the ROC calculation of multiple RA types. After outputting the calculation results, manually merge the five RA types.

##Data read and processed raw
library(readxl)
Methylation_raw_data <- read_excel("Methylation_raw_data.xlsx", 
                                   sheet = "Haplotype")
Haplotype_data_RA=Methylation_raw_data[,-c(1,3)]
Haplotype_name=Haplotype_data_RA$Haplotype
Haplotype_ID=colnames(Haplotype_data_RA)[2:302]
Haplotype_data_RA=as.data.frame(t(Haplotype_data_RA))
colnames(Haplotype_data_RA)=toupper(Haplotype_name)
Haplotype_data_RA=Haplotype_data_RA[-1,]
Haplotype_data_RA=apply(Haplotype_data_RA, 2, function(x){as.numeric(x)})
rownames(Haplotype_data_RA)=Haplotype_ID
Haplotype_data_RA=as.data.frame(Haplotype_data_RA)

##Output raw data of haplotype ratios for subsequent use
write.csv(Haplotype_data_RA,"Haplotype_raw_data.csv")

##Calculation of AUC for RA diagnosis alone for different Haplotype
{Haplotype_data_RA$group=ifelse(str_detect(rownames(Haplotype_data_RA),"HC"),"HC",
                              ifelse(str_detect(rownames(Haplotype_data_RA),"OA"),"OA","RA"))
Haplotype_data_RA$group=ifelse(str_detect(Haplotype_data_RA$group,"HC"),"0",
                             ifelse(str_detect(Haplotype_data_RA$group,"OA"),"0","1"))
Haplotype_data_RA$group=as.numeric(Haplotype_data_RA$group)
Haplotype_data_RA$group=as.factor(Haplotype_data_RA$group)

multigroup <- colnames(Haplotype_data_RA)[1:15]
rocRes <- lapply(multigroup, function(x) quiteROCFunc(Haplotype_data_RA, "group", x)$result)
rocResDat <- do.call(rbind, rocRes)
rocResDat=as.data.frame(rocResDat)
rocResDat$variable=multigroup 
rocResDat<- apply(rocResDat,2,as.character)
write.csv(rocResDat,file ="Figure3_raw_data6_Haplotype_RA.csv")}

##RF-/CCP-
identical(rownames(Haplotype_data_RA),rownames(data_raw_Roc_RA))
Haplotype_data_RA$RF=data_raw_Roc_RA$RF
Haplotype_data_RA$CCP=data_raw_Roc_RA$CCP
Haplotype_RF_CCP_seprate_rawdata=Haplotype_data_RA[Haplotype_data_RA$group==1,]
{Haplotype_RF_CCP_seprate_rawdata$group=ifelse(Haplotype_RF_CCP_seprate_rawdata$RF<20 & Haplotype_RF_CCP_seprate_rawdata$CCP<25,1,0)
  Haplotype_RF_CCP_seprate_rawdata$group=as.factor(Haplotype_RF_CCP_seprate_rawdata$group)
  table(danbeixing_RF_CCP_seprate_rawdata$group)
  rocRes <- lapply(multigroup, function(x) quiteROCFunc(Haplotype_RF_CCP_seprate_rawdata, "group", x)$result)
  rocResDat <- do.call(rbind, rocRes)
  rocResDat
  rocResDat=as.data.frame(rocResDat)
  rocResDat$variable=multigroup 
  rocResDat<- apply(rocResDat,2,as.character)
  write.csv(rocResDat,file ="Figure3_raw_data7_Haplotype_RF_neg_CCP_neg.csv")
}

##RF-/CCP+
{
  Haplotype_RF_negtive_data=Haplotype_RF_CCP_seprate_rawdata
  Haplotype_RF_negtive_data$group=ifelse(Haplotype_RF_negtive_data$RF<20 & Haplotype_RF_negtive_data$CCP>25,1,0)
  Haplotype_RF_negtive_data$group=as.factor(Haplotype_RF_negtive_data$group)
  rocRes3 <- lapply(multigroup, function(x) quiteROCFunc(Haplotype_RF_negtive_data, "group", x)$result)
  rocResDat3 <- do.call(rbind, rocRes3)
  rocResDat3
  rocResDat3=as.data.frame(rocResDat3)
  rocResDat3$variable=multigroup 
  rocResDat3<- apply(rocResDat3,2,as.character)
  write.csv(rocResDat3,file ="Figure3_raw_data8_Haplotype_RF_neg_CCP_pos.csv")
}

##RF+/CCP-
{
  Haplotype_CCP_negtive_data=Haplotype_RF_CCP_seprate_rawdata
  Haplotype_CCP_negtive_data$group=ifelse(Haplotype_CCP_negtive_data$RF>20 & Haplotype_CCP_negtive_data$CCP<25,1,0)
  Haplotype_CCP_negtive_data$group=as.factor(Haplotype_CCP_negtive_data$group)
  rocRes4 <- lapply(multigroup, function(x) quiteROCFunc(Haplotype_CCP_negtive_data, "group", x)$result)
  rocResDat4 <- do.call(rbind, rocRes4)
  rocResDat4
  rocResDat4=as.data.frame(rocResDat4)
  rocResDat4$variable=multigroup 
  rocResDat4<- apply(rocResDat4,2,as.character)
  write.csv(rocResDat4,file ="Figure3_raw_data9_Haplotype_RF_pos_CCP_neg.csv")
}

##RF+/CCP+
{
  Haplotype_RF_CCP_positive_data=Haplotype_RF_CCP_seprate_rawdata
  Haplotype_RF_CCP_positive_data$group=ifelse(Haplotype_RF_CCP_positive_data$RF>20 & Haplotype_RF_CCP_positive_data$CCP>25,1,0)
  Haplotype_RF_CCP_positive_data$group=as.factor(Haplotype_RF_CCP_positive_data$group)
  rocRes5 <- lapply(multigroup, function(x) quiteROCFunc( Haplotype_RF_CCP_positive_data, "group", x)$result)
  rocResDat5 <- do.call(rbind, rocRes5)
  rocResDat5
  rocResDat5=as.data.frame(rocResDat5)
  rocResDat5$variable=multigroup 
  rocResDat5<- apply(rocResDat5,2,as.character)
  write.csv(rocResDat5,file ="Figure3_raw_data10_Haplotype_RF_pos_CCP_pos.csv")
}

##3. ROC calculation for 3 RA comorbidities diagnosed by different CpGs, after outputting the calculation results, manually merge

##Process raw data
RA_comorbidities_data=cbind(data_raw_Roc_RA[,c(1:8,12)],clinical_meta[,c(18:20)])

##Interstitial lung disease
{multigroup <- colnames(RA_comorbidities_data)[c(1:8)]
rocRes <- lapply(multigroup, function(x) quiteROCFunc(RA_comorbidities_data, "lung (yes 1, no 0)", x)$result)
rocResDat <- do.call(rbind, rocRes)
rocResDat=as.data.frame(rocResDat)
rocResDat$variable=multigroup 
rocResDat<- apply(rocResDat,2,as.character)
write.csv(rocResDat,file ="Figure3_raw_data11_lung_ROC_ALL.csv")}

##high_blood_pressure
{rocRes <- lapply(multigroup, function(x) quiteROCFunc(RA_comorbidities_data, "hyper (yes 1, no 0)", x)$result)
rocResDat <- do.call(rbind, rocRes)
rocResDat
rocResDat=as.data.frame(rocResDat)
rocResDat$variable=multigroup 
rocResDat<- apply(rocResDat,2,as.character)
write.csv(rocResDat,file ="Figure3_raw_data12_high_blood_pressure_ROC_ALL.csv")}

##osteoporosis
{rocRes <- lapply(multigroup, function(x) quiteROCFunc(RA_comorbidities_data, "bone (yes 1, no 0)", x)$result)
rocResDat <- do.call(rbind, rocRes)
rocResDat=as.data.frame(rocResDat)
rocResDat$variable=multigroup 
rocResDat<- apply(rocResDat,2,as.character)
write.csv(rocResDat,file ="Figure3_raw_data13_osteoporosis_ROC_ALL.csv")}


##4. The ROC calculation of different haplotypes for diagnosing three RA comorbidities, after outputting the calculation results, manually merge

##Process raw data
RA_comorbidities_haplotype_data=cbind(Haplotype_data_RA[,c(1:15)],clinical_meta[,c(18:20)])

##Interstitial lung disease
{multigroup <- colnames(RA_comorbidities_haplotype_data)[c(1:15)]
  rocRes <- lapply(multigroup, function(x) quiteROCFunc(RA_comorbidities_haplotype_data, "lung (yes 1, no 0)", x)$result)
  rocResDat <- do.call(rbind, rocRes)
  rocResDat=as.data.frame(rocResDat)
  rocResDat$variable=multigroup 
  rocResDat<- apply(rocResDat,2,as.character)
  write.csv(rocResDat,file ="Figure3_raw_data14_lung_haplotype_ROC_ALL.csv")}

## hypertension
{rocRes <- lapply(multigroup, function(x) quiteROCFunc(RA_comorbidities_haplotype_data, "hyper (yes 1, no 0)", x)$result)
  rocResDat <- do.call(rbind, rocRes)
  rocResDat
  rocResDat=as.data.frame(rocResDat)
  rocResDat$variable=multigroup 
  rocResDat<- apply(rocResDat,2,as.character)
  write.csv(rocResDat,file ="Figure3_raw_data15_high_blood_pressure_haplotype_ROC_ALL.csv")}

##Osteoporosis
{rocRes <- lapply(multigroup, function(x) quiteROCFunc(RA_comorbidities_haplotype_data, "bone (yes 1, no 0)", x)$result)
  rocResDat <- do.call(rbind, rocRes)
  rocResDat=as.data.frame(rocResDat)
  rocResDat$variable=multigroup 
  rocResDat<- apply(rocResDat,2,as.character)
  write.csv(rocResDat,file ="Figure3_raw_data16_osteoporosis_haplotype_ROC_ALL.csv")}

###After manually combining all the above results, two files are generated, 
###including"AUC for CpGs diagnosis of different RA.csv"and "AUC for haplotype diagnosis of different RA.csv"

##Plot Figures 3D-3G separately
library(readr)
library(RColorBrewer)
library(ggplot2)

AUC_for_CpGs_diagnosis_of_different_RA <- read_csv("AUC for CpGs diagnosis of different RA.csv")
data1=AUC_for_CpGs_diagnosis_of_different_RA

P100=ggplot(data = data1, mapping = aes(x = subgroup, y = `auc(95%CI)...8`,group=variable)) + 
  geom_line(aes(colour=variable))+xlab("group")+ylab("AUC")+theme_classic()

AUC_for_haplotype_diagnosis_of_different_RA <- read_csv("AUC for haplotype diagnosis of different RA.csv")
View(AUC_for_haplotype_diagnosis_of_different_RA)
mycolors=c("black", "#A1DAB4", "#41B6C4", "#2C7FB8", "#253494",
           "blue" ,"#FECC5C" ,"#FD8D3C", "#F03B20" ,"#BD0026",
           "pink", "#FDB863", "grey" ,"#B2ABD2" ,"#5E3C99")

data2=AUC_for_haplotype_diagnosis_of_different_RA
colnames(data2)
P103=ggplot(data = data2, mapping = aes(x = subgroup, y = `auc(95%CI)...8`,group=variable)) + 
  geom_line(aes(colour= variable,linetype=variable))+xlab("group")+ylab("AUC")+theme_classic()+
  scale_color_manual(values= mycolors)+scale_linetype_manual(values=c("twodash", "dotted","solid","longdash","dotdash",
                                                                      "dashed","solid","twodash","twodash","twodash",
                                                                      "solid","solid","solid","solid","solid"
  ))

CpGs <- read_csv("AUC for CpGs diagnosis of RA-related comorbidities.csv")
data3=CpGs
P101=ggplot(data = data3, mapping = aes(x = subgroup, y = `auc(95%CI)...8`,group=variable)) + 
  geom_line(aes(colour=variable))+xlab("group")+ylab("AUC")+theme_classic()

haplotype <- read_csv("AUC for haplotype diagnosis of RA complications.csv")
data4=haplotype
P102=ggplot(data = data4, mapping = aes(x = subgroup, y = `auc(95%CI)...8`,group=variable)) + 
  geom_line(aes(colour= variable,linetype=variable))+xlab("group")+ylab("AUC")+theme_classic()+
  scale_color_manual(values= mycolors)+scale_linetype_manual(values=c("twodash", "dotted","solid","longdash","dotdash",
                                                                      "dashed","solid","twodash","twodash","twodash",
                                                                      "solid","solid","solid","solid","solid"
  ))

##Plot boxplots of all haplotypes
##Process raw data

Haplotype_data_RA1=Haplotype_data_RA[,c(1:15)]
Haplotype_data_RA1$ID=rownames(Haplotype_data_RA1)
library(reshape2)
Haplotype_data_RA2=melt(Haplotype_data_RA1)
Haplotype_data_RA2$group=ifelse(str_detect(Haplotype_data_RA2$ID,"OA"),"OA",
                   ifelse(str_detect(Haplotype_data_RA2$ID,"HC"),"HC",
                          " RA"))
Haplotype_data_RA2$group=as.factor(Haplotype_data_RA2$group)

data_CCCCCCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="CCCCCCC",]
data_TCCCCCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="TCCCCCC",]
data_TTTTTTT=Haplotype_data_RA2[Haplotype_data_RA2$variable=="TTTTTTT",]
data_TTCCCCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="TTCCCCC",]
data_CCCTCCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="CCCTCCC",]
data_TTTTTCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="TTTTTCC",]
data_CTCCCCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="CTCCCCC",]
data_CCTCCCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="CCTCCCC",]
data_TTTTCCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="TTTTCCC",]
data_TCCTCCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="TCCTCCC",]
data_TCTCCCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="TCTCCCC",]
data_TTTTTTC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="TTTTTTC",]
data_CCCCCCT=Haplotype_data_RA2[Haplotype_data_RA2$variable=="CCCCCCT",]
data_CCCCCTC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="CCCCCTC",]
data_TTTCTCC=Haplotype_data_RA2[Haplotype_data_RA2$variable=="TTTCTCC",]

##remove outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
data_TCCTCCC$value=remove_outliers(data_TCCTCCC$value)
data_TCTCCCC$value=remove_outliers(data_TCTCCCC$value)
data_TTTCTCC$value=remove_outliers(data_TTTCTCC$value)
data_TCCCCCC$value=remove_outliers(data_TCCCCCC$value)
data_CCTCCCC$value=remove_outliers(data_CCTCCCC$value)

library(ggstatsplot) 
P1=ggbetweenstats(data_CCCCCCC,x = "group",y = "value",color = "group",xlab = "CCCCCCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P2=ggbetweenstats(data_TCCCCCC,x = "group",y = "value",color = "group",xlab = "TCCCCCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P3=ggbetweenstats(data_TTTTTTT,x = "group",y = "value",color = "group",xlab = "TTTTTTT",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P4=ggbetweenstats(data_TTCCCCC,x = "group",y = "value",color = "group",xlab = "TTCCCCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P5=ggbetweenstats(data_CCCTCCC,x = "group",y = "value",color = "group",xlab = "CCCTCCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P6=ggbetweenstats(data_TTTTTCC,x = "group",y = "value",color = "group",xlab = "TTTTTCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P7=ggbetweenstats(data_CTCCCCC,x = "group",y = "value",color = "group",xlab = "CTCCCCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P8=ggbetweenstats(data_CCTCCCC,x = "group",y = "value",color = "group",xlab = "CCTCCCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P9=ggbetweenstats(data_TTTTCCC,x = "group",y = "value",color = "group",xlab = "TTTTCCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE,outlier.shape = NA)
P10=ggbetweenstats(data_TCCTCCC,x = "group",y = "value",color = "group",xlab = "TCCTCCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE,outlier.shape = NA)
P11=ggbetweenstats(data_TCTCCCC,x = "group",y = "value",color = "group",xlab = "TCTCCCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE,outlier.shape = NA)
P12=ggbetweenstats(data_TTTTTTC,x = "group",y = "value",color = "group",xlab = "TTTTTTC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P13=ggbetweenstats(data_CCCCCCT,x = "group",y = "value",color = "group",xlab = "CCCCCCT",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P14=ggbetweenstats(data_CCCCCTC,x = "group",y = "value",color = "group",xlab = "CCCCCTC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)
P15=ggbetweenstats(data_TTTCTCC,x = "group",y = "value",color = "group",xlab = "TTTCTCC",ylab = "methylation ratio",K=3,grouping.var="group",pairwise.comparisons = TRUE)



library(cowplot)
library(patchwork)
layout <- 
  'ABB
CDD
EFF
GGH
'

P1+P100+P3+P101+P6+P103+P102+plot_layout(design = layout)

ggsave(filename = "Figure3.pdf",width = 18.5,height = 18)

##Figure_S2

P2+P4+P5+P7+P8+P9+P10+P11+P12+P13+P14+P15
ggsave(P2+P4+P5+P7+P8+P9+P10+P11+P12+P13+P14+P15,
       filename = "FigureS2.pdf",width = 22,height = 10)

save.image("Figure3_FigureS2_Raw.Rdata")
