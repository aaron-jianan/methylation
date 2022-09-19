

library(readr)
AUC_for_CpGs_diagnosis_of_different_RA <- read_csv("AUC for CpGs diagnosis of different RA.csv")
data1=AUC_for_CpGs_diagnosis_of_different_RA

colnames(data1)
library(RColorBrewer)


library(ggplot2)
P=ggplot(data = data1, mapping = aes(x = subgroup, y = `auc(95%CI)...8`,group=variable)) + 
  geom_line(aes(colour=variable))+xlab("group")+ylab("AUC")+theme_classic()

library(readr)
AUC_for_haplotype_diagnosis_of_different_RA <- read_csv("AUC for CpGs diagnosis of different RA.csv")
View(AUC_for_haplotype_diagnosis_of_different_RA)
mycolors=c("black", "#A1DAB4", "#41B6C4", "#2C7FB8", "#253494",
           "blue" ,"#FECC5C" ,"#FD8D3C", "#F03B20" ,"#BD0026",
           "pink", "#FDB863", "grey" ,"#B2ABD2" ,"#5E3C99")

data2=AUC_for_haplotype_diagnosis_of_different_RA
colnames(data2)
library(RColorBrewer)

library(ggplot2)
P2=ggplot(data = data2, mapping = aes(x = subgroup, y = `auc(95%CI)...8`,group=variable)) + 
  geom_line(aes(colour= variable,linetype=variable))+xlab("group")+ylab("AUC")+theme_classic()+
  scale_color_manual(values= mycolors)+scale_linetype_manual(values=c("twodash", "dotted","solid","longdash","dotdash",
                                                                      "dashed","solid","twodash","twodash","twodash",
                                                                      "solid","solid","solid","solid","solid"
  ))

library(readr)
CpGs <- read_csv("AUC for CpGs diagnosis of RA-related comorbidities.csv")
View(CpGs)

data3=CpGs
library(RColorBrewer)

library(ggplot2)
P1=ggplot(data = data3, mapping = aes(x = subgroup, y = `auc(95%CI)...8`,group=variable)) + 
  geom_line(aes(colour=variable))+xlab("group")+ylab("AUC")+theme_classic()


library(readr)
haplotype <- read_csv("AUC for haplotype diagnosis of RA complications.csv")
View(haplotype)
data4=haplotype
library(RColorBrewer)
P3=ggplot(data = data4, mapping = aes(x = subgroup, y = `auc(95%CI)...8`,group=variable)) + 
  geom_line(aes(colour= variable,linetype=variable))+xlab("group")+ylab("AUC")+theme_classic()+
  scale_color_manual(values= mycolors)+scale_linetype_manual(values=c("twodash", "dotted","solid","longdash","dotdash",
                                                                      "dashed","solid","twodash","twodash","twodash",
                                                                      "solid","solid","solid","solid","solid"
  ))

library(cowplot)
P20=P
P21=P2
P22=P1
P23=P3

save(AUC_for_CpGs_diagnosis_of_different_RA,P6,AUC_for_haplotype_diagnosis_of_different_RA,P20,CpGs,
     P1,haplotype,P3,P20,P21,P22,file ="Figure3_Raw.Rdata")

load("Figure3_Raw.Rdata")
library(patchwork)
layout <- 
  'ABB
CDD
EFF
GGH
'

P1+P20+P3+P22+P6+P21+P23+plot_layout(design = layout)

ggsave(filename = "Figure3.pdf",width = 18.5,height = 18)

