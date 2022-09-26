##Load R packages, read and process raw data

library(ggstatsplot)
library(patchwork)
library(ggplot2)
library(readxl)
library(tidyverse)
data1 <- read_excel("Methylation_raw_data.xlsx", 
                                   sheet = "Site")
##draw
##cg15692052_average
{data2=read_excel("Methylation_raw_data.xlsx", 
                  sheet = "gene")
  group=ifelse(str_detect(colnames(data2),"HC"),"HC",
               ifelse(str_detect(colnames(data2),"OA"),"OA","RA"))
  group=as.factor(group)
  
  library(reshape2)
  data_m_cg15692052 <- melt(data2)
  data_m_cg15692052$group=data_m_cg15692052$variable
  data_m_cg15692052$group=ifelse(str_detect(data_m_cg15692052$group,"HC"),"HC",
                                 ifelse(str_detect(data_m_cg15692052$group,"OA"),"OA","RA"))
  
  data_m_cg15692052$group=as.factor(data_m_cg15692052$group)
  P1=ggbetweenstats(data_m_cg15692052,x = "group",y = "value",color = "group",xlab ="cg15692052",ylab = "methylation level",K=3,
                    grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")}

##46898006
{data_4689006=data1[6,]
data_4689006=data_4689006[,-c(1:6)]
group=ifelse(str_detect(colnames(data_4689006),"HC"),"HC",
             ifelse(str_detect(colnames(data_4689006),"OA"),"OA","RA"))
group=as.factor(group)

library(reshape2)
data_m_4689006 <- melt(data_4689006)
data_m_4689006$group=data_m_4689006$variable
data_m_4689006$group=ifelse(str_detect(data_m_4689006$group,"HC"),"HC",
                    ifelse(str_detect(data_m_4689006$group,"OA"),"OA","RA"))
P8=ggbetweenstats(data_m_4689006,x = "group",y = "value",color = "group",xlab ="46898006",ylab = "methylation level",K=3
                  ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")}

##46898116
{data_46898116=data1[1,]
  data_46898116=data_46898116[,-c(1:6)]
  group=ifelse(str_detect(colnames(data_46898116),"HC"),"HC",
               ifelse(str_detect(colnames(data_46898116),"OA"),"OA","RA"))
  group=as.factor(group)
  
  library(reshape2)
  data_m_46898116 <- melt(data_46898116)
  data_m_46898116$group=data_m_46898116$variable
  data_m_46898116$group=ifelse(str_detect(data_m_46898116$group,"HC"),"HC",
                               ifelse(str_detect(data_m_46898116$group,"OA"),"OA","RA"))
  P2=ggbetweenstats(data_m_46898116,x = "group",y = "value",color = "group",xlab ="46898116",ylab = "methylation level",K=3
                    ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")}

##46898066
{data_46898066=data1[2,]
  data_46898066=data_46898066[,-c(1:6)]
  group=ifelse(str_detect(colnames(data_46898066),"HC"),"HC",
               ifelse(str_detect(colnames(data_46898066),"OA"),"OA","RA"))
  group=as.factor(group)
  
  library(reshape2)
  data_m_46898066 <- melt(data_46898066)
  data_m_46898066$group=data_m_46898066$variable
  data_m_46898066$group=ifelse(str_detect(data_m_46898066$group,"HC"),"HC",
                               ifelse(str_detect(data_m_46898066$group,"OA"),"OA","RA"))
  P3=ggbetweenstats(data_m_46898066,x = "group",y = "value",color = "group",xlab ="46898066",ylab = "methylation level",K=3
                    ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")}

##46898048
{data_46898048=data1[3,]
  data_46898048=data_46898048[,-c(1:6)]
  group=ifelse(str_detect(colnames(data_46898048),"HC"),"HC",
               ifelse(str_detect(colnames(data_46898048),"OA"),"OA","RA"))
  group=as.factor(group)
  
  library(reshape2)
  data_m_46898048 <- melt(data_46898048)
  data_m_46898048$group=data_m_46898048$variable
  data_m_46898048$group=ifelse(str_detect(data_m_46898048$group,"HC"),"HC",
                               ifelse(str_detect(data_m_46898048$group,"OA"),"OA","RA"))
  P4=ggbetweenstats(data_m_46898048,x = "group",y = "value",color = "group",xlab ="46898048",ylab = "methylation level",K=3
                    ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")}

##46898042
{data_46898042=data1[4,]
  data_46898042=data_46898042[,-c(1:6)]
  group=ifelse(str_detect(colnames(data_46898042),"HC"),"HC",
               ifelse(str_detect(colnames(data_46898042),"OA"),"OA","RA"))
  group=as.factor(group)
  
  library(reshape2)
  data_m_46898042 <- melt(data_46898042)
  data_m_46898042$group=data_m_46898042$variable
  data_m_46898042$group=ifelse(str_detect(data_m_46898042$group,"HC"),"HC",
                               ifelse(str_detect(data_m_46898042$group,"OA"),"OA","RA"))
  P5=ggbetweenstats(data_m_46898042,x = "group",y = "value",color = "group",xlab ="46898042",ylab = "methylation level",K=3
                    ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")}

##46898024
{data_46898024=data1[5,]
  data_46898024=data_46898024[,-c(1:6)]
  group=ifelse(str_detect(colnames(data_46898024),"HC"),"HC",
               ifelse(str_detect(colnames(data_46898024),"OA"),"OA","RA"))
  group=as.factor(group)
  
  library(reshape2)
  data_m_46898024 <- melt(data_46898024)
  data_m_46898024$group=data_m_46898024$variable
  data_m_46898024$group=ifelse(str_detect(data_m_46898024$group,"HC"),"HC",
                               ifelse(str_detect(data_m_46898024$group,"OA"),"OA","RA"))
  P6=ggbetweenstats(data_m_46898024,x = "group",y = "value",color = "group",xlab ="46898024",ylab = "methylation level",K=3
                    ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")}

##46898004
{data_46898004=data1[7,]
  data_46898004=data_46898004[,-c(1:6)]
  group=ifelse(str_detect(colnames(data_46898004),"HC"),"HC",
               ifelse(str_detect(colnames(data_46898004),"OA"),"OA","RA"))
  group=as.factor(group)
  
  library(reshape2)
  data_m_46898004 <- melt(data_46898004)
  data_m_46898004$group=data_m_46898004$variable
  data_m_46898004$group=ifelse(str_detect(data_m_46898004$group,"HC"),"HC",
                               ifelse(str_detect(data_m_46898004$group,"OA"),"OA","RA"))
  P7=ggbetweenstats(data_m_46898004,x = "group",y = "value",color = "group",xlab ="46898004",ylab = "methylation level",K=3
                    ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")}

##puzzle
P_all=P1+P2+P3+P4+P5+P6+P7+P8
ggsave("Figure1.pdf",width = 14,
       height = 14.5)
save.image(file = "figure1_Raw_update.Rdata")

