library(readxl)
Methylation <- read_excel("Methylation_raw_data.xlsx", 
                          sheet = "Haplotype", col_types = c("text", 
                                                             "text", "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric"))
View(Methylation)

##Process raw data
data=Methylation
data1=data[data$Target=="cg15692052_35",]
data1=data1[,-c(1,3)]
data1$Haplotype=toupper(data1$Haplotype)
rownames(data1)=data1$Haplotype

library(reshape2)
library(tidyverse)
data2=melt(data1)
data2$group=ifelse(str_detect(data2$variable,"OA"),"OA",
                   ifelse(str_detect(data2$variable,"HC"),"HC",
                          " RA"))

data2$group=as.factor(data2$group)

rownames(data1)
data_CCCCCCC=data2[data2$Haplotype=="CCCCCCC",]
data_TCCCCCC=data2[data2$Haplotype=="TCCCCCC",]
data_TTTTTTT=data2[data2$Haplotype=="TTTTTTT",]
data_TTCCCCC=data2[data2$Haplotype=="TTCCCCC",]
data_CCCTCCC=data2[data2$Haplotype=="CCCTCCC",]
data_TTTTTCC=data2[data2$Haplotype=="TTTTTCC",]
data_CTCCCCC=data2[data2$Haplotype=="CTCCCCC",]
data_CCTCCCC=data2[data2$Haplotype=="CCTCCCC",]
data_TTTTCCC=data2[data2$Haplotype=="TTTTCCC",]
data_TCCTCCC=data2[data2$Haplotype=="TCCTCCC",]
data_TCTCCCC=data2[data2$Haplotype=="TCTCCCC",]
data_TTTTTTC=data2[data2$Haplotype=="TTTTTTC",]
data_CCCCCCT=data2[data2$Haplotype=="CCCCCCT",]
data_CCCCCTC=data2[data2$Haplotype=="CCCCCTC",]
data_TTTCTCC=data2[data2$Haplotype=="TTTCTCC",]

## remove outliers
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
library(ggplot2)
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

library(patchwork)
P2+P4+P5+P7+P8+P9+P10+P11+P12+P13+P14+P15
ggsave(P2+P4+P5+P7+P8+P9+P10+P11+P12+P13+P14+P15,
       filename = "FigureS2.pdf",width = 22,height = 10)
save.image(file = "FigureS2.Rdata")

