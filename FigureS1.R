##FigureS1-A

library(readr)
data_raw_Roc_RA <- read_csv("data_raw_Roc_RA.csv")
names_id=data_raw_Roc_RA$...1
data_raw_Roc_RA=data_raw_Roc_RA[,-1]
row.names(data_raw_Roc_RA)=names_id

##Correlation analysis
library(Hmisc)
cor_matr_Hmisc=rcorr(as.matrix(data_raw_Roc_RA))

##Correlation matrix
cor_matr_Hmisc$r
cor_matr_Hmisc$n
cor_matr_Hmisc$P

library(PerformanceAnalytics)
pdf(file="FigureS1A.pdf")
chart.Correlation(data_raw_Roc_RA,histogram = F,pch=19)
dev.off()

##FigureS1-B
library(readr)
Haplotype_raw_data <- read_csv("Haplotype_raw_data.csv")
Haplotype_raw_data=Haplotype_raw_data[,-1]
row.names(Haplotype_raw_data)=names_id

##Correlation analysis
library(Hmisc)
cor_matr_Hmisc=rcorr(as.matrix(Haplotype_raw_data))

##Correlation matrix
cor_matr_Hmisc$r
cor_matr_Hmisc$n
cor_matr_Hmisc$P

library(PerformanceAnalytics)
pdf(file="FigureS1B.pdf")
chart.Correlation(Haplotype_raw_data,histogram = F,pch=19)
dev.off()

save.image(file = "FigureS1.Rdata")



