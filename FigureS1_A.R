load("FigureS1_A.Rdata")

##Correlation analysis
colnames(data_raw)=c("46898116",'46898066','46898048',"46898042",
                     "46898024",	"46898006",	"46898004")

data_raw$cg15692052=combine_data_raw_merge_data$cg15692052
rm(mean)

library(Hmisc)
cor_matr_Hmisc=rcorr(as.matrix(data_raw))

##Correlation matrix
cor_matr_Hmisc$r
cor_matr_Hmisc$n
cor_matr_Hmisc$P

library(PerformanceAnalytics)
pdf(file="cor_matr_PerformanceAnalytics.pdf")
chart.Correlation(data_raw,histogram = F,pch=19)
dev.off()

save.image(file = "FigureS1_A.Rdata")



