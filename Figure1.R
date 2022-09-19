library(ggstatsplot)
library(patchwork)
library(ggplot2)
load("figure1_Raw.Rdata") 

P8=ggbetweenstats(data_raw_cg15692052_35_185,x = "group",y = "Percentage",color = "group",xlab ="46898006",ylab = "methylation percentage",K=3
                  ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")

P2=ggbetweenstats(data_raw_cg15692052_35_75,x = "group",y = "Percentage",color = "group",xlab ="46898116",ylab = "methylation percentage",K=3
                  ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")

P3=ggbetweenstats(data_raw_cg15692052_35_125,x = "group",y = "Percentage",color = "group",xlab ="46898066",ylab = "methylation percentage",K=3
                  ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")

P4=ggbetweenstats(data_raw_cg15692052_35_143,x = "group",y = "Percentage",color = "group",xlab ="46898048",ylab = "methylation percentage",K=3
                  ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")

P5=ggbetweenstats(data_raw_cg15692052_35_149,x = "group",y = "Percentage",color = "group",xlab ="46898042",ylab = "methylation percentage",K=3
                  ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")

P6=ggbetweenstats(data_raw_cg15692052_35_167,x = "group",y = "Percentage",color = "group",xlab ="46898024",ylab = "methylation percentage",K=3
                  ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")

P7=ggbetweenstats(data_raw_cg15692052_35_187,x = "group",y = "Percentage",color = "group",xlab ="46898004",ylab = "methylation percentage",K=3
                  ,grouping.var="group",pairwise.comparisons = TRUE,pairwise.display="significant")

P1=ggbetweenstats(data_m1,x = "variable",y = "value",color = "group",xlab ="cg15692052",ylab = "methylation level",K=3
                  ,nboot = 10,
                  messages = FALSE,
                  effsize.type = "unbiased", # type of effect size (unbiased = omega)
                  partial = FALSE, # partial omega or omega?
                  pairwise.comparisons = TRUE, # display results from pairwise comparisons
                  pairwise.display = "significant", # display only significant pairwise comparisons
                  pairwise.annotation = "p.value", # annotate the pairwise comparisons using p-values
                  # adjust p-values for multiple tests using this method
)
P_all=P1+P2+P3+P4+P5+P6+P7+P8
ggsave("cg15692052_all.boxplot.pdf",width = 14,
       height = 14.5)

