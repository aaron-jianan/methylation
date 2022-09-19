###Load the R package
library(pROC)
library(glmnet)
library(ggplot2)
library(tidyverse)

##load raw data
load("Figure2_Raw.Rdata")
class(data_raw_Roc_RA$RF)
data_raw_Roc_RA$CCP=as.numeric(data_raw_Roc_RA$CCP)

fit1 <- glm(group ~ `46898116` + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit1, 
                                newdata=data_raw_Roc_RA, 
                                type="response")

roc1 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$`46898116`) 
roc2 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$CCP)
roc3 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g1 <- ggroc(list(`46898116(0.7345)`=roc1, `CCP(0.9498)`=roc2, `combanation(0.9709)`=roc3
))

P1=g1+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+
 theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(legend.text = element_text( size = 8, face = "bold"))
P1


fit2 <- glm(group ~ `46898066` + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  
data_raw_Roc_RA$prob <- predict(fit2, 
                                newdata=data_raw_Roc_RA, 
                                type="response")

roc4 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$`46898066`) 
roc5 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g2 <- ggroc(list(`46898066(0.6934)`=roc4, `CCP(0.9498)`=roc2, `combanation(0.965)`=roc5
))

P2=g2+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(legend.text = element_text( size = 8, face = "bold"))
P2


fit3 <- glm(group ~ `46898048` + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit3, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc6 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$`46898048`) 
roc7 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g3 <- ggroc(list(`46898048(0.6921)`=roc6, `CCP(0.9498)`=roc2, `combanation(0.9641)`=roc7
))

P3=g3+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(legend.text = element_text( size = 8, face = "bold"))
P3


fit4 <- glm(group ~ `46898042` + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  
data_raw_Roc_RA$prob <- predict(fit4, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc8 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$`46898042`) 
roc9 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g4 <- ggroc(list(`46898042(0.6832)`=roc8, `CCP(0.9498)`=roc2, `combanation(0.9605)`=roc9
))

P4=g4+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(legend.text = element_text( size = 8, face = "bold"))
P4


fit5 <- glm(group ~ `46898024` + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit5, 
                                newdata=data_raw_Roc_RA, 
                                type="response")

roc10 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$`46898024`) 
roc11 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g5 <- ggroc(list(`46898024(0.6993)`=roc10, `CCP(0.9498)`=roc2, `combanation(0.965)`=roc11
))

P5=g5+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(legend.text = element_text( size = 8, face = "bold"))
P5


fit6 <- glm(group ~ `46898006` + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit6, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc12 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$`46898006`) 
roc13 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g6 <- ggroc(list(`46898006(0.6767)`=roc12, `CCP(0.9498)`=roc2, `combanation(0.9585)`=roc13
))

P6=g6+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(legend.text = element_text( size = 8, face = "bold"))
P6


fit7 <- glm(group ~ `46898004` + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  
data_raw_Roc_RA$prob <- predict(fit7, 
                                newdata=data_raw_Roc_RA, 
                                type="response")

roc14 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$`46898004`) 
roc15 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g7 <- ggroc(list(`46898004(0.6698)`=roc14, `CCP(0.9498)`=roc2, `combanation(0.9556)`=roc15
))

P7=g7+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(legend.text = element_text( size = 8, face = "bold"))
P7

fit8 <- glm(group ~ cg15692052 + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  
data_raw_Roc_RA$prob <- predict(fit8, 
                                newdata=data_raw_Roc_RA, 
                                type="response")

roc16 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$cg15692052) 
roc17 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g8 <- ggroc(list(`cg15692052(0.7057)`=roc16, `CCP(0.9498)`=roc2, `combanation(0.9648)`=roc17
))

P8=g8+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(legend.text = element_text( size = 8, face = "bold"))
P8


fit9 <- glm(group ~ `46898116` + RF,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit9, 
                                newdata=data_raw_Roc_RA, 
                                type="response")

roc18 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$`46898116`) 
roc19<- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$RF)
roc20<- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g9 <- ggroc(list(`46898116(0.7345)`=roc18, `RF(0.9213)`=roc19, `combanation(0.9256)`=roc20
))

P9=g9+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(legend.text = element_text( size = 8, face = "bold"))
P9


fit10 <- glm(group ~`46898066`  + RF,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit10, 
                                newdata=data_raw_Roc_RA, 
                                type="response")

roc21 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$`46898066`) 
roc22 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g10 <- ggroc(list(`46898066(0.6934)`=roc21, `RF(0.9213)`=roc19, `combanation(0.9215)`=roc22
))

P10=g10+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(legend.text = element_text( size = 8, face = "bold"))
P10

fit11 <- glm(group ~ cg15692052 + RF,
            data=data_raw_Roc_RA,
            family = binomial())  
data_raw_Roc_RA$prob <- predict(fit11, 
                                newdata=data_raw_Roc_RA, 
                                type="response")


roc23 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$cg15692052) 
roc24 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g11<- ggroc(list(`cg15692052_average(0.7057)`=roc23, `RF(0.9213)`=roc19, `combanation(0.9209)`=roc24
))

P11=g11+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,
                                 colour = "black"),
        axis.title = element_text(colour = "black",
                                  size = 10,
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 10,
                                 face = "bold"),
        text = element_text(size = 8,
                            color = "black",
                            family = "serif"),
        legend.background = element_rect(fill = "transparent"))+ 
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
theme(legend.text = element_text( size = 8, face = "bold"))
  
P11


library(patchwork)
P1+P2+P3+P4+P5+P6+P7+P8+P9+P10+P11
ggsave(filename = "Figure2.pdf",width =14,height = 9)
save.image(file = "Figure2_Raw.Rdata")





