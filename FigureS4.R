###Load R package
library(pROC)
library(glmnet)
library(ggplot2)
library(tidyverse)

##Process raw data
load("FigureS4_raw_data.Rdata")

data_raw_Roc_RA$CCP=as.numeric(data_raw_Roc_RA$CCP)
colnames(data_raw_Roc_RA)
data_raw_Roc_RA=data_raw_Roc_RA[,c(1:16,29:30)]
colnames(data_raw_Roc_RA)=toupper(colnames(data_raw_Roc_RA))
colnames(data_raw_Roc_RA)[16]="group"

##ccccccc
fit1 <- glm(group ~ CCCCCCC + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  
data_raw_Roc_RA$prob <- predict(fit1, 
                                newdata=data_raw_Roc_RA, 
                                type="response")

roc1 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$CCCCCCC) 
roc2 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$CCP)
roc3 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g1 <- ggroc(list(`CCCCCCC(0.7219)`=roc1, `CCP(0.9498)`=roc2, `combanation(0.9693)`=roc3
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


##tcccccc
fit2 <- glm(group ~ TCCCCCC + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  
data_raw_Roc_RA$prob <- predict(fit2, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc4 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$TCCCCCC) 
roc5 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g2 <- ggroc(list(`TCCCCCC(0.5579)`=roc4, `CCP(0.9498)`=roc2, `combanation(0.9572)`=roc5
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

#ttttttt
fit3 <- glm(group ~ TTTTTTT + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit3, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc6 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$TTTTTTT) 
roc7 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g3 <- ggroc(list(`TTTTTTT(0.6744)`=roc6, `CCP(0.9498)`=roc2, `combanation(0.9561)`=roc7
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


#ccctccc
fit4 <- glm(group ~ CCCTCCC + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit4, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc8 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$CCCTCCC) 
roc9 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g4 <- ggroc(list(`CCCTCCC(0.6013)`=roc8, `CCP(0.9498)`=roc2, `combanation(0.9556)`=roc9
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



#tttttcc
fit5 <- glm(group ~ TTTTTCC + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit5, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc10 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$TTTTTCC) 
roc11 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g5 <- ggroc(list(`TTTTTCC(0.6754)`=roc10, `CCP(0.9498)`=roc2, `combanation(0.9624)`=roc11
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


#ctccccc
fit6 <- glm(group ~ CTCCCCC + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit6, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc12 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$CTCCCCC) 
roc13 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g6 <- ggroc(list(`CTCCCCC(0.5901)`=roc12, `CCP(0.9498)`=roc2, `combanation(0.9627)`=roc13
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


#cctcccc

fit7 <- glm(group ~ CCTCCCC + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit7, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc14 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$CCTCCCC) 
roc15 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g7 <- ggroc(list(`CCTCCCC(0.56)`=roc14, `CCP(0.9498)`=roc2, `combanation(0.9549)`=roc15
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



#ttttccc
fit8 <- glm(group ~ TTTTCCC + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit8, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc16<- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$TTTTCCC) 
roc17 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g8 <- ggroc(list(`TTTTCCC(0.6409)`=roc16, `CCP(0.9498)`=roc2, `combanation(0.956)`=roc17
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



#tcctccc
fit9 <- glm(group ~ TCCTCCC + CCP,
            data=data_raw_Roc_RA,
            family = binomial())  

data_raw_Roc_RA$prob <- predict(fit9, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc18<- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$TCCTCCC) 
roc19 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g9 <- ggroc(list(`TCCTCCC(0.4846)`=roc18, `CCP(0.9498)`=roc2, `combanation(0.946)`=roc19
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

#tctcccc
fit10 <- glm(group ~ TCTCCCC + CCP,
             data=data_raw_Roc_RA,
             family = binomial())  

data_raw_Roc_RA$prob <- predict(fit10, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc20<- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$TCTCCCC) 
roc21 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g10 <- ggroc(list(`TCTCCCC(0.5419)`=roc20, `CCP(0.9498)`=roc2, `combanation(0.9511)`=roc21
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



#ttttttc
fit11 <- glm(group ~ TTTTTTC + CCP,
             data=data_raw_Roc_RA,
             family = binomial())  

data_raw_Roc_RA$prob <- predict(fit11, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc22<- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$ TTTTTTC) 
roc23 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g11 <- ggroc(list(`TTTTTTC(0.6079)`=roc22, `CCP(0.9498)`=roc2, `combanation(0.9583)`=roc23
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


#cccccct
fit12 <- glm(group ~ CCCCCCT + CCP,
             data=data_raw_Roc_RA,
             family = binomial())  

data_raw_Roc_RA$prob <- predict(fit12, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc24<- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$CCCCCCT) 
roc25 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g12 <- ggroc(list(`CCCCCCT(0.5873)`=roc24, `CCP(0.9498)`=roc2, `combanation(0.9562)`=roc25
))

P12=g12+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
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
P12

#ccccctc
fit13 <- glm(group ~ CCCCCTC + CCP,
             data=data_raw_Roc_RA,
             family = binomial())  

data_raw_Roc_RA$prob <- predict(fit13, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc26<- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$CCCCCTC) 
roc27 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g13 <- ggroc(list(`CCCCCTC(0.5791)`=roc26, `CCP(0.9498)`=roc2, `combanation(0.9585)`=roc27
))

P13=g13+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
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
P13

#tttctcc
fit14 <- glm(group ~ TTTCTCC + CCP,
             data=data_raw_Roc_RA,
             family = binomial())  

data_raw_Roc_RA$prob <- predict(fit14, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc28<- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$TTTCTCC) 
roc29 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)

g14 <- ggroc(list(`TTTCTCC(0.5821)`=roc28, `CCP(0.9498)`=roc2, `combanation(0.9673)`=roc29
))

P14=g14+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
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
P14

#TTCCCCC
fit15 <- glm(group ~ TTCCCCC + CCP,
             data=data_raw_Roc_RA,
             family = binomial())  

data_raw_Roc_RA$prob <- predict(fit15, 
                                newdata=data_raw_Roc_RA, 
                                type="response")
roc30<- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$TTCCCCC) 
roc31 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g15 <- ggroc(list(`TTCCCCC(0.6565)`=roc30, `CCP(0.9498)`=roc2, `combanation(0.9699)`=roc31
))

P15=g15+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
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
P15

#CCCCCCCC+RF
fit16 <- glm(group ~ CCCCCCC + RF,
             data=data_raw_Roc_RA,
             family = binomial())  
data_raw_Roc_RA$prob <- predict(fit1, 
                                newdata=data_raw_Roc_RA, 
                                type="response")

roc32 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$CCCCCCC) 
roc33 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$RF)
roc34 <- roc(data_raw_Roc_RA$group, data_raw_Roc_RA$prob)
g16 <- ggroc(list(`CCCCCCC(0.9213)`=roc32, `RF(0.9213)`=roc33, `combanation(0.9693)`=roc34
))

P16=g16+annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
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
P16



library(patchwork)
P1+P2+P3+P4+P5+P6+P7+P8+P9+P10+P11+P12+P13+P14+P15+P16
ggsave(filename = "FigureS4.pdf",width =13,height = 11)
save.image(file = "FigureS4.Rdata")  

