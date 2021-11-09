#--------R script written by JVR for chi-square analysis of biochemical and immunological data---------

library(caret)
library(glmnet)
library(mlbench)
library(psych)
library(arm)
library(MASS)
library(dplyr)
library(magrittr)
library(olsrr)
library(ggplot2)
library(cowplot)
library(ROCR)
library(oddsratio)
library(ResourceSelection)
library(pROC)
#Data reading 
myData=read.csv(file.choose(),header=T, row.names=1)
#Data summary
myData
colnames(myData)
rownames(myData)
str(myData)


#-------------------------------------------------------------------------------------
#MODEL1 - Univariate logistic regression model
#-------------------------------------------------------------------------------------
# set.seed(123) #not important. Will get same results if not used as well
model_1 <- glm(Ulcer.healing.status~ MFI.CD14..CD63 , data = myData, family = 'binomial', maxit=100)
summary(model_1)
# -----------------------------------------
# #ODDS RATIO
# -----------------------------------------
exp(cbind("Odds ratio" = coef(model_1), confint(model_1, level = 0.95)))
# -----------------------------------------
#Goodness of fit - model evaluation statistics
# -----------------------------------------
#Chisquare test
with(model_1, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))
#Hosmer lemeshow test
glm1_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(model_1), g=10)
print(glm1_hl)
# -----------------------------------------
#ROC and AUC
# -----------------------------------------
probabilities1 <- model_1 %>% predict(myData, type = "response")
predict_b1<-prediction(probabilities1, myData$Ulcer.healing.status)
perf_b1<-performance(predict_b1, measure="tpr", x.measure="fpr")
AUC_b1<-performance(predict_b1, measure="auc")
#AUC
AUC_b1@y.values
#Plot ROC Curve
tiff("MFI.CD14..CD11b.tiff", width = 700, height = 600)

par(mar=c(5,5,5,5), cex.axis=1.5)
par(pty="s")
roc1<-roc(myData$Ulcer.healing.status,probabilities1, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
          main="MFI CD14 CD11b", cex.lab=1.7, cex.main=1.4,type="o",
          xlab="False positive rate", ylab="True positive rate", font.lab=2)
dev.off()
# -----------------------------------------


#-------------------------------------------------------------------------------------
#MODEL2 - Bayesian logistic regression model
#-------------------------------------------------------------------------------------
model_bayes <- bayesglm(Ulcer.healing.status~.  , data = myData, family = 'binomial', maxit=100)
summary(model_bayes)

exp(model_bayes$coefficients)
exp(confint(model_bayes))

#Goodness of fit - model evaluation statistics
with(model_bayes, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))

#Hosmer lemeshow test
glm2_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(model_bayes), g=10)
print(glm2_hl)

# -----------------------------------------
#ROC and AUC
# -----------------------------------------
probabilities2 <- model_bayes %>% predict(myData, type = "response")
predict_b2<-prediction(probabilities2, myData$Ulcer.healing.status)
perf_b2<-performance(predict_b2, measure="tpr", x.measure="fpr")
AUC_b2<-performance(predict_b2, measure="auc")
#AUC
AUC_b2@y.values
#Plot ROC Curve
tiff("Model 1 - Bayes Model.tiff", width = 700, height = 600)

par(mar=c(5,5,5,5), cex.axis=1.5)
par(pty="s")
roc1<-roc(myData$Ulcer.healing.status,probabilities2, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
          main="Model 1 - Bayes model", cex.lab=1.7, cex.main=1.4,type="o",
          xlab="False positive rate", ylab="True positive rate", font.lab=2)
dev.off()
# -----------------------------------------


#-------------------------------------------------------------------------------------
#MODEL3 - Full logistic regression model - ALL Clinical parameters
#-------------------------------------------------------------------------------------
model_glm_fullclin <- glm(Ulcer.healing.status~Age                       + BMI                        + HbA1c.. +                   
                 FBS..mg.dl.          +       PPBS..mg.dl.          +             
                    Total.Cholesterol..mg.dl. + HDL.Cholesterol..mg.dl.  +  
                 LDL.Cholesterol..mg.dl.+     S..Triglycerides..mg.dl. + VLDL.Cholesterol..mg.dl.  +
                 Total.Chol.HDL.Chol  +     ESR              +        Alkaline.Phosphatase..IU.L. +
                   X..CD14.             +      X..CD15.High    +         X.CD15...medium,           
                   data = myData, family = 'binomial', maxit=100)
summary(model_glm_fullclin)

#ODDS Ratios
exp(model_glm_fullclin$coefficients)
exp(confint(model_glm_fullclin))

#Goodness of fit - model evaluation statistics
with(model_glm_fullclin, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))

#Hosmer lemeshow test
glm3_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(model_glm_fullclin), g=10)
print(glm3_hl)

# -----------------------------------------
#ROC and AUC
# -----------------------------------------
probabilities3 <- model_glm_fullclin %>% predict(myData, type = "response")
predict_b3<-prediction(probabilities3, myData$Ulcer.healing.status)
perf_b3<-performance(predict_b3, measure="tpr", x.measure="fpr")
AUC_b3<-performance(predict_b3, measure="auc")
#AUC
AUC_b3@y.values
#Plot ROC Curve
tiff("Model 2 - Clinical model 1.tiff", width = 700, height = 600)
par(mar=c(5,5,5,5), cex.axis=1.5)
par(pty="s")
roc1<-roc(myData$Ulcer.healing.status,probabilities3, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
          main="Model 2 - Clinical model 1", cex.lab=1.7, cex.main=1.4,type="o",
          xlab="False positive rate", ylab="True positive rate", font.lab=2)

dev.off()
# -----------------------------------------

#-------------------------------------------------------------------------------------
#MODEL4 - Reduced logistic regression model - Clinical parameters
#-------------------------------------------------------------------------------------
#parameters selected based on p value less than 0.25 from full model

model_glm_redclin <- glm(Ulcer.healing.status~ FBS..mg.dl.  +Total.Cholesterol..mg.dl. + 
                           LDL.Cholesterol..mg.dl.+S..Triglycerides..mg.dl.+VLDL.Cholesterol..mg.dl.+
                           X..CD14. ,           
                         data = myData, family = 'binomial', maxit=100)
summary(model_glm_redclin)

#ODDS Ratios
exp(model_glm_redclin$coefficients)
exp(confint(model_glm_redclin))

#Goodness of fit - model evaluation statistics
with(model_glm_redclin, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))

#Hosmer lemeshow test
glm4_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(model_glm_redclin), g=10)
print(glm4_hl)

# -----------------------------------------
#ROC and AUC
# -----------------------------------------
probabilities4 <- model_glm_redclin %>% predict(myData, type = "response")
predict_b4<-prediction(probabilities4, myData$Ulcer.healing.status)
perf_b4<-performance(predict_b4, measure="tpr", x.measure="fpr")
AUC_b4<-performance(predict_b4, measure="auc")
#AUC
AUC_b4@y.values
#Plot ROC Curve
tiff("Model 3 - Clinical model 2.tiff", width = 700, height = 600)
par(mar=c(5,5,5,5), cex.axis=1.5)
par(pty="s")
roc1<-roc(myData$Ulcer.healing.status,probabilities4, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
          main="Model 3 - Clinical model 2", cex.lab=1.7, cex.main=1.4,type="o",
          xlab="False positive rate", ylab="True positive rate", font.lab=2)
dev.off()
# -----------------------------------------

#-------------------------------------------------------------------------------------
#MODEL5 - Stepwise reduced logistic regression model - Clinical parameters
#-------------------------------------------------------------------------------------

  #REDUCED MODEL
  full.model1 <- glm(Ulcer.healing.status ~ 
                       Age                       + BMI                        + HbA1c.. +                   
                       FBS..mg.dl.          +       PPBS..mg.dl.          +             
                       Total.Cholesterol..mg.dl. + HDL.Cholesterol..mg.dl.  +  
                       LDL.Cholesterol..mg.dl.+     S..Triglycerides..mg.dl. + VLDL.Cholesterol..mg.dl.  +
                       Total.Chol.HDL.Chol  +     ESR              +        Alkaline.Phosphatase..IU.L. +
                       X..CD14.             +      X..CD15.High    +         X.CD15...medium , 
                     data = myData, family = 'binomial', maxit=100)
  
  #STEP MODEL
  s_model <- full.model1 %>% stepAIC(trace = FALSE)
  smod_coef<-coef(s_model)
  print(smod_coef)
  summary(s_model)
  
  #ODDS Ratios
  exp(s_model$coefficients)
  exp(confint(s_model))
  
  #Goodness of fit - model evaluation statistics
  with(s_model, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))
  
  #Hosmer lemeshow test
  glm5_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(s_model), g=10)
  print(glm5_hl)
  
  # -----------------------------------------
  #ROC and AUC
  # -----------------------------------------
  probabilities5 <- s_model %>% predict(myData, type = "response")
  predict_b5<-prediction(probabilities5, myData$Ulcer.healing.status)
  perf_b5<-performance(predict_b5, measure="tpr", x.measure="fpr")
  AUC_b5<-performance(predict_b5, measure="auc")
  #AUC
  AUC_b5@y.values
  #Plot ROC Curve
  tiff("Model 4 - Clinical model 3.tiff", width = 700, height = 600)
  par(mar=c(5,5,5,5), cex.axis=1.5)
  par(pty="s")
  roc1<-roc(myData$Ulcer.healing.status,probabilities5, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
            main="Model 4 - Clinical model 3", cex.lab=1.7, cex.main=1.4,type="o",
            xlab="False positive rate", ylab="True positive rate", font.lab=2)
  
  dev.off()
  # -----------------------------------------

  #-------------------------------------------------------------------------------------
  #MODEL6 - Full logistic regression model - ALL Immune parameters
  #-------------------------------------------------------------------------------------
  model_glm_full_immun <- glm(Ulcer.healing.status~
                                                                MFI.CD15..CD16   +   MFI.CD15..CD62L  +   MFI.CD15..CD63   +         
                                MFI.CD15..CD282     +        MFI.CD15..CD284    +         MFI.CD15..HLA.DR  +
                                MFI.CD15..CD11b     +        MFI.CD14..CD16     +         MFI.CD14..CD62L   +
                                MFI.CD14..CD63      +        MFI.CD14..CD282    +         MFI.CD14..CD284   +
                                MFI.CD14..CD36         +     MFI.CD14..HLA.DR     +       MFI.CD14..CD11b  ,
                              data = myData, family = 'binomial', maxit=100)
  summary(model_glm_full_immun)
  
  #ODDS Ratios
  exp(model_glm_full_immun$coefficients)
  exp(confint(model_glm_full_immun))
  
  #Goodness of fit - model evaluation statistics
  with(model_glm_full_immun, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))
  
  #Hosmer lemeshow test
  glm6_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(model_glm_full_immun), g=10)
  print(glm6_hl)
  
  # -----------------------------------------
  #ROC and AUC
  # -----------------------------------------
  probabilities6 <- model_glm_full_immun %>% predict(myData, type = "response")
  predict_b6<-prediction(probabilities6, myData$Ulcer.healing.status)
  perf_b6<-performance(predict_b6, measure="tpr", x.measure="fpr")
  AUC_b6<-performance(predict_b6, measure="auc")
  #AUC
  AUC_b6@y.values
  #Plot ROC Curve
  tiff("Model 5 - Immune model 1.tiff", width = 700, height = 600)
  par(mar=c(5,5,5,5), cex.axis=1.5)
  par(pty="s")
  roc1<-roc(myData$Ulcer.healing.status,probabilities6, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
            main="Model 5 - Immune model 1", cex.lab=1.7, cex.main=1.4,type="o",
            xlab="False positive rate", ylab="True positive rate", font.lab=2)
  
  dev.off()
  # -----------------------------------------

  #-------------------------------------------------------------------------------------
  #MODEL7 - Full logistic regression model - Selected Immune parameters
  #-------------------------------------------------------------------------------------
  model_glm_red_immun <- glm(Ulcer.healing.status~ MFI.CD15..CD11b     +        MFI.CD14..CD16     +         MFI.CD14..CD62L   +
                                MFI.CD14..CD63      +        MFI.CD14..CD282    +        
                                MFI.CD14..CD36         +     MFI.CD14..CD11b  ,
                            data = myData, family = 'binomial', maxit=100)
  summary(model_glm_red_immun)
  
  #ODDS Ratios
  exp(model_glm_red_immun$coefficients)
  exp(confint(model_glm_red_immun))
  
  #Goodness of fit - model evaluation statistics
  with(model_glm_red_immun, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))
  
  #Hosmer lemeshow test
  glm7_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(model_glm_red_immun), g=10)
  print(glm7_hl)
  
  # -----------------------------------------
  #ROC and AUC
  # -----------------------------------------
  probabilities7 <- model_glm_red_immun %>% predict(myData, type = "response")
  predict_b7<-prediction(probabilities7, myData$Ulcer.healing.status)
  perf_b7<-performance(predict_b7, measure="tpr", x.measure="fpr")
  AUC_b7<-performance(predict_b7, measure="auc")
  #AUC
  AUC_b7@y.values
  #Plot ROC Curve
  tiff("Model 6 - Immune model 2.tiff", width = 700, height = 600)
  par(mar=c(5,5,5,5), cex.axis=1.5)
  par(pty="s")
  roc1<-roc(myData$Ulcer.healing.status,probabilities7, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
            main="Model 6 - Immune model 2", cex.lab=1.7, cex.main=1.4,type="o",
            xlab="False positive rate", ylab="True positive rate", font.lab=2)
  
  dev.off()
  # -----------------------------------------
  
  
  #-------------------------------------------------------------------------------------
  #MODEL8 - Stepwise reduced logistic regression model - immune parameters
  #-------------------------------------------------------------------------------------
  
  #REDUCED MODEL
  full.model2 <- glm(Ulcer.healing.status ~ 
                       MFI.CD15..CD16   +   MFI.CD15..CD62L  +   MFI.CD15..CD63   +         
                       MFI.CD15..CD282     +        MFI.CD15..CD284    +         MFI.CD15..HLA.DR  +
                       MFI.CD15..CD11b     +        MFI.CD14..CD16     +         MFI.CD14..CD62L   +
                       MFI.CD14..CD63      +        MFI.CD14..CD282    +         MFI.CD14..CD284   +
                       MFI.CD14..CD36         +     MFI.CD14..HLA.DR     +       MFI.CD14..CD11b  , 
                     data = myData, family = 'binomial', maxit=100)
  
  #STEP MODEL
  s_model2 <- full.model2 %>% stepAIC(trace = FALSE)
  smod_coef2<-coef(s_model2)
  print(smod_coef2)
  summary(s_model2)
  
  #ODDS Ratios
  exp(s_model2$coefficients)
  exp(confint(s_model2))
  
  #Goodness of fit - model evaluation statistics
  with(s_model2, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))
  
  #Hosmer lemeshow test
  glm8_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(s_model2), g=10)
  print(glm8_hl)
  
  # -----------------------------------------
  #ROC and AUC
  # -----------------------------------------
  probabilities8 <- s_model2 %>% predict(myData, type = "response")
  predict_b8<-prediction(probabilities8, myData$Ulcer.healing.status)
  perf_b8<-performance(predict_b8, measure="tpr", x.measure="fpr")
  AUC_b8<-performance(predict_b8, measure="auc")
  #AUC
  AUC_b8@y.values
  #Plot ROC Curve
  tiff("Model 7 - Immune model 3.tiff", width = 700, height = 600)
  par(mar=c(5,5,5,5), cex.axis=1.5)
  par(pty="s")
  roc1<-roc(myData$Ulcer.healing.status,probabilities8, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
            main="Model 7 - Immune model 3", cex.lab=1.7, cex.main=1.4,type="o",
            xlab="False positive rate", ylab="True positive rate", font.lab=2)
  
  dev.off()
  # -----------------------------------------  
  
  #-------------------------------------------------------------------------------------
  #MODEL9 - Combined model 1 - Selected parameters of step models
  #-------------------------------------------------------------------------------------
  model_glm_comb <- glm(Ulcer.healing.status~ Total.Cholesterol..mg.dl.+LDL.Cholesterol..mg.dl.+X..CD14.+X..CD15.High+
                          MFI.CD14..CD63+MFI.CD14..CD282+MFI.CD14..CD11b,
                             data = myData, family = 'binomial', maxit=100)
  summary(model_glm_comb)
  
  #ODDS Ratios
  exp(model_glm_comb$coefficients)
  exp(confint(model_glm_comb))
  
  #Goodness of fit - model evaluation statistics
  with(model_glm_comb, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))
  
  #Hosmer lemeshow test
  glm9_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(model_glm_comb), g=10)
  print(glm9_hl)
  
  # -----------------------------------------
  #ROC and AUC
  # -----------------------------------------
  probabilities9 <- model_glm_comb %>% predict(myData, type = "response")
  predict_b9<-prediction(probabilities9, myData$Ulcer.healing.status)
  perf_b9<-performance(predict_b9, measure="tpr", x.measure="fpr")
  AUC_b9<-performance(predict_b9, measure="auc")
  #AUC
  AUC_b9@y.values
  #Plot ROC Curve
  tiff("Model 8 - Combined model 1.tiff", width = 700, height = 600)
  par(mar=c(5,5,5,5), cex.axis=1.5)
  par(pty="s")
  roc1<-roc(myData$Ulcer.healing.status,probabilities9, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
            main="Model 8 - Combined model 1", cex.lab=1.7, cex.main=1.4,type="o",
            xlab="False positive rate", ylab="True positive rate", font.lab=2)
  
  dev.off()
  # -----------------------------------------
  
  #-------------------------------------------------------------------------------------
  #MODEL10 - Stepwise reduced logistic regression model - combined parameters
  #-------------------------------------------------------------------------------------
  
  #REDUCED MODEL
  full.model3 <- glm(Ulcer.healing.status ~ 
                       Total.Cholesterol..mg.dl.+LDL.Cholesterol..mg.dl.+X..CD14.+X..CD15.High+
                       MFI.CD14..CD63+MFI.CD14..CD282+MFI.CD14..CD11b , 
                     data = myData, family = 'binomial', maxit=100)
  
  #STEP MODEL
  s_model3 <- full.model3 %>% stepAIC(trace = FALSE)
  smod_coef3<-coef(s_model3)
  print(smod_coef3)
  summary(s_model3)
  
  #ODDS Ratios
  exp(s_model3$coefficients)
  exp(confint(s_model3))
  
  #Goodness of fit - model evaluation statistics
  with(s_model3, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))
  
  #Hosmer lemeshow test
  glm10_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(s_model3), g=10)
  print(glm10_hl)
  
  # -----------------------------------------
  #ROC and AUC
  # -----------------------------------------
  probabilities10 <- s_model3 %>% predict(myData, type = "response")
  predict_b10<-prediction(probabilities10, myData$Ulcer.healing.status)
  perf_b10<-performance(predict_b10, measure="tpr", x.measure="fpr")
  AUC_b10<-performance(predict_b10, measure="auc")
  #AUC
  AUC_b10@y.values
  #Plot ROC Curve
  tiff("Model 9 - Combined model 2.tiff", width = 700, height = 600)
  par(mar=c(5,5,5,5), cex.axis=1.5)
  par(pty="s")
  roc1<-roc(myData$Ulcer.healing.status,probabilities10, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
            main="Model 9 - Combined model 2", cex.lab=1.7, cex.main=1.4,type="o",
            xlab="False positive rate", ylab="True positive rate", font.lab=2)
  
  dev.off()
  # -----------------------------------------  
  
  #-------------------------------------------------------------------------------------
  #MODEL11 - Selected parameters of step models - taking sig one's of model 9
  #-------------------------------------------------------------------------------------
  model_glm_comb2 <- glm(Ulcer.healing.status~ LDL.Cholesterol..mg.dl.+MFI.CD14..CD63,
                        data = myData, family = 'binomial', maxit=100)
  summary(model_glm_comb2)
  
  #ODDS Ratios
  exp(model_glm_comb2$coefficients)
  exp(confint(model_glm_comb2))
  
  #Goodness of fit - model evaluation statistics
  with(model_glm_comb2, pchisq(null.deviance-deviance, df.null-df.residual, lower.tail = F))
  
  #Hosmer lemeshow test
  glm11_hl<-hoslem.test(myData$Ulcer.healing.status,fitted(model_glm_comb2), g=10)
  print(glm11_hl)
  
  # -----------------------------------------
  #ROC and AUC
  # -----------------------------------------
  probabilities11 <- model_glm_comb2 %>% predict(myData, type = "response")
  predict_b11<-prediction(probabilities11, myData$Ulcer.healing.status)
  perf_b11<-performance(predict_b11, measure="tpr", x.measure="fpr")
  AUC_b11<-performance(predict_b11, measure="auc")
  #AUC
  AUC_b11@y.values
  #Plot ROC Curve
  tiff("Model 10 - Combined model 3.tiff", width = 700, height = 600)
  par(mar=c(5,5,5,5), cex.axis=1.5)
  par(pty="s")
  roc1<-roc(myData$Ulcer.healing.status,probabilities11, plot=TRUE, col="blue", lwd=2,legacy.axes=TRUE,
            main="Model 10 - Combined model 3", cex.lab=1.7, cex.main=1.4,type="o",
            xlab="False positive rate", ylab="True positive rate", font.lab=2)
  
  dev.off()
  # -----------------------------------------
  
  #----------------------------------------------------------------------------------
  #PLOT RESULTS - MODEL FIT PLOTS
  #----------------------------------------------------------------------------------
  #theme definition
  JtSNEtheme<-theme(panel.background = element_rect(fill = "white", colour = "black"), 
                    legend.direction = 'vertical', 
                    legend.position = c(0.975, 0),
                    legend.justification=c(1,-0.45),
                    legend.background = element_blank(),
                    legend.key = element_blank(),
                    plot.title = element_text(size = 16, face = "bold", hjust=0.5),
                    axis.text = element_text(size = 12, colour = "black"),
                    axis.title=element_text(size=14,face="bold"),
                    legend.text = element_text(size = 12),
                    legend.title = element_text(size=14,face="bold"),
                    axis.ticks.length=unit(.25, "cm"))
  
  #Plot Results of model 9 
  p1<-data.frame(Probability_of_Healing=probabilities1, Healing_status=as.numeric(myData$Ulcer.healing.status), Org_Healing_status=as.factor(myData$Ulcer.healing.status))
  p1$Patient_ID<-factor(rownames(p1))
  p1<-p1[order(p1$Probability_of_Healing, decreasing = FALSE),]
  p1$Rank<-1:nrow(p1)
  a1<-ggplot(data=p1, aes(x=Rank, y=Probability_of_Healing, color=Org_Healing_status))+
    geom_point(alpha=1, shape=3, stroke=2)+
    scale_color_manual(labels = c("Healed ulcers", "Non-healed ulcers"), values = c("#00AFBB", "#FC4E07"))+
    ylim(0,1)+
    xlab("Patient Index")+
    ylab("Predicted probability of ulcer healing")+
    labs(title="CD14 CD63", colour="Healing status")+
    JtSNEtheme
  a1
  ggsave("CD14 CD63 - 01.tiff", width = 5, height = 6)
  a2<-a1+geom_point(aes(y=Healing_status), colour="grey50")
   a2
   
 #for saving vector files
 pdf("CD63 - Uni - Model fit.pdf",
       width=6,height=6,
       bg="white",
       colormodel="srgb")
 plot(a2)
   dev.off()     
   
   ####-----------------------------THE END-----------------------------####