#--------R script written by JVR for chi-square analysis of histopathological data uploaded on ---------

library(vcd)
library(grid)
library(corrplot)
#Data reading 
myData=read.csv(file.choose(),header=T, row.names=1)
#Data summary
myData
colnames(myData)
rownames(myData)
str(myData)
#Factorise all variables
col_names <- names(myData)
myData[,col_names] <- lapply(myData[,col_names] , factor)
# myData$Granulation.Tissue<-as.factor(myData$Granulation.Tissue)

# check structure of dataframe after factorising
str(myData)

#------Chi-square analysis--------------

# ----------------------Healing status vs all variables-------------------------

#chi square test - Granulation tissue
var_1 <- table(myData$Ulcer.healing.status, myData$Granulation.Tissue) #order of variables specified does not matter!
print(var_1)
print(chisq.test(var_1))
print(chisq.test(var_1, correct=FALSE))
print(fisher.test(var_1))

#chi square test - Epidermal Thickening
var_2 <- table(myData$Ulcer.healing.status, myData$Epidermal.Thickening)
print(var_2)
print(chisq.test(var_2))
print(chisq.test(var_2, correct=FALSE))
print(fisher.test(var_2))

#chi square test - Vessel Thickening
var_3 <- table(myData$Ulcer.healing.status, myData$Vessels.Thickening)
print(var_3)
print(chisq.test(var_3))
print(chisq.test(var_3, correct=FALSE))
print(fisher.test(var_3))

#chi square test - Endothelial nuclear hypertrophy
var_4 <- table(myData$Ulcer.healing.status, myData$Endothelial.Nuclear.hypertrophy)
print(var_4)
print(chisq.test(var_4))
print(chisq.test(var_4, correct=FALSE))
print(fisher.test(var_4))

#chi square test - Collagen
var_5 <- table(myData$Ulcer.healing.status, myData$Collagen)
print(var_5)
print(chisq.test(var_5))
print(chisq.test(var_5, correct=FALSE))
print(fisher.test(var_5))

#chi square test - Fibrin
var_6 <- table(myData$Ulcer.healing.status, myData$Fibrin)
print(var_6)
print(chisq.test(var_6))
print(chisq.test(var_6, correct=FALSE))
print(fisher.test(var_6))

#chi square test - Karyorrhectic debri
var_7 <- table(myData$Ulcer.healing.status, myData$karyorrhectic.debris)
print(var_7)
print(chisq.test(var_7))
print(chisq.test(var_7, correct=FALSE))
print(fisher.test(var_7))

#chi square test - Koilocytosis
var_8 <- table(myData$Ulcer.healing.status, myData$Koilocytosis.Koilocytotic.changes)
print(var_8)
print(chisq.test(var_8))
print(chisq.test(var_8, correct=FALSE))
print(fisher.test(var_8))

#chi square test - Ulcer stage
var_9 <- table(myData$Ulcer.healing.status, myData$Ulcer.stage)
print(var_9)
print(chisq.test(var_9))
print(chisq.test(var_9, correct=FALSE))
print(fisher.test(var_9))

#chi square test - Neutrophils
var_10 <- table(myData$Ulcer.healing.status, myData$Neutrophils)
print(var_10)
print(chisq.test(var_10))
print(chisq.test(var_10, correct=FALSE))
print(fisher.test(var_10))

#chi square test - Lymphocytes
var_11 <- table(myData$Ulcer.healing.status, myData$Lymphocytes.Plasmocytes)
print(var_11)
print(chisq.test(var_11))
print(chisq.test(var_11, correct=FALSE))
print(fisher.test(var_11))



# ----------------------Stage vs all variables-------------------------
#chi square test - Granulation tissue
var_1_stg <- table(myData$Ulcer.stage, myData$Granulation.Tissue) #order of variables specified does not matter!
print(var_1_stg)
print(chisq.test(var_1_stg))
print(chisq.test(var_1_stg, correct=FALSE))
print(fisher.test(var_1_stg))

#chi square test - Epidermal Thickening
var_2_stg <- table(myData$Ulcer.stage, myData$Epidermal.Thickening)
print(var_2_stg)
print(chisq.test(var_2_stg))
print(chisq.test(var_2_stg, correct=FALSE))
print(fisher.test(var_2_stg))

#chi square test - Vessel Thickening
var_3_stg <- table(myData$Ulcer.stage, myData$Vessels.Thickening)
print(var_3_stg)
print(chisq.test(var_3_stg))
print(chisq.test(var_3_stg, correct=FALSE))
print(fisher.test(var_3_stg))

#chi square test - Endothelial nuclear hypertrophy
var_4_stg <- table(myData$Ulcer.stage, myData$Endothelial.Nuclear.hypertrophy)
print(var_4_stg)
print(chisq.test(var_4_stg))
print(chisq.test(var_4_stg, correct=FALSE))
print(fisher.test(var_4_stg))

#chi square test - Collagen
var_5_stg <- table(myData$Ulcer.stage, myData$Collagen)
print(var_5_stg)
print(chisq.test(var_5_stg))
print(chisq.test(var_5_stg, correct=FALSE))
print(fisher.test(var_5_stg))

#chi square test - Fibrin
var_6_stg <- table(myData$Ulcer.stage, myData$Fibrin)
print(var_6_stg)
print(chisq.test(var_6_stg))
print(chisq.test(var_6_stg, correct=FALSE))
print(fisher.test(var_6_stg))

#chi square test - Karyorrhectic debri
var_7_stg <- table(myData$Ulcer.stage, myData$karyorrhectic.debris)
print(var_7_stg)
print(chisq.test(var_7_stg))
print(chisq.test(var_7_stg, correct=FALSE))
print(fisher.test(var_7_stg))

#chi square test - Koilocytosis
var_8_stg <- table(myData$Ulcer.stage, myData$Koilocytosis.Koilocytotic.changes)
print(var_8_stg)
print(chisq.test(var_8_stg))
print(chisq.test(var_8_stg, correct=FALSE))
print(fisher.test(var_8_stg))

#chi square test - Neutrophils
var_10_stg <- table(myData$Ulcer.stage, myData$Neutrophils)
print(var_10_stg)
print(chisq.test(var_10_stg))
print(chisq.test(var_10_stg, correct=FALSE))
print(fisher.test(var_10_stg))

#chi square test - Lymphocytes
var_11_stg <- table(myData$Ulcer.stage, myData$Lymphocytes.Plasmocytes)
print(var_11_stg)
print(chisq.test(var_11_stg))
print(chisq.test(var_11_stg, correct=FALSE))
print(fisher.test(var_11_stg))
#------------------------------
#PLOT RESULTS
#------------------------------

#mosaic plot
tiff("Filtered patients 2 - Gran tissue vs stage.tiff", width = 700, height = 600)
mosaic(~ Ulcer.stage + Granulation.Tissue,
       direction = c("h", "v"),
       data = myData,
       shade = TRUE
)
dev.off()
#Corrplot

#Pearson's residuals
var_1_stgres <- chisq.test(var_1_stg)
tiff("Filtered patients 2 - Gran tissue vs stage_Corr.tiff", width = 700, height = 600)
corrplot(var_1_stgres$residuals, is.cor = FALSE)
dev.off()
# Contibution in percentage (%)
contrib <- 100*var_1_stgres$residuals^2/var_1_stgres$statistic
round(contrib, 3)
corrplot(contrib, is.cor = FALSE)

####-----------------------------THE END-----------------------------####
