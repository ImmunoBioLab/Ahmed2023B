library(openxlsx)
library(magrittr)
library(stringr)
library(ggplot2)
library(showtext)


setwd("Z:/Groups/Peder Olofsson/VSS/Bioconductor/HepatoCohort")

#Create parameter vectors
categorical <- c("Sex", "Regular_smoker")
numerical <- c("Age",  "BMI", "BSA_DuBois", "Insulin", "ApoB", "ApoA1", "Lp(a)", "Height", "Weight", "Waist_circumference", "Hip_circumference", "Fasting_time", "ClinData_Glucose", 
               "ClinData_TG", "ClinData_Cholesterol", "ClinData_HDL-cholesterol", "ClinData_LDL-cholesterol", "ClinData_HbA1c", "ALAT")


#Read data
dataDf <- readr::read_tsv("PNPLA3.Osman.HuEx.expression.asap.liver.txt") %>% as.data.frame()
rownames(dataDf) <- paste0("P", dataDf[[1]])
dataDf %<>% .[, c(2:ncol(dataDf))]


#Read metadata
metadata <- read.xlsx("ProcessedDataset.xlsx", sheet = "ProcessedMetadata", rowNames = TRUE)
metadata %<>% .[rownames(.) %in% rownames(dataDf),]

dataDf %<>% .[match(rownames(metadata), rownames(.)),]


#Evaluate gene expression:
#MS vs non-MS without statin treatment
noStat <- metadata %>% .$corStatin == "No Statins"

msStat <- processDataset(dataDf[noStat,], metadata[noStat,], "MS", list(categorical, numerical))

#Export
exportDataset(msStat, "Alt_MSvsNoMS_noStatins.xlsx")


#MS vs non-MS without statin treatment: excluding diabetes
noStatNoDiab <- metadata$corStatin == "No Statins" & metadata$Diabetes == "No diabetes"
noStatNoDiab[is.na(noStatNoDiab) == TRUE] <- FALSE

msStatNoDiab <- processDataset(dataDf[noStatNoDiab,], metadata[noStatNoDiab,], "MS", list(categorical, numerical))

#Export
exportDataset(msStatNoDiab, "Alt_MSvsNoMS_noStatins_noDiabetes.xlsx")


#Statin vs non-Statin:
#Whole dataset:
whStat <- processDataset(dataDf, metadata, "corStatin", list(categorical, numerical))

#Export
exportDataset(whStat, "Alt_StatVsNoStat_wholeDataset.xlsx")


#In MS
ms <- metadata %>% .$MS == "MS"

statMS <- processDataset(dataDf[ms,], metadata[ms,], "corStatin", list(categorical, numerical))

#Export
exportDataset(statMS, "Alt_StatVsNoStat_MS.xlsx")


#In no-MS
nonMS <- metadata %>% .$MS == "Non-MS"

statNonMS <- processDataset(dataDf[nonMS,], metadata[nonMS,], "corStatin", list(categorical, numerical))

#Export
exportDataset(statNonMS, "Alt_StatVsNoStat_non-MS.xlsx")