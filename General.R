library(openxlsx)
library(magrittr)
library(stringr)
library(ggplot2)
library(showtext)
font_add(family = "docktrin", regular = "C:/Users/vlasha/Downloads/docktrin/docktrin.ttf")
showtext_auto(enable = TRUE, record = TRUE)


setwd("Z:/Groups/Peder Olofsson/VSS/Bioconductor/HepatoCohort")

toRemove <- list()

#Create sets of interesting genes
cholGenes <- c("SREBF2", "HMGCR", "HMGCS1", "LDLR", "PCSK9", "NPC1L1", "CETP", "LCAT", "ABCG1") %>% factor(., levels = .)
tgGenes <- c("SREBF1", "ACLY", "ACACA", "ACACB", "SCD", "FAS", "MLXIPL", "NR1H2", "AGPAT2", "ANGPTL4", "FABP4", "MTTP") %>% factor(., levels = .)

#Create parameter vectors
categorical <- c("Sex", "Regular_smoker")
numerical <- c("Age",  "BMI", "BSA_DuBois", "Insulin", "ApoB", "ApoA1", "Lp(a)", "Height", "Weight", "Waist_circumference", "Hip_circumference", "Fasting_time", "ClinData_Glucose", 
                "ClinData_TG", "ClinData_Cholesterol", "ClinData_HDL-cholesterol", "ClinData_LDL-cholesterol", "ClinData_HbA1c", "ALAT")


#Combine all files and create sensible data structure
allData <- read.xlsx("Data_2022-10-31/20220404 All data set updated.xlsx")
addGenes <- read.xlsx("Data_2022-10-31/220919_Additional genes PNPLA3 project.xlsx", sheet = "Usable")

allData %<>% dplyr::left_join(., addGenes, "ASAP_ID")

mTORgenes <- read.xlsx("Data_2022-10-31/mTOR genes_OA_221011.xlsx", rowNames = TRUE)
mTORgenes %<>% t() %>% as.data.frame() %>%
  dplyr::mutate(ASAP_ID = rownames(.))
mTORgenes %<>% apply(., 2, as.numeric) 
mTOR <- colnames(mTORgenes) %>% factor(., levels = .)

allData %<>% dplyr::left_join(., mTORgenes, "ASAP_ID", copy = TRUE)

addMeds <- read.xlsx("Data_2022-10-31/Current medications.xlsx")
colnames(addMeds)[1] <- "ASAP_ID"

allData %<>% dplyr::left_join(., addMeds, "ASAP_ID")

SNIPs <- read.xlsx("Data_2022-10-31/All SNIPS.xlsx", sheet = "All SNIPS")

allData %<>% dplyr::left_join(., SNIPs, "LIVER.ID")

ALAT <- read.xlsx("Data_2022-10-31/ALAT_ASAP_osman.xlsx")
ALAT$ALAT %<>% as.numeric()

allData %<>% dplyr::left_join(., ALAT, "ASAP_ID")


#Marking cohort based on criteria: MS/non-MS
#Waist circumferense 101 for men and 90 for women
#BP 130/85 or higher or taking blood pressure medication
#TG above 150 mg/dl (1.7 mmol/L)
#HDL less than 40 mg/dl (1.03 mmol/l) (men) or under 50 mg/dl (women)
#Fasting blood glucose (sugar) level greater than 5.6 mmol/L (100 mg/dl) or are taking glucose-lowering medications

metadata <- allData %>% .[, colnames(allData) %in% c("ASAP_ID", "Sex", "Age", "BMI", "BSA_DuBois", "Hypercholesterolemia", "Statins", "Regular_smoker", "Insulin",
                                                     "ApoB", "ApoA1", "Lp(a)", "blood_pressure_systolic_1",  "blood_pressure_diastolic_1",
                                                     "blood_pressure_systolic_2", "blood_pressure_diastolic_2", "Height", "Weight", "Waist_circumference",
                                                     "Hip_circumference", "Fasting_time", "ClinData_Glucose", "ClinData_TG", "ClinData_Cholesterol",
                                                     "ClinData_HDL-cholesterol", "ClinData_LDL-cholesterol", "ClinData_HbA1c_(IFCC)", "ClinData_HbA1c", "ASAP.ID",
                                                     "LIVER.ID", "rs738409", "Statins_OA", "ID_ASMADA", "Statins_ASMADA", "GENEPOOL_ID", "GENEPOOL_Current_Medication__Medication,_type",
                                                     "GENEPOOL_Current_Medication__Preparation_and_dose", "zocord", "Lipitor", "crestor", "simvastatin", "rs641738", "rs780094", "rs58542926",
                                                     "rs738409.y", "ALAT")]
rownames(metadata) <- str_c("P", metadata$ASAP_ID, sep = "")
colnames(metadata)[colnames(metadata) == "rs738409.y"] <- "rs738409"
colnames(metadata)[colnames(metadata) == "rs738409.x"] <- "rs738409"

metadata[["blood_pressure_diastolic"]] <- metadata %>% .[, c("blood_pressure_diastolic_1", "blood_pressure_diastolic_2")] %>% rowMeans(., na.rm = TRUE)
metadata[["blood_pressure_systolic"]] <- metadata %>% .[, c("blood_pressure_systolic_1", "blood_pressure_systolic_2")] %>% rowMeans(., na.rm = TRUE)

metadata %<>% .[, !colnames(metadata) %in% c("blood_pressure_diastolic_1", "blood_pressure_diastolic_2", "blood_pressure_systolic_1", "blood_pressure_systolic_2")]

metadata[["Regular_smoker"]][metadata$Regular_smoker == 0] <- "Non-smoker"
metadata[["Regular_smoker"]][metadata$Regular_smoker == 1] <- "Smoker"

exprDf <- allData %>% .[, !colnames(allData) %in% c(colnames(metadata), "rs738409.y", "rs738409.x", "blood_pressure_systolic_1", "blood_pressure_diastolic_1", 
                                                    "blood_pressure_systolic_2", "blood_pressure_diastolic_2")] %>% apply(., 2, as.double) %>% 
  .[, !colnames(.) %in% c("PNPLA3.1", "TM6SF2.1", "SREBF2.1", "PLIN2_A", "FAS_A", "SOAT2.1", "ACAT2_A", "SCAP_A", "HMGCR_A", 
                          "CIDEA_A", "ACACA_A", "CES2_A", "SREBF1.1", "ABCG1_A",  "ABCG1_B",  "ABCG1_C", "SYP.1")]
rownames(exprDf) <- str_c("P", metadata$ASAP_ID, sep = "")

#Gather diabetes data
metadata[["GENEPOOL_Current_Medication__Preparation_and_dose"]][metadata$`GENEPOOL_Current_Medication__Medication,_type` == "None"] <- "None"
diabetes <- metadata %>% .[, c("GENEPOOL_Current_Medication__Medication,_type", "GENEPOOL_Current_Medication__Preparation_and_dose")] %>% 
  apply(., 1, function(patient) str_c(patient, collapse = " "))
diabetes %<>% str_to_lower(.) %>% sapply(., function(patient) {
  sapply(c("insulin", "metformin", "glucophage", "mindiab"), function(medication) str_detect(patient, medication))
}) %>% t()
metadata$Diabetes <- apply(diabetes, 1, function(patient) sum(patient) > 0)
metadata$Diabetes <- ifelse(metadata$Diabetes == TRUE, "Diabetes", "No diabetes") %>% factor(., levels = c("No diabetes", "Diabetes"))

summary(metadata$Diabetes)


#Gather statin data
metadata$corStatin <- metadata %>% .[, colnames(.) %in% c("zocord", "Lipitor", "crestor", "simvastatin")] %>% apply(., 1, function(patient) sum(patient, na.rm = TRUE) > 0)
metadata$corStatin <- ifelse(metadata$corStatin == TRUE, "Statins", "No Statins") %>% factor(., levels = c("Statins", "No Statins"))

summary(metadata$corStatin)


#Clean the dataset:
#Remove donors with NA in statins
toRemove[["Statins"]] <- !(is.na(metadata$Statins) & metadata$corStatin == "No Statins")
names(toRemove[["Statins"]]) <- rownames(metadata)

metadata %<>% .[toRemove[["Statins"]],]
exprDf %<>% .[toRemove[["Statins"]],]


#Remove patients who were fasted for less then 6h
toRemove[["Fasting"]] <- metadata$Fasting_time >= 6
toRemove[["Fasting"]][is.na(toRemove[["Fasting"]])] <- FALSE
names(toRemove[["Fasting"]]) <- rownames(metadata)

metadata %<>% .[toRemove[["Fasting"]],]
exprDf %<>% .[toRemove[["Fasting"]],]


#Remove donors above ALAT threshold
toRemove[["ALAT"]] <- metadata$ALAT <= 1.1
toRemove[["ALAT"]][is.na(toRemove[["ALAT"]])] <- FALSE
names(toRemove[["ALAT"]]) <- rownames(metadata)

metadata %<>% .[toRemove[["ALAT"]],]
exprDf %<>% .[toRemove[["ALAT"]],]


#Create a selection matrix:
selectionMatrix <- matrix(nrow = nrow(metadata), ncol = 5)
colnames(selectionMatrix) <- c("Waist_circumference", "blood_pressure_diastolic", "ClinData_TG", "ClinData_HDL-cholesterol", "ClinData_Glucose")
rownames(selectionMatrix) <- rownames(metadata)

#Waist_circumference: male >= 100, female >= 90
selectionMatrix[,1] <- sapply(seq_along(metadata$ASAP_ID), function(index) {
  if(metadata$Sex[index] == "M") {
    ifelse(metadata[index, "Waist_circumference"] >= 100, TRUE, FALSE)
  } else {
    ifelse(metadata[index, "Waist_circumference"] >= 90, TRUE, FALSE)
  }
})

#Blood pressure: systolic >= 130 or diastolic >= 85
selectionMatrix[,2] <- sapply(seq_along(metadata$ASAP_ID), function(index) {
  if(is.na(metadata[index, "blood_pressure_systolic"]) == TRUE && is.na(metadata[index, "blood_pressure_diastolic"]) == TRUE) {
    NA
  } else {
    if(metadata[index, "blood_pressure_systolic"] >= 130 && metadata[index, "blood_pressure_diastolic"] >= 85) {
      TRUE
    } else {
      FALSE
    } 
  }
})

#TG: > 1.7
selectionMatrix[,3] <- sapply(seq_along(metadata$ASAP_ID), function(index) {
  ifelse(metadata[index, "ClinData_TG"] > 1.7, TRUE, FALSE)
})

#HDL cholesterol: male < 1.03, female < 1.293
selectionMatrix[,4] <- sapply(seq_along(metadata$ASAP_ID), function(index) {
  if(metadata$Sex[index] == "M") {
    ifelse(metadata[index, "ClinData_HDL-cholesterol"] < 1.03, TRUE, FALSE)
  } else {
    ifelse(metadata[index, "ClinData_HDL-cholesterol"] < 1.293, TRUE, FALSE)
  }
})

#Glucose: > 5.6
selectionMatrix[,5] <- sapply(seq_along(metadata$ASAP_ID), function(index) {
  ifelse(metadata[index, "ClinData_Glucose"] > 5.6, TRUE, FALSE)
})

apply(selectionMatrix, 2, summary)


#Mark the donors according to the selection matrix
metadata$MS <- rowSums(selectionMatrix, na.rm = TRUE) > 2
metadata$MS <- ifelse(metadata$MS == TRUE, ifelse(selectionMatrix[, "Waist_circumference"] == TRUE, "MS", "semi-MS"), "Non-MS") %>% factor(., levels = c("Non-MS", "semi-MS", "MS"))


#Record
wb <- createWorkbook()
addWorksheet(wb, "OriginalMetadata")
writeDataTable(wb, "OriginalMetadata", metadata, rowNames = TRUE)
addWorksheet(wb, "OriginalExpression")
writeDataTable(wb, "OriginalExpression", as.data.frame(exprDf), rowNames = TRUE)
addWorksheet(wb, "SelectionMatrix")
writeDataTable(wb, "SelectionMatrix", as.data.frame(selectionMatrix))
saveWorkbook(wb, file.path(getwd(), "Dataset.xlsx"))


#Remove donors without waist circumference data
toRemove[["NoWaist"]] <- metadata %>% .[, "Waist_circumference"] %>% is.na(.) %>% !.
names(toRemove[["NoWaist"]]) <- rownames(metadata)

metadata %<>% .[toRemove[["NoWaist"]],]
exprDf %<>% .[toRemove[["NoWaist"]],]
selectionMatrix %<>% .[toRemove[["NoWaist"]],]


#Remove healthy donors with big waist + one more significant criteria with at least one NA
toRemove[["NAs"]] <- metadata %>% .[.$MS == "Non-MS",] %>% rownames()
toRemove[["NAs"]] <- selectionMatrix %>% .[rownames(.) %in% toRemove[["NAs"]],] %>% .[.[, "Waist_circumference"] == TRUE,]
toRemove[["NAs"]] <- apply(toRemove[["NAs"]], 1, function(patient) sum(is.na(patient)) > 1) %>% .[. == TRUE] %>% names()


summary(metadata$MS)


#Record the Removed table
toRemove %<>% .[sapply(., function(reason) length(reason) > 0)]

toExport <- data.frame(Patient = names(toRemove[[1]]),
                       Statins = NA,
                       Fasting = NA,
                       ALAT = NA,
                       NoWaist = NA
                       )
for(i in seq_along(toRemove)) {
  toExport[toExport$Patient %in% names(toRemove[[i]]), i+1] <- toRemove[[i]]
}

wb <- createWorkbook()
addWorksheet(wb, "Removed")
writeDataTable(wb, "Removed", toExport)
saveWorkbook(wb, file.path(getwd(), "Removed.xlsx"))


#Find Non-MS patients that fullfill none of the MS criteria - not MS in any way
nonMS <- rownames(metadata)[metadata$MS == "Non-MS"]

healthyMatrix <- selectionMatrix %>% .[rownames(.) %in% nonMS,] %>% .[apply(., 1, function(patient) sum(patient) == 0),]

#Assign
metadata$SuperHealth <- metadata$MS %>% factor(., levels = c("SH", "Non-MS", "semi-MS", "MS"))
metadata[rownames(metadata) %in% rownames(healthyMatrix), "SuperHealth"] <- "SH"

summary(metadata$SuperHealth)

#Record
wb <- createWorkbook()
addWorksheet(wb, "ProcessedMetadata")
writeDataTable(wb, "ProcessedMetadata", metadata, rowNames = TRUE)
addWorksheet(wb, "ProcessedExpression")
writeDataTable(wb, "ProcessedExpression", as.data.frame(exprDf), rowNames = TRUE)
addWorksheet(wb, "SelectionMatrix")
writeDataTable(wb, "SelectionMatrix", as.data.frame(selectionMatrix))
saveWorkbook(wb, file.path(getwd(), "ProcessedDataset.xlsx"))


#Evaluate gene expression:
#MS vs non-MS without statin treatment
noStat <- metadata %>% .$corStatin == "No Statins"

msStat <- processDataset(exprDf[noStat,], metadata[noStat,], "MS", list(categorical, numerical))

#Export
exportDataset(msStat, "MSvsNoMS_noStatins.xlsx")

#Plot
plotBox(msStat, geneNames = c("ACLY", "ACACA", "FAS", "SCD"), "MS", color = c("slategray1","slategray3", "slategray4")) %>%
  pcrSaveGraph(., getwd(), "Plots", "MSvsNoMS_noStatins_FA.tiff", 4, 12)


#MS vs non-MS without statin treatment: excluding diabetes
noStatNoDiab <- metadata$corStatin == "No Statins" & metadata$Diabetes == "No diabetes"
noStatNoDiab[is.na(noStatNoDiab) == TRUE] <- FALSE

msStatNoDiab <- processDataset(exprDf[noStatNoDiab,], metadata[noStatNoDiab,], "MS", list(categorical, numerical))

#Export
exportDataset(msStatNoDiab, "MSvsNoMS_noStatins_noDiabetes.xlsx")


#Statin vs non-Statin:
#Whole dataset:
whStat <- processDataset(exprDf, metadata, "corStatin", list(categorical, numerical))

#Export
exportDataset(whStat, "StatVsNoStat_wholeDataset.xlsx")


#In MS
ms <- metadata %>% .$MS == "MS"

statMS <- processDataset(exprDf[ms,], metadata[ms,], "corStatin", list(categorical, numerical))

#Export
exportDataset(statMS, "StatVsNoStat_MS.xlsx")


#In no-MS
nonMS <- metadata %>% .$MS == "Non-MS"

statNonMS <- processDataset(exprDf[nonMS,], metadata[nonMS,], "corStatin", list(categorical, numerical))

#Export
exportDataset(statNonMS, "StatVsNoStat_non-MS.xlsx")