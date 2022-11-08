library(openxlsx)
library(magrittr)
library(stringr)

setwd("Z:/Groups/Peder Olofsson/VSS/Bioconductor/HepatoCohort")

toRemove <- list()

#Create sets of interesting genes
cholGenes <- c("SREBF2", "HMGCR", "HMGCS1", "LDLR", "PCSK9", "NPC1L1", "CETP", "LCAT", "ABCG1") %>% factor(., levels = .)
tgGenes <- c("SREBF1", "ACLY", "ACACA", "ACACB", "SCD", "FAS", "MLXIPL", "NR1H2", "AGPAT2", "ANGPTL4", "FABP4", "MTTP") %>% factor(., levels = .)


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

exprDf <- allData %>% .[, !colnames(allData) %in% c(colnames(metadata), "rs738409.y", "rs738409.x", "blood_pressure_systolic_1", "blood_pressure_diastolic_1", 
                                                    "blood_pressure_systolic_2", "blood_pressure_diastolic_2")] %>% apply(., 2, as.double) %>% 
  .[, !colnames(.) %in% c("PNPLA3.1", "TM6SF2.1", "SREBF2.1", "PLIN2_A", "FAS_A", "SOAT2.1", "ACAT2_A", "SCAP_A", "HMGCR_A", 
                          "CIDEA_A", "ACACA_A", "CES2_A", "SREBF1.1", "ABCG1_A",  "ABCG1_B",  "ABCG1_C", "SYP.1")]
rownames(exprDf) <- str_c("P", metadata$ASAP_ID, sep = "")


#Clean the dataset:
#Remove donors with NA in statins
toRemove[["Statins"]] <- !is.na(metadata$Statins)

metadata %<>% .[toRemove[["Statins"]],]
exprDf %<>% .[toRemove[["Statins"]],]


#Remove patients who were fasted for less then 6h
toRemove[["Fasting"]] <- metadata$Fasting_time >= 6
toRemove[["Fasting"]][is.na(toRemove[["Fasting"]])] <- FALSE

metadata %<>% .[toRemove[["Fasting"]],]
exprDf %<>% .[toRemove[["Fasting"]],]


#Remove donors above ALAT threshold
toRemove[["ALAT"]] <- metadata$ALAT <= 1.1
toRemove[["ALAT"]][is.na(toRemove[["ALAT"]])] <- FALSE

metadata %<>% .[toRemove[["ALAT"]],]
exprDf %<>% .[toRemove[["ALAT"]],]


#Create a selection matrix:
selectionMatrix <- matrix(nrow = nrow(metadata), ncol = 5)
colnames(selectionMatrix) <- c("Waist_circumference", "blood_pressure_diastolic", "ClinData_TG", "ClinData_HDL-cholesterol", "ClinData_Glucose")

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
  if(is.na(metadata[index, "blood_pressure_systolic"]) == TRUE || is.na(metadata[index, "blood_pressure_diastolic"]) == TRUE) {
    NA
  } else {
    if(metadata[index, "blood_pressure_systolic"] >= 130 || metadata[index, "blood_pressure_diastolic"] >= 85) {
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


#Remove donors, where NA is in one of crucial criteria
selectionMatrix <- matrix(nrow = nrow(metadata), ncol = 5)
colnames(selectionMatrix) <- c("Waist_circumference", "blood_pressure_diastolic", "ClinData_TG", "ClinData_HDL-cholesterol", "ClinData_Glucose")

#Remove NAs everywhere
toRemove[["NAs"]] <- metadata %>% .[, colnames(selectionMatrix)] %>% is.na() %>% apply(., 1, function(patient) sum(patient) < 1)

metadata %<>% .[toRemove[["NAs"]],]
exprDf %<>% .[toRemove[["NAs"]],]
selectionMatrix %<>% .[toRemove[["NAs"]],]

summary(metadata$MS)


#Look at gene expression:
#MS vs non-MS without statin treatment
noStat <- metadata %>% .$Statins == "0"

msStat <- calculateExpression(exprDf[noStat,], metadata[noStat,], "MS")

#Export
exportData(msStat, metadata[noStat,], "MSvsNoMS_noStatins.xlsx")


#Statin vs non-Statin in MS
rownames(metadata) %v% rownames(exprDf)

ms <- metadata %>% .[.$MS == "MS",] %>% rownames(.)

metaMS <- metadata %>% .[rownames(.) %in% ms,]
exprMS <- exprDf %>% .[rownames(.) %in% ms,]

metaMS$Statins <- ifelse(metaMS$Statins == 1, "Statins", "No-Statins")

stNst <- calculateExpression(exprMS, metaMS, "Statins")

#Export
exportData(stNst, metaMS, "MS_StatVsNoStat.xlsx")