#Functions

#Merge lists into data.frames where names of list components form a new column (with name from newColName)
pcrMerge <- function(object, newColName = "Experiment", renewIds = TRUE, saveAllCols = TRUE, oblColumns = NULL, newId = TRUE) {
  #Merges list of data.frames into one data.frame:
  #object - named list where names have meaning;
  #newColName - name for a new column that will be created from list names;
  #renewIds - if Id column is present in components - make a new one;
  #saveAllCols - preserve all columns from all data.frames;
  #oblColumns - character vector of all column names to preserve - the rest are removed.
  
  if(is.null(names(object)) == TRUE) {
    stop("Merging impossible - object lacks names!")
  } else {
    #Manage columns
    if(is.null(oblColumns) == TRUE) {
      oblColumns <- lapply(object, function(experiment) colnames(experiment))
      
      #Keep all columns or only the common ones
      if(saveAllCols == TRUE) {
        oblColumns %<>% unlist() %>% unique()
      } else {
        allColNames <- oblColumns %>% unlist() %>% unique()
        
        namesToRemove <- lapply(oblColumns, function(experiment) allColNames %>% .[!. %in% experiment]) %>% unlist()
        
        oblColumns %<>% unlist() %>% unique() %>% .[!. %in% namesToRemove]
      }
    } else {
      #Remove unwanted columns
      object %<>% lapply(., function(experiment) experiment %<>% .[, colnames(.) %in% oblColumns])
    }
    
    #Check if any data.frames lack columns
    namesToAdd <- lapply(object, function(experiment) oblColumns %>% .[!. %in% colnames(experiment)]) %>% unlist() %>% length()
    
    if(namesToAdd > 0) {
      object %<>% lapply(., function(experiment) {
        absentColumns <- sum(!oblColumns %in% colnames(experiment))
        
        if(absentColumns > 0) {
          for(i in 1:absentColumns) {
            experiment <- cbind(experiment, NA)
          }
          colnames(experiment) <- c(colnames(experiment)[!colnames(experiment) %in% "NA"], 
                                    oblColumns[!oblColumns %in% colnames(experiment)])
        }
        
        return(experiment)
      })
    }
    
    #Extract list component names
    experimentNames <- names(object)
    
    #Add component names to the data.frames
    object <- lapply(seq_along(object), function(index) {
      object[[index]] %<>% dplyr::mutate(., NewCol = names(object)[index]) 
    })
    object %<>% do.call("rbind", .)
    colnames(object)[colnames(object) == "NewCol"] <- newColName
    
    #Put new column at the front
    object %<>% .[,c(newColName, colnames(object)[c(1:ncol(object)-1)])]
    
    if("Id" %in% colnames(object)) {
      if(newId == TRUE) {
        object$Id <- c(1:nrow(object))
      }
      
      allCols <- colnames(object) %>% .[!. %in% "Id"]
      object %<>% .[,c("Id", allCols)]
    }
    
    return(object)
  }
}

#Format stat labels for graphs
setGeneric("pcrFormatStatLabels", function(object, formatType = NULL, statThresholds = c("*" = 0.05, "**" = 0.01)) standardGeneric("pcrFormatStatLabels"))
setMethod("pcrFormatStatLabels", signature(object = "vector"), function(object, formatType = NULL, 
                                                                        statThresholds = c("*" = 0.05, "**" = 0.01)) {
  #Turn vector into numeric
  if(str_detect(object[1], "p = ") == TRUE) {
    object %<>% str_remove(., "p = ") %>% as.double()
  }
  
  if(is.null(formatType) == TRUE || !formatType %in% c("Stars", "Numeric", "Vs p")) {
    cat("Choose stat label format:", "\n")
    formatType <- select.list(c("Stars", "Numeric", "Vs p"))
  }
  
  if(formatType == "Stars") {
    #Add non-signifncant threshold if necessary
    if(!"ns" %in% names(statThresholds)) {
      statThresholds <- c("ns" = 1, statThresholds)
    }
    
    #Find the lowest possible threshold for each value and get appropriate labels
    statLabels <- sapply(object, function(value) {
      sapply(statThresholds, function(pValue) value <= pValue)
    }) %>% apply(., 2, sum) %>%
      names(statThresholds)[.]
  } else if(formatType == "Numeric") {
    #Format labels in a scientific manner
    statLabels <- ifelse(object > 0.01,
                         round(object, 2),
                         format(object, digits = 2, scientific = TRUE) 
    )
    statLabels <- str_c("p = ", statLabels)
  } else if(formatType == "Vs p") {
    #Add non-signifncant threshold if necessary
    if(!"ns" %in% names(statThresholds)) {
      statThresholds <- c("ns" = 1, statThresholds)
    }
    
    #Find the lowest possible threshold for each value and get appropriate labels
    statLabels <- sapply(object, function(value) {
      sapply(statThresholds, function(pValue) value < pValue)
    }) %>% apply(., 2, sum)
    statLabels <- str_c("p < ", statThresholds[statLabels])
    statLabels[statLabels == "p < 1"] <- "ns"
  }
  
  return(statLabels)
})
setMethod("pcrFormatStatLabels", signature(object = "data.frame"), function(object, formatType = NULL, 
                                                                            statThresholds = c("*" = 0.05, "**" = 0.01)) {
  #object - dat.frame, an output of pcrMakeStatLabels
  object$p <- object$Label
  object$Label <- pcrFormatStatLabels(object$Label, formatType = formatType, statThresholds = statThresholds)
  
  return(object)
})

#Are all instances of vector X present in vector Y?
setGeneric("%v%", function(x, y) standardGeneric("%v%"))
setMethod("%v%", signature(x = "vector", y = "vector"), function(x, y) {
  all(x %in% y)
})
setMethod("%v%", signature(x = "vector", y = "list"), function(x, y) {
  all(sapply(y, function(z) x %v% z))
})

#Calcualte stuff
calculateExpression <- function(exprMat, metaDf, compCol) {
  sumList <- lapply(colnames(exprMat), function(gene) {
    exprGene <- exprMat[, gene, drop = FALSE]
    colnames(exprGene) <- "Value"
    
    exprGene <- cbind(metaDf[, compCol, drop = FALSE], exprGene)
    colnames(exprGene)[colnames(exprGene) == compCol] <- "Comparison"
    
    sumDf <- dplyr::group_by_at(exprGene, "Comparison") %>%
      dplyr::summarise(Mean = mean(Value, na.rm = TRUE),
                       SD = sd(Value, na.rm = TRUE),
                       Var = var(Value, na.rm = TRUE),
                       Count = sum(complete.cases(Value))
      )
    
    sumDf %<>% dplyr::mutate(SEM = .$SD/sqrt(.$Count))
    
    statDf <- ggpubr::compare_means(Value ~ Comparison, exprGene)
    
    list(Expression = exprGene,
         Summary = sumDf,
         Stats = statDf
    )
  })
  names(sumList) <- colnames(exprMat)
  
  return(sumList)
}

#Export data
exportData <- function(dataList, metaDf, fileName) {
  exprMat <- lapply(dataList, function(gene) gene[["Expression"]]$Value) %>%
    do.call("cbind", .)
  colnames(exprMat) <- names(dataList)
  rownames(exprMat) <- rownames(metaDf)
  
  sumData <- lapply(dataList, function(gene) gene$Summary) %>%
    pcrMerge(., "Gene")
  
  statData <- lapply(dataList, function(gene) gene$Stats) %>%
    pcrMerge(., "Gene")
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "Metadata")
  writeData(wb, 1, metaDf)

  addWorksheet(wb, "Expression")
  writeData(wb, 2, exprMat, rowNames = TRUE)
  
  addWorksheet(wb, "Results")
  writeData(wb, 3, sumData)
  
  addWorksheet(wb, "Stats")
  writeData(wb, 4, statData)
  
  saveWorkbook(wb, file.path(getwd(), fileName))
}
