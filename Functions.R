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

#Calculate expression
calculateExpression <- function(exprMat, metaDf, compCol, geneNames = NULL) {
  if(is.null(geneNames) == FALSE) {
    if(!geneNames %v% colnames(exprMat)) {
      stop("Gene(s) are absent from the exprMat!")
    } else {
     exprMat %<>% .[, geneNames, drop = FALSE] 
    }
  }
  
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
  
  exprMat <- lapply(sumList, function(gene) gene[["Expression"]]$Value) %>%
    do.call("cbind", .)
  colnames(exprMat) <- names(sumList)
  rownames(exprMat) <- rownames(metaDf)
  
  sumData <- lapply(sumList, function(gene) gene$Summary) %>%
    pcrMerge(., "Gene")
  
  statData <- lapply(sumList, function(gene) gene$Stats) %>%
    pcrMerge(., "Gene")
  
  statData$p.adj <- p.adjust(statData$p, "fdr")
  
  sumList <- list(Expression = exprMat,
                  Summary = sumData,
                  Stats = statData
                  )
  
  return(sumList)
}

#Calculate clinical parameters
calculateParameters <- function(metaDf, compCol = NULL, categorical = FALSE, parameters = NULL) {
  #Otherwise compare_means misbehaves
  parameters %<>% str_remove(., "\\(")
  parameters %<>% str_remove(., "\\)")
  colnames(metaDf) %<>% str_remove(., "\\(")
  colnames(metaDf) %<>% str_remove(., "\\)")
  parameters %<>% str_replace(., "-", "_")
  parameters %<>% str_remove(., "ClinData_")
  colnames(metaDf) %<>% str_replace(., "-", "_")
  colnames(metaDf) %<>% str_remove(., "ClinData_")

  if(is.null(compCol) == TRUE) {
    stop("Choose a column to compare by!")
  } else if(!compCol %in% colnames(metaDf)) {
    stop("Chosen comparison column is not present in metaDf!")
  } else if(!parameters %v% colnames(metaDf)) {
    stop("Chosen parameter(s) are not in the colnames of metaDf!")
  } else {
    colnames(metaDf)[colnames(metaDf) == compCol] <- "Comparison"
    
    if(categorical == FALSE) {
      sumDf <- dplyr::group_by_at(metaDf, "Comparison") %>%
        dplyr::summarise_at(., .vars = parameters, list(Mean = mean, SD = sd, Var = var), na.rm = TRUE)
      
      statDf <- metaDf[, c("Comparison", parameters)]
      
      statDf <- lapply(parameters, function(parameter) {
        print(parameter)
        ggpubr::compare_means(as.formula(paste(parameter, "~", "Comparison")), data = metaDf)
      }) %>% do.call("rbind", .)
      
      statDf$p.adj <- p.adjust(statDf$p, "fdr")
      
      sumList <- list(ParamSummary = sumDf,
                      ParamStats = statDf)  
    } else {
      sumList <- lapply(parameters, function(parameter) table(metadata[, c(parameter, compCol)])) %>% do.call("rbind", .)
    }
    
    return(sumList)
  }
}

#Calculate whole dataset
processDataset <- function(exprMat, metaDf, compCol = NULL, parameters = list(categorical = NULL, numerical = NULL), geneNames = NULL) {
  cat("Processing expression...", "\n")
  exprSum <- calculateExpression(exprMat = exprMat, metaDf = metaDf, compCol = compCol, geneNames = geneNames)
  
  cat("Processing numeric clinical parameters...", "\n")
  numParams <- calculateParameters(metaDf = metaDf, compCol = compCol, categorical = FALSE, parameters = parameters[[2]])
  
  cat("Processing categorical clinical parameters...", "\n")
  catParams <- calculateParameters(metaDf = metaDf, compCol = compCol, categorical = TRUE, parameters = parameters[[1]])
  
  dataList <- list(Metadata = metaDf,
           Expression = exprSum$Expression,
           ExprSummary = exprSum$Summary,
           ExprStats = exprSum$Stats,
           ParamSummary = numParams$ParamSummary,
           ParamStats = numParams$ParamStats,
           ParamTally = catParams
           )
  
  return(dataList)
}

#Calculate stuff
calculateExpression2 <- function(exprMat, metaDf, compCols) {
  sumList <- lapply(colnames(exprMat), function(gene) {
    exprGene <- exprMat[, gene, drop = FALSE]
    colnames(exprGene) <- "Value"
    
    exprGene <- cbind(metaDf[, compCols, drop = FALSE], exprGene)
    exprGene <- exprGene[, compCols] %>% tidyr::unite(., "CondGroup", sep = " ") %>%
      cbind(exprGene, .)
    
    sumDf <- dplyr::group_by_at(exprGene, "CondGroup") %>%
      dplyr::summarise(Mean = mean(Value, na.rm = TRUE),
                       SD = sd(Value, na.rm = TRUE),
                       Var = var(Value, na.rm = TRUE),
                       Count = sum(complete.cases(Value))
      )
    
    sumDf %<>% dplyr::mutate(SEM = .$SD/sqrt(.$Count))
    
    statDf <- ggpubr::compare_means(Value ~ CondGroup, exprGene)
    
    list(Expression = exprGene,
         Summary = sumDf,
         Stats = statDf
    )
  })
  names(sumList) <- colnames(exprMat)
  
  return(sumList)
}

#Export expression data
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

#Export all data
exportDataset <- function(dataList, fileName) {
  wb <- createWorkbook()
  
  addWorksheet(wb, "Metadata")
  writeData(wb, 1, dataList$Metadata)
  
  addWorksheet(wb, "Expression")
  writeData(wb, 2, dataList$Expression, rowNames = TRUE)
  
  addWorksheet(wb, "Expression Results")
  writeData(wb, 3, dataList$ExprSummary)
  
  addWorksheet(wb, "Expression Stats")
  writeData(wb, 4, dataList$ExprStats)
  
  addWorksheet(wb, "Parameter Results")
  writeData(wb, 5, dataList$ParamSummary)
  
  addWorksheet(wb, "Parameter Stats")
  writeData(wb, 6, dataList$ParamStats)
  
  addWorksheet(wb, "Catagorical Parameters")
  writeData(wb, 7, dataList$ParamTally, rowNames = TRUE)
  
  saveWorkbook(wb, file.path(getwd(), fileName))
}

#Boxplot
plotBox <- function(dataList, geneNames, compCol = NULL, selectConds = NULL, color = c("white", "tomato"), customLevels = NULL, statY = 1.02) {
  if(geneNames %v% names(dataList)) {
    stop("Gene(s) requested are absent from the dataList!")
  } else {
    dataList[["Expression"]] %<>% apply(., 2, function(gene) {
      data.frame(Comparison = dataList$Metadata[[compCol]],
                 Value = gene
                 )
    })
    dataList[["ExprStats"]] %<>% split(., .$Gene)
    
    plotList <- lapply(geneNames, function(gene) {
      exprDf <- dataList$Expression[[gene]]
      statDf <- dataList$ExprStats[[gene]]
      
      if(is.null(customLevels) == FALSE) {
        exprDf$Comparison %<>% factor(., levels = customLevels)
      }
      
      if(is.null(selectConds) == FALSE) {
        if(!selectConds %v% exprDf$Comparison) {
          stop("Desired conditions are absent from the dataList!")
        } else {
          exprDf %<>% .[.$Comparison %in% selectConds,]
          statDf %<>% .[.$group2 %in% selectConds,]
        }
      }
      
      if(is.null(color) == TRUE) {
        color <- viridis::viridis(length(unique(exprDf$Comparison)))
      } else {
        if(length(unique(exprDf$Comparison)) != length(color)) {
          stop("Wrong number of colors supplied!")
        }
      }
      
      #Get max value per condition
      craneY <- exprDf %>% split(., .$Comparison) %>% sapply(., function(comparison) max(comparison$Value, na.rm = TRUE)) %>% sort()
      craneY <- craneY*statY
      
      statY <- statY + 0.03
      
      #First make cranes between neighbours
      statDf$lvl1 <- levels(exprDf$Comparison) %>% match(statDf$group1, .)
      statDf$lvl2 <- levels(exprDf$Comparison) %>% match(statDf$group2, .)
      statDf$Neighbourhood <- abs(statDf$lvl1 - statDf$lvl2)
      
      statDf %<>% .[order(.$Neighbourhood),]
      
      cranes <- matrix(nrow = nrow(statDf)*2, ncol = 3) %>% data.frame()
      colnames(cranes) <- c("x", "y", "group")
      counter <- 0
      
      for(i in seq(1, nrow(cranes), 2)) {
        counter <- counter + 1
        cranes[i, "x"] <- statDf[counter, "group1"]
        cranes[i+1, "x"] <- statDf[counter, "group2"]
        
        if(statDf[counter, "Neighbourhood"] < 2) {
          cranes[i, "y"] <- max(craneY[c(cranes[i, "x"], cranes[i+1, "x"])])
          cranes[i+1, "y"] <- max(craneY[c(cranes[i, "x"], cranes[i+1, "x"])])
        } else {
          if(cranes[i, "x"] %in% cranes$x || cranes[i+1, "x"] %in% cranes$x) {
            cranes[i, "y"] <-  max(cranes[cranes$x %in% c(cranes[i, "x"], cranes[i+1, "x"]), "y"], na.rm = TRUE)*statY
            cranes[i+1, "y"] <- cranes[i, "y"]
          }
        }

        cranes[i, "group"] <- i
        cranes[i+1, "group"] <- i
      }
      cranes$x %<>% factor(., levels = levels(exprDf$Comparison))
      
      statDf$p.adj %<>% pcrFormatStatLabels(., "Numeric")
      statDf$x <- cranes %>% split(., .$group) %>%
        sapply(., function(group) as.double(group$x) %>% mean())
      statDf$y <- cranes %>% split(., .$group) %>%
        sapply(., function(group) unique(group$y)*1.015)
      
      limits <- c(0, ceiling(max(cranes$y)))
      
      Plot <- ggplot() +
        geom_boxplot(data = exprDf, aes(x = Comparison, y = Value, fill = Comparison), outlier.shape = NA) +
        geom_point(data = exprDf, aes(x = Comparison, y = Value), 
                   position = position_jitter(0.3), size = 1, alpha = 0.5, fill = "snow", color = "grey65", shape = 21) +
        scale_fill_manual(values = color) +
        geom_line(data = cranes, aes(x = x, y = y, group = group), size = 0.5, color = "gray40") +
        geom_text(data = statDf, aes(x = x, y = y, label = p.adj), size = 15, family = "docktrin") +
        scale_y_continuous(expand = c(0, 0), limits = limits) +
        labs(y = bquote(log[2] ~ .(gene) ~ "mRNA level")) +
        theme_classic() +
        theme(axis.line = element_line(size = 0.75, color = "grey35"), axis.ticks = element_line(size = 0.75, color = "grey35"),
              axis.title.x = element_blank(), axis.title.y = element_text(size = 80, family = "docktrin", margin = margin(r = -10)),
              axis.text.x = element_text(size = 70, family = "docktrin"), axis.text.y = element_text(size = 60, family = "docktrin"),
              legend.position = "none"
        )
      
      return(Plot)
    })
    
    return(plotList)
  }
}

#Export the graph as a tiff file
pcrSaveGraph <- function(plots, ncol = length(plots), nrow = 1, expDir, wFolder, tiffName, height = 8, width = 20, res = 600, units = "in") {
  tiffAddress <- file.path(expDir, wFolder, tiffName)
  
  cat("Creating tiff file, address:", tiffAddress, "\n", "\n")
  
  ggpubr::ggarrange(plotlist = plots, nrow = nrow, ncol = ncol)
  ggsave(tiffAddress, height = height, width = width, units = units, dpi = res)
}
