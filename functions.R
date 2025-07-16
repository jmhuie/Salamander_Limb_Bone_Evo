
# Load GMM ----------------------------------------------------------------
library(SlicerMorphR)
library(jsonlite)
read_SlicerLand <- function(path = NULL) {
  files <- list.files(path)
  if (length(grep(".txt", files)) >= 1) {
    files <- files[-c(grep(".txt", files))]
  }
  if (length(grep("json",files[1])) > 0) {
    npoint = nrow(SlicerMorphR::read.markups.json(paste0(path,files[1])))
  }
  if (length(grep("fcsv",files[1])) > 0) {
    npoint = nrow(SlicerMorphR::read.markups.fcsv(paste0(path,files[1])))
  }
  land = array(dim = c(npoint, 3, length(files)))
  for (i in 1:length(files)) {
    if (length(grep("json",files[i])) > 0) {
      land[,,i] <- SlicerMorphR::read.markups.json(paste0(path,files[i]))
    } 
    if (length(grep("fcsv",files[i])) > 0) {
      land[,,i] <- SlicerMorphR::read.markups.fcsv(paste0(path,files[i]))
    } 
  }
  dimnames(land)[[1]] <- 1:npoint
  dimnames(land)[[2]] <- c("x","y","z")
  dimnames(land)[[3]] <- stringr::str_extract(files, "[^_]*_[^_]*")
  return(land)
}

# batch corHMM ------------------------------------------------------------

library(corHMM)
library(MASS)
RunModels <- function(phy, data, dual, name){
  
  # if polymorphic and ambiguous == FALSE, then species get trapped in one of two states without dual transitions so we add them in.
  if(dual == TRUE){
    ER <- getStateMat4Dat(data, "ER", dual = T)$rate.mat
    SYM <- getStateMat4Dat(data, "SYM", dual = T)$rate.mat
    ARD <- getStateMat4Dat(data, "ARD", dual = T)$rate.mat    
    RC <- getRateCatMat(2)
    ER2 <- getFullMat(list(ER, ER), RC)
    SYM2 <- getFullMat(list(SYM, SYM), RC)
    ARD2 <- getFullMat(list(ARD, ARD), RC)
  }else{
    ER <- getStateMat4Dat(data, "ER")$rate.mat
    SYM <- getStateMat4Dat(data, "SYM")$rate.mat
    ARD <- getStateMat4Dat(data, "ARD")$rate.mat
    RC <- getRateCatMat(2)
    ER2 <- getFullMat(list(ER, ER), RC)
    SYM2 <- getFullMat(list(SYM, SYM), RC)
    ARD2 <- getFullMat(list(ARD, ARD), RC)
  }
  
  RateCats <- c(1,1,1,2,2,2)
  Models <- list(ER, SYM, ARD, ER2, SYM2, ARD2)
  CorRes_i <- vector("list", length(Models))
  # fit the models
  for(i in 1:length(Models)){
    CorRes_i[[i]] <- corHMM(phy, data, RateCats[i], rate.mat = Models[[i]], nstarts = 10, n.cores = 11)
  }
  save(CorRes_i, file = name)
}

# saves the results table and returns it
getResultsTable <- function(file){
  load(file)
  # the results table
  obj <- CorRes_i
  ResTable <- matrix(0, length(obj), 5)
  rownames(ResTable) <- rep(c("ER", "SYM", "ARD", "ER/ER", "SYM/SYM", "ARD/ARD"), length(obj)/6)
  colnames(ResTable) <- c("k.rate", "AICc", "AICcWt", "MeanRate", "ASR")
  count <- 1
  for(i in 1:length(obj)){
    CorRes_i <- obj[[i]]
    Model <- CorRes_i$solution
    AICc <- round(CorRes_i$AICc, 2)
    AICwt <- 0
    MeanRate <- round(mean(CorRes_i$solution, na.rm = TRUE),2)
    Model[is.na(Model)] <- 0
    ASR <-  round(CorRes_i$states[1,][which.max(CorRes_i$states[1,])], 2)
    ASR <-  paste(names(ASR), paste(ASR*100, "%", sep = ""))
    k.rate <- max(CorRes_i$index.mat, na.rm = TRUE)
    diag(Model) <- -rowSums(Model)
    Eq <- c(Null(Model))/sum(Null(Model))
    EntUncond <- sum(Eq*-log2(Eq))
    EntCond <- mean(rowSums(CorRes_i$states * -log2(CorRes_i$states)))
    MutInfo <- round(EntUncond - EntCond, 2)
    #PropInfo <- round(MutInfo/EntUncond*100, 2)
    ResTable[count,] <- c(k.rate, AICc, AICwt, MeanRate, ASR)
    count <- count + 1
  }
  AICcs <- as.numeric(ResTable[,2])
  ResTable[,3] <- round(exp(-0.5 * AICcs - min(AICcs))/sum(exp(-0.5 * AICcs - min(AICcs))),2 )
  ResTable[ResTable[,3] < 0.01,3] <- "<0.01"
  ResTable[ResTable[,4] < 0.01,4] <- "<0.01"
  
  table.name <- paste("ResTable_", gsub(".Rsave", ".csv", file), sep = "")
  #write.csv(ResTable, file = table.name)
  return(ResTable)
}

# batch hOUwie ------------------------------------------------------------

c3func <- function(data, phy, models, nSim = nSim, path = NULL) {
  if (is.null(path)) {
    path <- paste0("hOUwie_", colnames(data)[length(colnames(data))], "_", nSim)
  }
  dir.create(path, showWarnings = FALSE) # Create directory if it doesn't exist
  hOUwie_dataprep <- data
  
  for (i in 1:length(models)) {
    model <- names(models)[i]
    if (length(grep(paste0(model, ".RDS"), list.files(path))) == 0) {
      print(paste0("Starting: ", model))
      
      # Retry logic for `hOUwie`
      success <- FALSE
      while (!success) {
        tryCatch({
          if (length(grep("cd", model)) == 1) {
            fit <- hOUwie(phy, hOUwie_dataprep, rate.cat = 1, discrete_model = "SYM", 
                          continuous_model = models[[i]], nSim = nSim)
          } else {
            fit <- hOUwie(phy, hOUwie_dataprep, rate.cat = 2, null.model = TRUE, 
                          discrete_model = "SYM", continuous_model = models[[i]], 
                          nSim = nSim, sample_nodes = TRUE, adaptive_sampling = TRUE)
          }
          saveRDS(fit, file = paste0(path, "/" ,model, ".RDS"))
          success <- TRUE # Exit loop if no error occurs
          print(paste0(model, " completed"))
        }, error = function(e) {
          print(paste0("Error in model: ", model, " - Retrying..."))
        })
      }
    } else {
      print(paste0(model, " already completed"))
    }
  }
}
