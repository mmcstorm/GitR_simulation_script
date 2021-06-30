#################### LET OP: Dit script bevat andere mis500 variabelen!! Zie meeting Mariska 14-06-2021

# MainSimulationScriptSLURM (following the code of the manual)
args <- commandArgs(TRUE) #SLURM command
args <- as.numeric(args)

RowOfDesign <- args[1]
Replication <- args[2]

############################# Simulation Design  #############################
factors <- c(2,4,6,8)					             #number of latent variables
nobs <- c(200,400,800)                      #sample size
##Create the simulation design matrix (full factorial)
Design_mixed <- expand.grid(factors = factors, nobs = nobs)

#load packages
library(lavaan)
library(usethis)

# load functions
source("MIXED_all_functions_script.R")

################################ Simulation start (1 cell) ##########################
# initialize values
tmp <- proc.time()
nvarp <- 6
fact <- Design_mixed[RowOfDesign,1]
nvar <- nvarp*fact

input<- cbind(1, 0)
mis_W_withoutC <- as.matrix(input)
mis_W_withC <- as.matrix(input)
mis_P_withoutC <- as.matrix(input)
mis_P_withC <- as.matrix(input)

# Generate data
# set a random number seed to be able to replicate the result exactly
set.seed((Replication + 1000)*RowOfDesign)
SimDat <- do.call(MyDataGeneration, Design_mixed[RowOfDesign,] )

#### WLS - Model without specified cross-loadings (with silent check) ####

fit_WLS_withoutC <- try(cfa( model <- model_withoutC(fact), 
                             data=SimDat, std.lv=TRUE, 
                             ordered = colnames(SimDat[indexes_ord(fact)]),
                             estimator="WLSMV"),silent=TRUE)
if(inherits(fit_WLS_withoutC, "try-error")) {
  mis_W_withoutC[1,2] <- 1
  MyResult_WLS_withoutC_est <- matrix(NA,
                                      nrow = 1,
                                      ncol = length(ColnamesGeneratorEst("WLS",
                                                                         "withoutC",
                                                                         fact)))
  MyResult_WLS_withoutC_err <- matrix(NA, 
                                      nrow = 1, 
                                      ncol = length(ColnamesGeneratorSE("WLS", 
                                                                        "withoutC", 
                                                                        fact)))
  MyResult_WLS_withoutC_FI<- matrix(NA, nrow = 1, ncol = 5)
} else { 
  
  # parameter estimates
  index <- which(fit_WLS_withoutC@ParTable$free != 0)
  MyAnalysisResult_WLS_withoutC_est <- fit_WLS_withoutC@ParTable$est[index]
  
  MyResult_WLS_withoutC_est <- matrix(MyAnalysisResult_WLS_withoutC_est,
                                      nrow = 1,
                                      ncol = length(ColnamesGeneratorEst("WLS",
                                                                         "withoutC",
                                                                         fact)))
  colnames(MyResult_WLS_withoutC_est) <- ColnamesGeneratorEst("WLS", "withoutC", fact)
  # standard errors
  MyAnalysisResult_WLS_withoutC_err <- fit_WLS_withoutC@ParTable$se[index]
  MyResult_WLS_withoutC_err <- matrix(MyAnalysisResult_WLS_withoutC_err, 
                                      nrow = 1, 
                                      ncol = length(ColnamesGeneratorSE("WLS", 
                                                                        "withoutC", 
                                                                        fact)))
  colnames(MyResult_WLS_withoutC_err) <- ColnamesGeneratorSE("WLS", 
                                                             "withoutC", 
                                                             fact)
  
  ### FITINDICES
  FI_WLS_withoutC <- fitMeasures(fit_WLS_withoutC, 
                                 c("chisq","df", 
                                   "pvalue", "cfi",
                                   "srmr"))
  MyResult_WLS_withoutC_FI<- matrix(FI_WLS_withoutC, nrow = 1, ncol = 5)
  colnames(MyResult_WLS_withoutC_FI) <- c("chisq","df", 
                                          "pvalue", "cfi",
                                          "srmr")
}

#### WLS - Model with specified cross-loadings (with silent check) ####
fit_WLS_withC <- try(cfa( model <- model_withC(fact), 
                          data=SimDat, std.lv=TRUE, 
                          ordered = colnames(SimDat[indexes_ord(fact)]),
                          estimator="WLSMV"),silent=TRUE)
if(inherits(fit_WLS_withC, "try-error")) {
  mis_W_withC[1,2] <- 1
  MyResult_WLS_withC_est <- matrix(NA,
                                      nrow = 1,
                                      ncol = length(ColnamesGeneratorEst("WLS",
                                                                         "withC",
                                                                         fact)))
  MyResult_WLS_withC_err <- matrix(NA, 
                                      nrow = 1, 
                                      ncol = length(ColnamesGeneratorSE("WLS", 
                                                                        "withC", 
                                                                        fact)))
  MyResult_WLS_withC_FI<- matrix(NA, nrow = 1, ncol = 5)
} else {
  # parameter estimates
  index <- which(fit_WLS_withC@ParTable$free != 0)
  MyAnalysisResult_WLS_withC_est <- fit_WLS_withC@ParTable$est[index]
  MyResult_WLS_withC_est <- matrix(MyAnalysisResult_WLS_withC_est, 
                                   nrow = 1, 
                                   ncol = length(ColnamesGeneratorEst("WLS", 
                                                                      "withC", 
                                                                      fact)))
  colnames(MyResult_WLS_withC_est) <- ColnamesGeneratorEst("WLS", "withC", fact)
  
  # standard errors 
  MyAnalysisResult_WLS_withC_err <- fit_WLS_withC@ParTable$se[index]
  MyResult_WLS_withC_err <- matrix(MyAnalysisResult_WLS_withC_err, 
                                   nrow = 1, 
                                   ncol = length(ColnamesGeneratorSE("WLS", 
                                                                     "withC", 
                                                                     fact)))
  colnames(MyResult_WLS_withC_err) <- ColnamesGeneratorSE("WLS", "withC", fact)
  
  ### FITINDICES
  FI_WLS_withC <- fitMeasures(fit_WLS_withC, c("chisq","df", 
                                               "pvalue", "cfi",
                                               "srmr"))
  MyResult_WLS_withC_FI <- matrix(FI_WLS_withC, nrow = 1, ncol = 5)
  colnames(MyResult_WLS_withC_FI) <- c("chisq","df", 
                                       "pvalue", "cfi",
                                       "srmr")
}

#### PML - Model without specified cross-loadings (with silent check) #### 
fit_PML_withoutC <- try(cfa( model <- model_withoutC(fact), 
                             data=SimDat, std.lv=TRUE, 
                             ordered = colnames(SimDat[indexes_ord(fact)]),
                             estimator="PML"),silent=TRUE)
if(inherits(fit_PML_withoutC, "try-error")) {
  mis_P_withoutC[1,2] <- 1
  MyResult_PML_withoutC_est <- matrix(NA, 
                                      nrow = 1, 
                                      ncol = length(ColnamesGeneratorEst("PML", 
                                                                         "withoutC", 
                                                                         fact)))
  MyResult_PML_withoutC_err <- matrix(NA, 
                                      nrow = 1, 
                                      ncol = length(ColnamesGeneratorSE("PML", 
                                                                        "withoutC", 
                                                                        fact)))
  MyResult_PML_withoutC_FI<- matrix(NA, nrow = 1, ncol = 5)
} else { 
  
  # parameter estimates
  index <- which(fit_PML_withoutC@ParTable$free != 0)
  MyAnalysisResult_PML_withoutC_est <- fit_PML_withoutC@ParTable$est[index]
  MyResult_PML_withoutC_est <- matrix(MyAnalysisResult_PML_withoutC_est, 
                                      nrow = 1, 
                                      ncol = length(ColnamesGeneratorEst("PML", 
                                                                         "withoutC", 
                                                                         fact)))
  colnames(MyResult_PML_withoutC_est) <- ColnamesGeneratorEst("PML", "withoutC", fact)
  
  # standard errors
  MyAnalysisResult_PML_withoutC_err <- fit_PML_withoutC@ParTable$se[index]
  MyResult_PML_withoutC_err <- matrix(MyAnalysisResult_PML_withoutC_err, 
                                      nrow = 1, 
                                      ncol = length(ColnamesGeneratorSE("PML", 
                                                                        "withoutC", 
                                                                        fact)))
  colnames(MyResult_PML_withoutC_err) <- ColnamesGeneratorSE("PML", "withoutC", fact)
  
  ### FITINDICES
  FI_PML_withoutC <- fitMeasures(fit_PML_withoutC, 
                                 c("chisq.scaled","df.scaled", 
                                   "pvalue.scaled", "cfi.scaled",
                                   "srmr"))
  MyResult_PML_withoutC_FI<- matrix(FI_PML_withoutC, nrow = 1, ncol = 5)
  colnames(MyResult_PML_withoutC_FI) <- c("chisq.scaled", "df.scaled", 
                                          "pvalue.scaled", "cfi.scaled",
                                          "srmr")
}

#### PML - Model with specified cross-loadings (with silent check) #### 
fit_PML_withC <- try(cfa( model <- model_withC(fact), 
                          data=SimDat, std.lv=TRUE, 
                          ordered = colnames(SimDat[indexes_ord(fact)]),
                          estimator="PML"),silent=TRUE)
if(inherits(fit_PML_withC, "try-error")) {
  mis_P_withC[1,2] <- 1
  MyResult_PML_withC_est <- matrix(MyAnalysisResult_PML_withC_est, 
                                   nrow = 1, 
                                   ncol = length(ColnamesGeneratorEst("PML", 
                                                                      "withC", 
                                                                      fact)))
  MyResult_PML_withC_err <- matrix(MyAnalysisResult_PML_withC_err, 
                                   nrow = 1, 
                                   ncol = length(ColnamesGeneratorSE("PML", 
                                                                     "withC", 
                                                                     fact)))
  MyResult_PML_withC_FI <- matrix(FI_PML_withC, nrow = 1, ncol = 5)
} else {
  # parameter estimates
  index <- which(fit_PML_withC@ParTable$free != 0)
  MyAnalysisResult_PML_withC_est <- fit_PML_withC@ParTable$est[index]
  MyResult_PML_withC_est <- matrix(MyAnalysisResult_PML_withC_est, 
                                   nrow = 1, 
                                   ncol = length(ColnamesGeneratorEst("PML", 
                                                                      "withC", 
                                                                      fact)))
  colnames(MyResult_PML_withC_est) <- ColnamesGeneratorEst("PML", "withC", fact)
  
  # standard errors 
  MyAnalysisResult_PML_withC_err <- fit_PML_withC@ParTable$se[index]
  MyResult_PML_withC_err <- matrix(MyAnalysisResult_PML_withC_err, 
                                   nrow = 1, 
                                   ncol = length(ColnamesGeneratorSE("PML", 
                                                                     "withC", 
                                                                     fact)))
  colnames(MyResult_PML_withC_err) <- ColnamesGeneratorSE("PML", "withC", fact)
  
  ### FITINDICES
  FI_PML_withC <- fitMeasures(fit_PML_withC, c("chisq.scaled","df.scaled", 
                                               "pvalue.scaled", "cfi.scaled",
                                               "srmr"))
  MyResult_PML_withC_FI <- matrix(FI_PML_withC, nrow = 1, ncol = 5)
  colnames(MyResult_PML_withC_FI) <- c("chisq.scaled", "df.scaled", 
                                       "pvalue.scaled", "cfi.scaled",
                                       "srmr")
}

################################ Simulation all cells  ###############################
setwd("/exports/fsw/mmcstorm/Analysis/finalmixedrep151to175")

# Write output of one cell of the design
# Save results

write.csv(MyResult_WLS_withoutC_est, 
          file = paste("WLS_withoutC_est", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_WLS_withoutC_err,
          file = paste("WLS_withoutC_err", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_WLS_withoutC_FI,
          file =paste("WLS_FI_withoutC", "Row", RowOfDesign,"Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_WLS_withC_est,
          file = paste("WLS_withC_est", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_WLS_withC_err,
          file = paste("WLS_withC_err", "Row", RowOfDesign,"Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_WLS_withC_FI,
          file =paste("WLS_FI_withC", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_PML_withoutC_est,
          file = paste("PML_withoutC_est", "Row", RowOfDesign ,"Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_PML_withoutC_err,
          file = paste("PML_withoutC_err", "Row", RowOfDesign ,"Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_PML_withoutC_FI,
          file =paste("PML_FI_withoutC", "Row", RowOfDesign,"Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_PML_withC_est,
          file = paste("PML_withC_est", "Row", RowOfDesign,"Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_PML_withC_err,
          file = paste("PML_withC_err", "Row", RowOfDesign,"Rep", Replication, ".csv" , sep =""))

write.csv(MyResult_PML_withC_FI,
          file =paste("PML_FI_withC", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))

#save the time to run the analyses of K data sets in one cell of the design.
#time <- proc.time() - tmp
#save(time, file =paste("Time", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))

# see whether all replications are ok
setwd("/exports/fsw/mmcstorm/Simdata/mixed/silent_check")
write.csv(mis_W_withoutC, 
          file = paste("mis_W_withoutC", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))

write.csv(mis_W_withC, 
          file = paste("mis_W_withC", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))

write.csv(mis_P_withoutC, 
          file = paste("mis_P_withoutC", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))

write.csv(mis_P_withC, 
          file = paste("mis_P_withC", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))

setwd("/exports/fsw/mmcstorm/Simdata/mixed")
# save data
write.csv(SimDat, 
          file = paste("Simulated_Data", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))