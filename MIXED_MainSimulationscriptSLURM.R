# MainSimulationScriptSLURM (following the code of the manual)
args <- commandArgs(TRUE) #SLURM command
args <- as.numeric(args)

RowOfDesign <- args[1]
Replication <- args[2]

############################# Simulation Design  #############################
factors <- c(2,4,6,8) 					            #number of latent variables
nobs <- c(200,400,800)                      #sample size
##Create the simulation design matrix (full factorial)
Design_mixed <- expand.grid(factors = factors,nobs = nobs)

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
  
  input<- cbind(seq(1,550), rep(0,550))
  mis500<- as.matrix(input)
 
  # Generate data
  # set a random number seed to be able to replicate the result exactly
  set.seed((Replication + 1000)*RowOfDesign)
  SimDat <- do.call(MyDataGeneration, Design_mixed[RowOfDesign,] )
  
#### WLS - Model without specified cross-loadings (with silent check) ####
  
  fit_WLS_withoutC <- try(cfa( model <- model_withoutC(fact), 
                               data=SimDat, std.lv=TRUE, 
                               ordered = colnames(SimDat[indexes_ord(fact)]),
                               estimator="WLS"),silent=TRUE)
  if(inherits(fit_WLS_withoutC, "try-error")) {
    mis500[Replication,1] <- 1        
  } else { 
    
    # parameter estimates
    index <- which(fit_WLS_withoutC@ParTable$free != 0)
    MyAnalysisResult_WLS_withoutC_est <- fit_WLS_withoutC@ParTable$est[index]
    MyResult_WLS_withoutC_est <- matrix(NA, 
                                        nrow = Replication, 
                                        ncol = length(ColnamesGeneratorEst("WLS", 
                                                                           "withoutC", 
                                                                           fact)))
    colnames(MyResult_WLS_withoutC_est) <- ColnamesGeneratorEst("WLS", "withoutC", fact)
    names(MyAnalysisResult_WLS_withoutC_est) <- ColnamesGeneratorEst("WLS", "withoutC", fact)
    MyResult_WLS_withoutC_est[Replication, ] <- MyAnalysisResult_WLS_withoutC_est
    
    # standard errors
    MyAnalysisResult_WLS_withoutC_err <- fit_WLS_withoutC@ParTable$se[index]
    MyResult_WLS_withoutC_err <- matrix(NA, nrow = Replication, 
                                        ncol = length(ColnamesGeneratorSE("WLS", 
                                                                          "withoutC", 
                                                                          fact)))
    colnames(MyResult_WLS_withoutC_err) <- ColnamesGeneratorSE("WLS", 
                                                               "withoutC", 
                                                               fact)
    names(MyAnalysisResult_WLS_withoutC_err) <- ColnamesGeneratorSE("WLS", "withoutC", fact)
    MyResult_WLS_withoutC_err[Replication, ] <- MyAnalysisResult_WLS_withoutC_err
    
    ### FITINDICES
    MyResult_WLS_withoutC_FI<- matrix(NA, nrow = Replication, ncol = 5)
    colnames(MyResult_WLS_withoutC_FI) <- c("chisq","df", 
                                            "pvalue", "cfi",
                                            "srmr")
    FI_WLS_withoutC <- fitMeasures(fit_WLS_withoutC, 
                                   c("chisq","df", 
                                     "pvalue", "cfi",
                                     "srmr"))
    MyResult_WLS_withoutC_FI[Replication,] <- FI_WLS_withoutC
  }
  
#### WLS - Model with specified cross-loadings (with silent check) ####
    fit_WLS_withC <- try(cfa( model <- model_withC(fact), 
                              data=SimDat, std.lv=TRUE, 
                              ordered = colnames(SimDat[indexes_ord(fact)]),
                              estimator="WLS"),silent=TRUE)
    if(inherits(fit_WLS_withC, "try-error")) {
      mis500[Replication,2] <- 1        
    } else {
      # parameter estimates
      index <- which(fit_WLS_withC@ParTable$free != 0)
      MyResult_WLS_withC_est <- matrix(NA, 
                                       nrow = Replication, 
                                       ncol = length(ColnamesGeneratorEst("WLS", 
                                                                          "withC", 
                                                                          fact)))
      colnames(MyResult_WLS_withC_est) <- ColnamesGeneratorEst("WLS", "withC", fact)
      MyAnalysisResult_WLS_withC_est <- fit_WLS_withC@ParTable$est[index]
      names(MyAnalysisResult_WLS_withC_est) <- ColnamesGeneratorEst("WLS", "withC", fact)
      MyResult_WLS_withC_est[Replication, ] <- MyAnalysisResult_WLS_withC_est
      
      # standard errors 
      MyAnalysisResult_WLS_withC_err <- fit_WLS_withC@ParTable$se[index]
      MyResult_WLS_withC_err <- matrix(NA, 
                                       nrow = Replication, 
                                       ncol = length(ColnamesGeneratorSE("WLS", 
                                                                         "withC", 
                                                                         fact)))
      colnames(MyResult_WLS_withC_err) <- ColnamesGeneratorSE("WLS", "withC", fact)
      names(MyAnalysisResult_WLS_withC_err) <- ColnamesGeneratorSE("WLS", "withC", fact)
      MyResult_WLS_withC_err[Replication, ] <- MyAnalysisResult_WLS_withC_err
      
      ### FITINDICES
      FI_WLS_withC <- fitMeasures(fit_WLS_withC, c("chisq","df", 
                                                   "pvalue", "cfi",
                                                   "srmr"))
      MyResult_WLS_withC_FI <- matrix(NA, nrow = Replication, ncol = 5)
      colnames(MyResult_WLS_withC_FI) <- c("chisq","df", 
                                           "pvalue", "cfi",
                                           "srmr")
      MyResult_WLS_withC_FI[Replication,] <- FI_WLS_withC
    }
    
    
 #### PML - Model without specified cross-loadings (with silent check) #### 
    fit_PML_withoutC <- try(cfa( model <- model_withoutC(fact), 
                                 data=SimDat, std.lv=TRUE, 
                                 ordered = colnames(SimDat[indexes_ord(fact)]),
                                 estimator="PML"),silent=TRUE)
    if(inherits(fit_PML_withoutC, "try-error")) {
      mis500[Replication,1] <- 1        
    } else { 
      
      # parameter estimates
      index <- which(fit_PML_withoutC@ParTable$free != 0)
      MyResult_PML_withoutC_est <- matrix(NA, 
                                          nrow = Replication, 
                                          ncol = length(ColnamesGeneratorEst("PML", 
                                                                             "withoutC", 
                                                                             fact)))
      colnames(MyResult_PML_withoutC_est) <- ColnamesGeneratorEst("PML", "withoutC", fact)
      MyAnalysisResult_PML_withoutC_est <- fit_PML_withoutC@ParTable$est[index]
      names(MyAnalysisResult_PML_withoutC_est) <- ColnamesGeneratorEst("PML", "withoutC", fact)
      MyResult_PML_withoutC_est[Replication, ] <- MyAnalysisResult_PML_withoutC_est
      
      # standard errors
      MyAnalysisResult_PML_withoutC_err <- fit_PML_withoutC@ParTable$se[index]
      MyResult_PML_withoutC_err <- matrix(NA, 
                                          nrow = Replication, 
                                          ncol = length(ColnamesGeneratorSE("PML", 
                                                                            "withoutC", 
                                                                            fact)))
      colnames(MyResult_PML_withoutC_err) <- ColnamesGeneratorSE("PML", "withoutC", fact)
      names(MyAnalysisResult_PML_withoutC_err) <- ColnamesGeneratorSE("PML", "withoutC", fact)
      MyResult_PML_withoutC_err[Replication, ] <- MyAnalysisResult_PML_withoutC_err
      
      ### FITINDICES
      FI_PML_withoutC <- fitMeasures(fit_PML_withoutC, 
                                     c("chisq.scaled","df.scaled", 
                                       "pvalue.scaled", "cfi.scaled",
                                       "srmr"))
      MyResult_PML_withoutC_FI<- matrix(NA, nrow = Replication, ncol = 5)
      colnames(MyResult_PML_withoutC_FI) <- c("chisq.scaled", "df.scaled", 
                                              "pvalue.scaled", "cfi.scaled",
                                              "srmr")
      MyResult_PML_withoutC_FI[Replication,] <- FI_PML_withoutC
    }
    
#### PML - Model with specified cross-loadings (with silent check) #### 
    fit_PML_withC <- try(cfa( model <- model_withC(fact), 
                              data=SimDat, std.lv=TRUE, 
                              ordered = colnames(SimDat[indexes_ord(fact)]),
                              estimator="PML"),silent=TRUE)
    if(inherits(fit_PML_withC, "try-error")) {
      mis500[Replication,2] <- 1        
    } else {
      # parameter estimates
      index <- which(fit_PML_withC@ParTable$free != 0)
      MyResult_PML_withC_est <- matrix(NA, 
                                       nrow = Replication, 
                                       ncol = length(ColnamesGeneratorEst("PML", 
                                                                          "withC", 
                                                                          fact)))
      colnames(MyResult_PML_withC_est) <- ColnamesGeneratorEst("PML", "withC", fact)
      MyAnalysisResult_PML_withC_est <- fit_PML_withC@ParTable$est[index]
      names(MyAnalysisResult_PML_withC_est) <- ColnamesGeneratorEst("PML", "withC", fact)
      MyResult_PML_withC_est[Replication, ] <- MyAnalysisResult_PML_withC_est
      
      # standard errors 
      MyAnalysisResult_PML_withC_err <- fit_PML_withC@ParTable$se[index]
      MyResult_PML_withC_err <- matrix(NA, 
                                       nrow = Replication, 
                                       ncol = length(ColnamesGeneratorSE("PML", 
                                                                         "withC", 
                                                                         fact)))
      colnames(MyResult_PML_withC_err) <- ColnamesGeneratorSE("PML", "withC", fact)
      names(MyAnalysisResult_PML_withC_err) <- ColnamesGeneratorSE("PML", "withC", fact)
      MyResult_PML_withC_err[Replication, ] <- MyAnalysisResult_PML_withC_err
      
      ### FITINDICES
      FI_PML_withC <- fitMeasures(fit_PML_withC, c("chisq.scaled","df.scaled", 
                                                   "pvalue.scaled", "cfi.scaled",
                                                   "srmr"))
      
      MyResult_PML_withC_FI <- matrix(NA, nrow = Replication, ncol = 5)
      colnames(MyResult_PML_withC_FI) <- c("chisq.scaled", "df.scaled", 
                                           "pvalue.scaled", "cfi.scaled",
                                           "srmr")
      MyResult_PML_withC_FI[Replication,] <- FI_PML_withC
    }
  
################################ Simulation all cells  ###############################
  setwd("/exports/fsw/mmcstorm/Analysis/mixed")

  # Write output of one cell of the design
  # Save results
  
  write.csv(MyResult_WLS_withoutC_est, 
            file = paste("WLS_withoutC_est", "Row", RowOfDesign, ".csv" , sep =""))
  
  write.csv(MyResult_WLS_withoutC_err,
            file = paste("WLS_withoutC_err", "Row", RowOfDesign, ".csv" , sep =""))
  
  write.csv(MyResult_WLS_withoutC_FI,
            file =paste("WLS_FI_withoutC", "Row", RowOfDesign,".csv" , sep =""))
  
  write.csv(MyResult_WLS_withC_est,
            file = paste("WLS_withC_est", "Row", RowOfDesign, ".csv" , sep =""))
  
  write.csv(MyResult_WLS_withC_err,
            file = paste("WLS_withC_err", "Row", RowOfDesign,".csv" , sep =""))
  
  write.csv(MyResult_WLS_withC_FI,
            file =paste("WLS_FI_withC", "Row", RowOfDesign, ".csv" , sep =""))
  
  write.csv(MyResult_PML_withoutC_est,
            file = paste("PML_withoutC_est", "Row", RowOfDesign ,".csv" , sep =""))
  
  write.csv(MyResult_PML_withoutC_err,
            file = paste("PML_withoutC_err", "Row", RowOfDesign ,".csv" , sep =""))
  
  write.csv(MyResult_PML_withoutC_FI,
            file =paste("PML_FI_withoutC", "Row", RowOfDesign,".csv" , sep =""))
  
  write.csv(MyResult_PML_withC_est,
            file = paste("PML_withC_est", "Row", RowOfDesign,".csv" , sep =""))
  
  write.csv(MyResult_PML_withC_err,
            file = paste("PML_withC_err", "Row", RowOfDesign,".csv" , sep =""))
  
  write.csv(MyResult_PML_withC_FI,
            file =paste("PML_FI_withC", "Row", RowOfDesign, ".csv" , sep =""))
  
  #save the time to run the analyses of K data sets in one cell of the design.
  time <- proc.time() - tmp
  save(time, file =paste("Time", "Row", RowOfDesign, ".csv" , sep =""))
  
  # see whether all replications are ok
  write.csv(mis500, 
            file = paste("mis500", "Row", RowOfDesign, ".csv" , sep =""))
  
  setwd("/exports/fsw/mmcstorm/SimData/mixed")
  # save data
  write.csv(SimDat, 
            file = paste("Simulated_Data", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))