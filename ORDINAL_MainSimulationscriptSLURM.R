# MainSimulationScriptSLURM (following the code of the manual)
args <- commandArgs(TRUE) #SLURM command
args <- as.numeric(args)

RowOfDesign <- args[1]
Replication <- args[2]
############################# Simulation Design  #############################
factors <- c(2,4,6,8) 					            #number of latent variables
nobs <- c(200,400,800)                      #sample size
ncat <- c(2,4)                              #number of categories

##Create the simulation design matrix (full factorial)

Design <- expand.grid(factors = factors, nobs = nobs, ncat = ncat)

#load packages
library(lavaan)
library(usethis)

# load prepared functions
source("ORDINAL_all_functions_script.R")

################################ Simulation start (1 cell) ##########################
  ### Fit indices
  ### Parameter estimates
  
  #create a loop over the replications k = 1 to K:
  tmp <- proc.time()
  nvarp <- 6
  fact <- Design[RowOfDesign,1]
  nvar <- nvarp*fact

  input<- cbind(seq(1,550), rep(0,550))
  mis500<- as.matrix(input)
  
    # Generate data
    # set a random number seed to be able to replicate the result exactly
    set.seed((Replication + 1000)*RowOfDesign)
    SimDat <- do.call(MyDataGeneration, Design[RowOfDesign,] )

  ####### WLS model without specified cross-loadings ########
    fit1_W <- try(cfa( model <- model_withoutC(fact), 
                       data=SimDat, std.lv=TRUE, 
                       ordered=c(colnames(SimDat)),
                       estimator="WLSMV"),silent=TRUE)
    if(inherits(fit1_W, "try-error")) {
      mis500[Replication,1] <- 1        
    } else { 
      # parameter estimates
      index <- which(fit1_W@ParTable$free != 0)
      MyAnalysisResult_WLS1est <- fit1_W@ParTable$est[index]
      MyResult_WLS_withoutC_est <- matrix(MyAnalysisResult_WLS1est, 
                                    nrow = 1, 
                                    ncol = length(ColnamesGeneratorEst("WLS",
                                                                       "withoutC", 
                                                                       fact,
                                                                       Design[RowOfDesign,3])))
      colnames(MyResult_WLS_withoutC_est) <- ColnamesGeneratorEst("WLS","withoutC",
                                                            fact,Design[RowOfDesign,3])
      
      # standard errors 
      MyAnalysisResult_WLS1err <- fit1_W@ParTable$se[index]
      MyResult_WLS_withoutC_err <- matrix(MyAnalysisResult_WLS1err, 
                                    nrow = 1, 
                                    ncol = length(ColnamesGeneratorSE("WLS", 
                                                                      "withoutC", 
                                                                      fact, 
                                                                      Design[RowOfDesign,3])))
      colnames(MyResult_WLS_withoutC_err) <- ColnamesGeneratorSE("WLS","withoutC",
                                                           fact,Design[RowOfDesign,3])
      
      ### FITINDICES
      FI_WLS_withoutC <- fitMeasures(fit1_W, c("chisq.scaled","df.scaled", 
                                         "pvalue.scaled", "cfi.scaled", "srmr"))
      MyResult_WLS_withoutC_FI <- matrix(FI_WLS_withoutC, nrow = 1, ncol = 5)
      colnames(MyResult_WLS_withoutC_FI) <- c("chisq.scaled", "df.scaled", 
                                        "pvalue.scaled", "cfi.scaled",
                                        "srmr")
    }
    
    ####### WLS model with specified cross-loadings ########
    fit2_W <- try(cfa( model <- model_withC(fact), 
                       data=SimDat, std.lv=TRUE, 
                       ordered=c(colnames(SimDat)),
                       estimator="WLSMV"),silent=TRUE)
    if(inherits(fit2_W, "try-error")) {
      mis500[Replication,1] <- 1        
    } else { 
      # parameter estimates
      index <- which(fit2_W@ParTable$free != 0)
      MyAnalysisResult_WLS2est <- fit2_W@ParTable$est[index]
      MyResult_WLS_withC_est <- matrix(MyAnalysisResult_WLS2est, 
                                    nrow = 1, 
                                    ncol = length(ColnamesGeneratorEst("WLS",
                                                                       "withC", 
                                                                       fact,
                                                                       Design[RowOfDesign,3])))
      colnames(MyResult_WLS_withC_est) <- ColnamesGeneratorEst("WLS","withC",
                                                            fact,Design[RowOfDesign,3])
      
      
      # standard errors 
      MyAnalysisResult_WLS2err <- fit2_W@ParTable$se[index]
      MyResult_WLS_withC_err <- matrix(MyAnalysisResult_WLS2err, 
                                    nrow = 1, 
                                    ncol = length(ColnamesGeneratorSE("WLS", 
                                                                      "withC", 
                                                                      fact, 
                                                                      Design[RowOfDesign,3])))
      colnames(MyResult_WLS_withC_err) <- ColnamesGeneratorSE("WLS","withC",
                                                           fact,Design[RowOfDesign,3])
      
      ### FITINDICES
      FI_WLS_withC <- fitMeasures(fit2_W, c("chisq.scaled","df.scaled", 
                                         "pvalue.scaled", "cfi.scaled", "srmr"))
      MyResult_WLS_withC_FI <- matrix(FI_WLS_withC, nrow = 1, ncol = 5)
      colnames(MyResult_WLS_withC_FI) <- c("chisq.scaled", "df.scaled", 
                                        "pvalue.scaled", "cfi.scaled",
                                        "srmr")
    }
    
    ###### PML model without specified cross-loadings ###### 
    fit1_P <- try(cfa( model <- model_withoutC(fact), 
                       data=SimDat, std.lv=TRUE, 
                       ordered=c(colnames(SimDat)),
                       estimator="PML"),silent=TRUE)
    if(inherits(fit1_P, "try-error")) {
      mis500[Replication,1] <- 1        
    } else { 
      # parameter estimates
      index <- which(fit1_P@ParTable$free != 0)
      MyAnalysisResult_PML1est <- fit1_P@ParTable$est[index]
      MyResult_PML_withoutC_est <- matrix(MyAnalysisResult_PML1est, 
                                    nrow = 1, 
                                    ncol = length(ColnamesGeneratorEst("PML",
                                                                       "withoutC", 
                                                                       fact,
                                                                       Design[RowOfDesign,3])))
      colnames(MyResult_PML_withoutC_est) <- ColnamesGeneratorEst("PML","withoutC",
                                                            fact,Design[RowOfDesign,3])
           
      # standard errors 
      MyAnalysisResult_PML1err <- fit1_P@ParTable$se[index]
      MyResult_PML_withoutC_err <- matrix(MyAnalysisResult_PML1err, 
                                    nrow = 1, 
                                    ncol = length(ColnamesGeneratorSE("PML", 
                                                                      "withoutC", 
                                                                      fact, 
                                                                      Design[RowOfDesign,3])))
      colnames(MyResult_PML_withoutC_err) <- ColnamesGeneratorSE("PML","withoutC",
                                                           fact,Design[RowOfDesign,3])
      
      ### FITINDICES
      FI_PML_withoutC <- fitMeasures(fit1_P, c("chisq.scaled","df.scaled", 
                                         "pvalue.scaled", "cfi.scaled", "srmr"))
      MyResult_PML_withoutC_FI <- matrix(FI_PML_withoutC, nrow = 1, ncol = 5)
      colnames(MyResult_PML_withoutC_FI) <- c("chisq.scaled", "df.scaled", 
                                        "pvalue.scaled", "cfi.scaled",
                                        "srmr")
    }
    
    
    ###### PML model with specified cross-loadings #####
    fit2_P <- try(cfa( model <- model_withC(fact), 
                       data=SimDat, std.lv=TRUE, 
                       ordered=c(colnames(SimDat)),
                       estimator="PML"),silent=TRUE)
    if(inherits(fit2_P, "try-error")) { 
      mis500[Replication,1] <- 1        
    } else { 
      # parameter estimates
      index <- which(fit2_P@ParTable$free != 0)
      MyAnalysisResult_PML2est <- fit2_P@ParTable$est[index]
      MyResult_PML_withC_est <- matrix(MyAnalysisResult_PML2est, 
                                       nrow = 1, 
                                       ncol = length(ColnamesGeneratorEst("PML",
                                                                          "withC", 
                                                                          fact,
                                                                          Design[RowOfDesign,3])))
      colnames(MyResult_PML_withC_est) <- ColnamesGeneratorEst("PML","withC",fact,Design[RowOfDesign,3])
      
      # standard errors
      MyAnalysisResult_PML2err <- fit2_P@ParTable$se[index]
      MyResult_PML_withC_err <- matrix(MyAnalysisResult_PML2err, 
                                    nrow = 1, 
                                    ncol = length(ColnamesGeneratorSE("PML", 
                                                                      "withC", 
                                                                      fact, 
                                                                      Design[RowOfDesign,3])))
      colnames(MyResult_PML_withC_err) <- ColnamesGeneratorSE("PML","withC",
                                                           fact,Design[RowOfDesign,3])
      
      ### FITINDICES
      FI_PML_withC <- fitMeasures(fit2_P, c("chisq.scaled","df.scaled", 
                                         "pvalue.scaled", "cfi.scaled", "srmr"))
      MyResult_PML_withC_FI <- matrix(FI_PML_withC, nrow = 1, ncol = 5)
      colnames(MyResult_PML_withC_FI) <- c("chisq.scaled", "df.scaled", 
                                        "pvalue.scaled", "cfi.scaled",
                                        "srmr")
    }
    
 ################################ Simulation all cells  ###############################
setwd("/exports/fsw/mmcstorm/Analysis/R20ordinal1to800")
  # save the results
    write.csv(MyResult_WLS_withoutC_est, 
              file = paste("WLS_withoutC_est", "Row", RowOfDesign, "Rep", Replication,".csv" , sep =""))
    
    write.csv(MyResult_WLS_withoutC_err,
              file = paste("WLS_withoutC_err", "Row", RowOfDesign,"Rep", Replication,".csv" , sep =""))
    
    write.csv(MyResult_WLS_withoutC_FI,
              file =paste("WLS_FI_withoutC", "Row", RowOfDesign,"Rep", Replication, ".csv" , sep =""))
    
    write.csv(MyResult_WLS_withC_est,
              file = paste("WLS_withC_est", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))
    
    write.csv(MyResult_WLS_withC_err,
              file = paste("WLS_withC_err", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))
    
    write.csv(MyResult_WLS_withC_FI,
              file =paste("WLS_FI_withC", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))
    
    write.csv(MyResult_PML_withoutC_est,
              file = paste("PML_withoutC_est", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))
    
    write.csv(MyResult_PML_withoutC_err,
              file = paste("PML_withoutC_err", "Row", RowOfDesign,"Rep", Replication,".csv" , sep =""))
    
    write.csv(MyResult_PML_withoutC_FI,
              file =paste("PML_FI_withoutC", "Row", RowOfDesign,"Rep", Replication, ".csv" , sep =""))
    
    write.csv(MyResult_PML_withC_est,
              file = paste("PML_withC_est", "Row", RowOfDesign,"Rep", Replication, ".csv" , sep =""))
    
    write.csv(MyResult_PML_withC_err,
              file = paste("PML_withC_err", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))
    
    write.csv(MyResult_PML_withC_FI,
              file = paste("PML_FI_withC", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))
    
    #save the time to run the analyses of K data sets in one cell of the design.
    time <- proc.time() - tmp
    save(time, file =paste("Time", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))
    
    # see whether all replications are ok
    setwd("/exports/fsw/mmcstorm/Simdata/silent_check")
    write.csv(mis500, 
              file = paste("mis500", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))
    
    # save data
    write.csv(SimDat, 
              file = paste("Simulated_Data", "Row", RowOfDesign, "Rep", Replication, ".csv" , sep =""))
  

setwd("/exports/fsw/mmcstorm/Simdata/ordinal")
# create folder data
save(SimDat, file =paste("Data", "Row", RowOfDesign, "Rep", Replication ,".Rdata" , sep =""))