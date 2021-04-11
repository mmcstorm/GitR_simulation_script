# MainSimulationScript for running on my computer (following the code of the manual)

# load packages
library(lavaan)
library(usethis)

#args <- commandArgs(TRUE) # SLURM language
#args <- as.numeric(args)

#RowOfDesign <- args[1]
#Replication <- args[2]

RowOfDesign <- 1 #24
Replication <- 1

############################# Simulation Design  #############################
factors <- c(2,4,6,8) 					            #number of latent variables
nobs <- c(200,400,800)                      #sample size
ncat <- c(2,4)                              #number of categories

##Create the simulation design matrix (full factorial)
# Design is a data.frame with all possible combinations of the factor levels
# Each row of the design matrix represents a cell of your simulation design
# Design <- expand.grid(factors = factors, nobs = nobs, ncat = ncat)

Design <- expand.grid(factors = factors, nobs = nobs, ncat = ncat)

# prepared functions
source("MyDataGeneration.R")
source("Method_New.R")
source("Method_Old.R")
source("ColnamesGenerators.R")

################################ Simulation start (1 cell) ##########################
MySimulationCell<- function(Design = Design, RowOfDesign, K){
  
  # initialize values
  nvarp <- 6
  fact <- Design[RowOfDesign,1]
  nvar <- nvarp*fact
  
  # Create matrices to store the results:
  ## WLS_CS
  MyResult_WLS_CS_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("WLS",
  "CS", fact,Design[RowOfDesign,3])))
  colnames(MyResult_WLS_CS_est) <- ColnamesGeneratorEst("WLS","CS",fact,Design[RowOfDesign,3])
  
  MyResult_WLS_CS_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("WLS", "CS", fact, Design[RowOfDesign,3])))
  colnames(MyResult_WLS_CS_err) <- ColnamesGeneratorSE("WLS","CS",fact,Design[RowOfDesign,3])
  
  MyResult_WLS_CS_FI <- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_WLS_CS_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "rmsea.scaled")
  
  ## WLS_MS
  MyResult_WLS_MS_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("WLS",
                                                                                 "MS", fact,Design[RowOfDesign,3])))
  colnames(MyResult_WLS_MS_est) <- ColnamesGeneratorEst("WLS","MS",fact,Design[RowOfDesign,3])
  
  MyResult_WLS_MS_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("WLS", "MS", fact, Design[RowOfDesign,3])))
  colnames(MyResult_WLS_MS_err) <- ColnamesGeneratorSE("WLS","MS",fact,Design[RowOfDesign,3])
  
  MyResult_WLS_MS_FI <- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_WLS_MS_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "srmr")
  
  ## PML_CS
  MyResult_PML_CS_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("PML",
                                                                                 "CS", fact,Design[RowOfDesign,3])))
  colnames(MyResult_PML_CS_est) <- ColnamesGeneratorEst("PML","CS",fact,Design[RowOfDesign,3])
  
  MyResult_PML_CS_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("PML", "CS", fact, Design[RowOfDesign,3])))
  colnames(MyResult_PML_CS_err) <- ColnamesGeneratorSE("PML","CS",fact,Design[RowOfDesign,3])
  
  MyResult_PML_CS_FI <- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_PML_CS_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "srmr")
  
  ## PML_MS
  MyResult_PML_MS_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("PML",
                                                                                 "MS", fact,Design[RowOfDesign,3])))
  colnames(MyResult_PML_MS_est) <- ColnamesGeneratorEst("PML","MS",fact,Design[RowOfDesign,3])
  
  MyResult_PML_MS_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("PML", "MS", fact, Design[RowOfDesign,3])))
  colnames(MyResult_PML_MS_err) <- ColnamesGeneratorSE("PML","MS",fact,Design[RowOfDesign,3])
  
  MyResult_PML_MS_FI <- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_PML_MS_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "srmr")

  #create a loop over the replications k = 1 to K:
  tmp <- proc.time()
  for (k in 1:K){
    # Generate data
    # set a random number seed to be able to replicate the result exactly
    set.seed((k + 1000)*RowOfDesign)
    SimDat <- do.call(MyDataGeneration, Design[RowOfDesign,] )

    ## WLS correctly specified model
    fit1_W <- Method_old_CS(SimDat, fact = Design[RowOfDesign,1])

    # parameter estimates
    index <- which(fit1_W@ParTable$free != 0)
    MyAnalysisResult_WLS1est <- fit1_W@ParTable$est[index]
    MyResult_WLS_CS_est[k, ] <- MyAnalysisResult_WLS1est
    
    # standard errors 
    MyAnalysisResult_WLS1err <- fit1_W@ParTable$se[index]
    MyResult_WLS_CS_err[k, ] <- MyAnalysisResult_WLS1err
    
    ### FITINDICES
    FI_WLS_CS <- fitMeasures(fit1_W, c("chisq.scaled", "df.scaled", 
                                       "pvalue.scaled", "cfi.scaled",
                                       "srmr"))
    MyResult_WLS_CS_FI[k,] <- FI_WLS_CS
  
    ## WLS misspecified model
    fit2_W <- Method_old_MS(SimDat, fact = Design[RowOfDesign,1])
    
    # parameter estimates 
    index <- which(fit2_W@ParTable$free != 0)
    MyAnalysisResult_WLS2est <- fit2_W@ParTable$est[index]
    names(MyAnalysisResult_WLS2est) <- ColnamesGeneratorEst("WLS", 
                                                            "MS", 
                                                            fact,
                                                         Design[RowOfDesign,3])
    MyResult_WLS_MS_est[k, ] <- MyAnalysisResult_WLS2est

    # standard errors
    MyAnalysisResult_WLS2err <- fit2_W@ParTable$se[index]
    names(MyAnalysisResult_WLS2err) <- ColnamesGeneratorSE("WLS", 
                                                           "MS", 
                                                           fact, 
                                                           Design[RowOfDesign,3])
    MyResult_WLS_MS_err[k, ] <- MyAnalysisResult_WLS2err

    
    ### FITINDICES
    FI_WLS_MS <-fitMeasures(fit2_W, c("chisq.scaled","df.scaled", 
                                      "pvalue.scaled", "cfi.scaled",
                                      "srmr"))
    MyResult_WLS_MS_FI[k,] <- FI_WLS_MS
    
    ## PML correctly specified model
    fit1_P <- Method_new_CS(SimDat, fact = Design[RowOfDesign,1])
    summary(fit1_P)
    # parameter estimates
    index <- which(fit1_P@ParTable$free != 0)
    MyAnalysisResult_PML1est <- fit1_P@ParTable$est[index]
    names(MyAnalysisResult_PML1est) <- ColnamesGeneratorEst("PML", 
                                                            "CS", 
                                                            fact, 
                                                            Design[RowOfDesign,3])
    MyResult_PML_CS_est[k, ] <- MyAnalysisResult_PML1est
    
    # standard errors
    MyAnalysisResult_PML1err <- fit1_P@ParTable$se[index]
    names(MyAnalysisResult_PML1err) <- ColnamesGeneratorSE("PML", 
                                                           "CS", 
                                                           fact, 
                                                           Design[RowOfDesign,3])
    MyResult_PML_CS_err[k, ] <- MyAnalysisResult_PML1err
    
    ### FITINDICES
    FI_PML_CS <- fitMeasures(fit1_P, 
                             c("chisq.scaled","df.scaled", 
                               "pvalue.scaled", "cfi.scaled",
                               "srmr"))
    MyResult_PML_CS_FI[k,] <- FI_PML_CS
    
    #PML misspecified model
    fit2_P <- Method_new_MS(SimDat, fact = Design[RowOfDesign,1])
    summary(fit2_P)
    # parameter estimates
    index <- which(fit2_P@ParTable$free != 0)
    MyAnalysisResult_PML2est <- fit2_P@ParTable$est[index]
    names(MyAnalysisResult_PML2est) <- ColnamesGeneratorEst("PML", 
                                                            "MS", 
                                                            fact, 
                                                            Design[RowOfDesign,3])
    MyResult_PML_MS_est[k, ] <- MyAnalysisResult_PML2est
    
    # standard errors 
    MyAnalysisResult_PML2err <- fit2_P@ParTable$se[index]
    names(MyAnalysisResult_PML2err) <- ColnamesGeneratorSE("PML", 
                                                           "MS", 
                                                           fact, 
                                                           Design[RowOfDesign,3])
    MyResult_PML_MS_err[k, ] <- MyAnalysisResult_PML2err
  
    ### FITINDICES
    FI_PML_MS <- fitMeasures(fit2_P, c("chisq.scaled","df.scaled", 
                                       "pvalue.scaled", "cfi.scaled",
                                       "srmr"))
    MyResult_PML_MS_FI[k,] <- FI_PML_MS
    
  
  #save the time to run the analyses of K data sets in one cell of the design.
  time <- proc.time() - tmp
  
  }
  # save all relevant results
  return(list(MyResult_WLS_CS_est = MyResult_WLS_CS_est, 
              MyResult_WLS_CS_err = MyResult_WLS_CS_err, 
              MyResult_WLS_CS_FI = MyResult_WLS_CS_FI,
              MyResult_WLS_MS_est = MyResult_WLS_MS_est,
              MyResult_WLS_MS_err = MyResult_WLS_MS_err,
              MyResult_WLS_MS_FI = MyResult_WLS_MS_FI,
              MyResult_PML_CS_est = MyResult_PML_CS_est,
              MyResult_PML_CS_err = MyResult_PML_CS_err,
              MyResult_PML_CS_FI = MyResult_PML_CS_FI,
              MyResult_PML_MS_est = MyResult_PML_MS_est,
              MyResult_PML_MS_err = MyResult_PML_MS_err,
              MyResult_PML_MS_FI = MyResult_PML_MS_FI,
              time = time))
}

# collect data
MyResult_onecell <- MySimulationCell(Design, RowOfDesign = 1, K = 2)
MyResult_onecell

################################ Simulation all cells  ###############################
TotalCells <- nrow(Design)
for (i in 1:TotalCells){
  Row <- i
  MyResult <- MySimulationCell(Design = Design, RowOfDesign = Row, K = 2) #10!

  # Write output of one cell of the design
  # Save WLS results
  MyResult_WLS_CS_est <- MyResult[1]
  #save(MyResult1, 
       #file = paste("WLS_CS_est", "Row", Row,".Rdata" , sep = ""))
  write.csv(MyResult_WLS_CS_est, 
            file = paste("WLS_CS_est", "Row", Row,".csv" , sep = ""))

  MyResult_WLS_CS_err <- MyResult[2]
  # save(MyResult2,
  #     file = paste("WLS_CS_err", "Row", Row,".Rdata" , sep =""))
  write.csv(MyResult_WLS_CS_err,
       file = paste("WLS_CS_err", "Row", Row,".csv" , sep =""))

  MyResult_WLS_CS_FI <- MyResult[3]
  # save(MyResult3,
  #      file =paste("FI_WLS_CS", "Row", Row, ".Rdata" , sep =""))
  write.csv(MyResult_WLS_CS_FI,
        file =paste("FI_WLS_CS", "Row", Row, ".csv" , sep =""))
  
  MyResult_WLS_MS_est <- MyResult[4]
  # save(MyResult4,
  #      file = paste("WLS_MS_est", "Row", Row,".Rdata" , sep =""))
  write.csv(MyResult_WLS_MS_est,
            file = paste("WLS_MS_est", "Row", Row,".csv" , sep =""))

  MyResult_WLS_MS_err <- MyResult[5]
  # save(MyResult5,
  #      file = paste("WLS_MS_err", "Row", Row,".Rdata" , sep =""))
  write.csv(MyResult_WLS_MS_err,
       file = paste("WLS_MS_err", "Row", Row,".csv" , sep =""))

  MyResult_WLS_MS_FI <- MyResult[6]
  # save(MyResults6,
  #      file =paste("FI_WLS_MS", "Row", Row, ".Rdata" , sep =""))
  write.csv(MyResult_WLS_MS_FI,
      file =paste("FI_WLS_MS", "Row", Row, ".csv" , sep =""))

  # Save PML results
  MyResult_PML_CS_est <- MyResult[7]
  # save(MyResult7,
  #      file = paste("PML_CS_est", "Row", Row,".Rdata" , sep =""))
  write.csv(MyResult_PML_CS_est,
      file = paste("PML_CS_est", "Row", Row,".csv" , sep =""))

  MyResult_PML_CS_err <- MyResult[8]
  # save(MyResult8,
  #      file = paste("PML_CS_err", "Row", Row,".Rdata" , sep =""))
  write.csv(MyResult_PML_CS_err,
        file = paste("PML_CS_err", "Row", Row,".csv" , sep =""))
       
  MyResult_PML_CS_FI <- MyResult[9]
  # save(MyResults9,
  #      file =paste("FI_PML_CS", "Row", Row, ".Rdata" , sep =""))
  write.csv(MyResult_PML_CS_FI,
        file =paste("FI_PML_CS", "Row", Row, ".csv" , sep =""))

  MyResult_PML_MS_est <- MyResult[10]
  # save(MyResult10,
  #      file = paste("PML_MS_est", "Row", Row,".Rdata" , sep =""))
  write.csv(MyResult_PML_MS_est,
       file = paste("PML_MS_est", "Row", Row,".csv" , sep =""))
  
  MyResult_PML_MS_err <- MyResult[11]
  # save(MyResult11,
  #      file = paste("PML_MS_err", "Row", Row,".Rdata" , sep =""))
  write.csv(MyResult_PML_MS_err,
         file = paste("PML_MS_err", "Row", Row,".csv" , sep =""))

  MyResult_PML_MS_FI <- MyResult[12]
  # save(MyResults12,
  #    file =paste("FI_PML_MS", "Row", Row, ".Rdata" , sep =""))
  write.csv(MyResult_PML_MS_FI,
        file =paste("FI_PML_MS", "Row", Row, ".csv" , sep =""))

# Save time
  time <- MyResult[13]
  save(time, file =paste("Time", "Row", Row, ".Rdata" , sep =""))
  #write.csv(MyResults13, 
            #file =paste("Time", "Row", Row, ".csv" , sep =""))
}

