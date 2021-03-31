# MainSimulationScript for running on my computer (following the code of the manual)

# load packages
library(lavaan)
library(usethis)

args <- commandArgs(TRUE) # SLURM language
args <- as.numeric(args)

RowOfDesign <- args[1]
Replication <- args[2]

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
source("Method_new.R")
source("Method_old.R")
source("ColnamesGenerators.R")

################################ Simulation start (1 cell) ##########################
MySimulationCell<- function(Design = Design, RowOfDesign = 1, K = 2){
  # Input arguments:
  # Design = designmatrix
  # RowOfDesign: number that refers to the row of the design matrix = one cell
  # K: Total number of replications = number of data sets generated in one cell
  # Create matrix or dataframe to store the results:
  MyResult <- matrix(NA, nrow = K, ncol = 4)
  
  # initialize values
  nvarp <- 6
  fact <- Design[RowOfDesign,1]
  nvar <- nvarp*fact
  
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
    names(MyAnalysisResult_WLS1est) <- ColnamesGeneratorEst("WLS", "CS", fact)
    
    # standard errors 
    MyAnalysisResult_WLS1err <- fit1_W@ParTable$se[index]
    names(MyAnalysisResult_WLS1err) <- ColnamesGeneratorSE("WLS", "CS", fact)
    
    ### FITINDICES
    FI_WLS_CS <- fitMeasures(fit1_W, c("chisq.scaled","df.scaled", 
                                       "pvalue.scaled", "cfi.scaled",
                                       "rmsea.scaled","srmr.scaled", 
                                       "nfi.scaled"))
  
    ## WLS misspecified model
    fit2_W <- Method_old_MS(SimDat, fact = Design[RowOfDesign,1])
    
    # parameter estimates 
    index <- which(fit2_W@ParTable$free != 0)
    MyAnalysisResult_WLS2est <- fit2_W@ParTable$est[index]
    names(MyAnalysisResult_WLS2est) <- ColnamesGeneratorEst("WLS", "MS", fact)
    
    # standard errors
    MyAnalysisResult_WLS2err <- fit2_W@ParTable$se[index]
    names(MyAnalysisResult_WLS2err) <- ColnamesGeneratorSE("WLS", "MS", fact)
    
    ### FITINDICES
    FI_WLS_MS <-fitMeasures(fit2_W, c("chisq.scaled","df.scaled", 
                                      "pvalue.scaled", "cfi.scaled",
                                      "rmsea.scaled","srmr.scaled", 
                                      "nfi.scaled"))
    
    ## PML correctly specified model
    fit1_P <- Method_new_CS(SimDat, fact = Design[RowOfDesign,1])
    
    # parameter estimates
    index <- which(fit1_P@ParTable$free != 0)
    MyAnalysisResult_PML1est <- fit1_P@ParTable$est[index]
    names(MyAnalysisResult_PML1est) <- ColnamesGeneratorEst("PML", "CS", fact)
    
    # standard errors
    MyAnalysisResult_PML1err <- fit1_P@ParTable$se[index]
    names(MyAnalysisResult_PML1err) <- ColnamesGeneratorSE("PML", "CS", fact)
    
    ### FITINDICES
    FI_PML_CS <- fitMeasures(fit1_P, 
                             c("chisq.scaled","df.scaled", 
                               "pvalue.scaled", "cfi.scaled",
                               "rmsea.scaled","srmr.scaled", 
                               "nfi.scaled"))
    fitMeasures(fit1_W)
    #PML misspecified model
    fit2_P <- Method_new_MS(SimDat, fact = Design[RowOfDesign,1])
    
    # parameter estimates
    index <- which(fit2_P@ParTable$free != 0)
    MyAnalysisResult_PML2est <- fit2_P@ParTable$est[index]
    names(MyAnalysisResult_PML2est) <- ColnamesGeneratorEst("PML", "MS", fact)
    
    # standard errors 
    MyAnalysisResult_PML2err <- fit2_P@ParTable$se[index]
    names(MyAnalysisResult_PML2err) <- ColnamesGeneratorSE("PML", "MS", fact)
  
    ### FITINDICES
    FI_PML_MS <- fitMeasures(fit2_P, c("chisq.scaled","df.scaled", 
                                       "pvalue.scaled", "cfi.scaled",
                                       "rmsea.scaled","srmr.scaled", 
                                       "nfi.scaled"))
  }
  
  #save the time to run the analyses of K data sets in one cell of the design.
  time <- proc.time() - tmp
  # save all relevant results
  return(list(MyAnalysisResult_WLS1est, 
              MyAnalysisResult_WLS1err, 
              FI_WLS_CS, 
              MyAnalysisResult_WLS2est, 
              MyAnalysisResult_WLS2err, 
              FI_WLS_MS, 
              MyAnalysisResult_PML1est, 
              MyAnalysisResult_PML1err, 
              FI_PML_CS, 
              MyAnalysisResult_PML2est, 
              MyAnalysisResult_PML2err, 
              FI_PML_MS,
              time))
}

# collect data
#Row <- 1
MyResult_onecell <- MySimulationCell(Design, RowOfDesign = 24, K = 1)
MyResult_onecell

################################ Simulation all cells  ###############################
TotalCells <- nrow(Design)
for (i in 1:TotalCells){
  Row <- i
  MyResult <- MySimulationCell(Design = Design, RowOfDesign = Row, K = 1) #10!
  
  # Write output of one cell of the design
  # Save WLS results
  MyResult1 <- unlist(MyResult[1])
  matrix1 == MyResult1
  save(MyResult1, 
       file = paste("WLS_CS_est", "Row", Row,".Rdata" , sep = ""))
  
  MyResult2 <- unlist(MyResult[2])
  save(MyResult2, 
      file = paste("WLS_CS_err", "Row", Row,".Rdata" , sep =""))
  
  MyResult3 <- unlist(MyResult[3])
  save(MyResult3, 
       file = paste("WLS_MS_est", "Row", Row,".Rdata" , sep =""))
  
  MyResult4 <- unlist(MyResult[4])
  save(MyResult4, 
       file = paste("WLS_MS_err", "Row", Row,".Rdata" , sep =""))
  
  # Save PML results
  MyResult5 <- unlist(MyResult[5])
  save(MyResult5, 
       file = paste("PML_CS_est", "Row", Row,".Rdata" , sep =""))
  
  MyResult6 <- unlist(MyResult[6])
  save(MyResult6, 
       file = paste("PML_CS_err", "Row", Row,".Rdata" , sep =""))
  
  MyResult7 <- unlist(MyResult[7])
  save(MyResult7, 
       file = paste("PML_MS_est", "Row", Row,".Rdata" , sep =""))
  
  MyResult8 <- unlist(MyResult[8])
  save(MyResult8, 
       file = paste("PML_MS_err", "Row", Row,".Rdata" , sep =""))

  # Save time 
  save(time, file =paste("Time", "Row", Row, ".Rdata" , sep =""))

  # Save Fit indices
  MyResult9 <- unlist(matrixf1)
save(MyResult9, 
     file =paste("FI_WLS_CS", "Row", Row, ".Rdata" , sep =""))

  MyResults10 <- unlist(matrixf2)
  save(MyResults10, 
     file =paste("FI_WLS_MS", "Row", Row, ".Rdata" , sep =""))
  
  MyResults11 <- unlist(matrixf3)
save(MyResults11, 
     file =paste("FI_PML_CS", "Row", Row, ".Rdata" , sep =""))
  
  MyResults12 <- unlist(matrixf4)
save(MyResults12, 
     file =paste("FI_PML_MS", "Row", Row, ".Rdata" , sep =""))
}