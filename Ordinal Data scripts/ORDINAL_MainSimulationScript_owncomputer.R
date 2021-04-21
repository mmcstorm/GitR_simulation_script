# MainSimulationScript for running on my computer (following the code of the manual)

# load packages
library(lavaan)
library(usethis)

RowOfDesign <- 1 
Replication <- 1

############################# Simulation Design  #############################
factors <- c(2,4,6,8) 					            #number of latent variables
nobs <- c(200,400,800)                      #sample size
ncat <- c(2,4)                              #number of categories

##Create the simulation design matrix (full factorial)
Design <- expand.grid(factors = factors, nobs = nobs, ncat = ncat)

# load functions
source("ORDINAL_all_functions_script.R")

################################ Simulation start (1 cell) ##########################
MySimulationCell<- function(Design = Design, RowOfDesign, Z){
  
  # initialize values
  nvarp <- 6
  fact <- Design[RowOfDesign,1]
  nvar <- nvarp*fact

  # Create matrices to store the results:
  ## WLS_withoutC
  MyResult_WLS_withoutC_est <- matrix(NA, 
                                nrow = Z, 
                                ncol = length(ColnamesGeneratorEst("WLS",
                                                                   "withoutC", 
                                                                   fact,                                      Design[RowOfDesign,3])))
  colnames(MyResult_WLS_withoutC_est) <- ColnamesGeneratorEst("WLS","withoutC",fact,Design[RowOfDesign,3])

  MyResult_WLS_withoutC_err <- matrix(NA, 
                                nrow = Z, 
                                ncol = length(ColnamesGeneratorSE("WLS", 
                                                                  "withoutC", 
                                                                   fact, 
                                                                   Design[RowOfDesign,3])))
  colnames(MyResult_WLS_withoutC_err) <- ColnamesGeneratorSE("WLS","withoutC",fact,Design[RowOfDesign,3])
  
  MyResult_WLS_withoutC_FI <- matrix(NA, nrow = Z, ncol = 5)
  colnames(MyResult_WLS_withoutC_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "rmsea.scaled")
  
  ## WLS_withC
  MyResult_WLS_withC_est <- matrix(NA, 
                                nrow = Z, 
                                ncol = length(ColnamesGeneratorEst("WLS",
                                                                   "withC", 
                                                                    fact,
                                                                    Design[RowOfDesign,3])))
  colnames(MyResult_WLS_withC_est) <- ColnamesGeneratorEst("WLS","withC",fact,Design[RowOfDesign,3])
  
  MyResult_WLS_withC_err <- matrix(NA, 
                                nrow = Z, 
                                ncol = length(ColnamesGeneratorSE("WLS", 
                                                                  "withC", 
                                                                   fact, 
                                                                   Design[RowOfDesign,3])))
  colnames(MyResult_WLS_withC_err) <- ColnamesGeneratorSE("WLS","withC",fact,Design[RowOfDesign,3])
  
  MyResult_WLS_withC_FI <- matrix(NA, nrow = Z, ncol = 5)
  colnames(MyResult_WLS_withC_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "rmsea.scaled")
  
  ## PML_withoutC
  MyResult_PML_withoutC_est <- matrix(NA, 
                                nrow = Z, 
                                ncol = length(ColnamesGeneratorEst("PML",
                                                                   "withoutC", 
                                                                    fact,
                                                                    Design[RowOfDesign,3])))
  
  colnames(MyResult_PML_withoutC_est) <- ColnamesGeneratorEst("PML","withoutC",fact,Design[RowOfDesign,3])
  
  MyResult_PML_withoutC_err <- matrix(NA, 
                                nrow = Z, 
                                ncol = length(ColnamesGeneratorSE("PML", 
                                                                  "withoutC", 
                                                                   fact, 
                                                                   Design[RowOfDesign,3])))
  colnames(MyResult_PML_withoutC_err) <- ColnamesGeneratorSE("PML","withoutC",fact,Design[RowOfDesign,3])
  
  MyResult_PML_withoutC_FI <- matrix(NA, nrow = Z, ncol = 5)
  colnames(MyResult_PML_withoutC_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "rmsea.scaled")
  
  ## PML_withC
  MyResult_PML_withC_est <- matrix(NA, 
                                nrow = Z, 
                                ncol = length(ColnamesGeneratorEst("PML",
                                                                   "withC", 
                                                                    fact,
                                                                    Design[RowOfDesign,3])))
  colnames(MyResult_PML_withC_est) <- ColnamesGeneratorEst("PML","withC",fact,Design[RowOfDesign,3])
  
  MyResult_PML_withC_err <- matrix(NA, 
                                nrow = Z, 
                                ncol = length(ColnamesGeneratorSE("PML",
                                                                  "withC",
                                                                   fact, 
                                                                   Design[RowOfDesign,3])))
  colnames(MyResult_PML_withC_err) <- ColnamesGeneratorSE("PML","withC",fact,Design[RowOfDesign,3])
  
  MyResult_PML_withC_FI <- matrix(NA, nrow = Z, ncol = 5)
  colnames(MyResult_PML_withC_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "srmr")
  
  #create a loop over the replications z = 1 to Z:
  tmp <- proc.time()

  input<- cbind(seq(1,550), rep(0,550))
  mis500<- as.matrix(input)
  
  for (z in 1:Z){
    # Generate data
    # set a random number seed to be able to replicate the result exactly
    set.seed((z + 1000)*RowOfDesign)
    SimDat <- do.call(MyDataGeneration, Design[RowOfDesign,] )
    
    ## WLS model without specified cross-loadings
    fit1_W <- try(cfa( model <- model_withoutC(fact), 
                       data=SimDat, std.lv=TRUE, 
                       ordered=c(colnames(SimDat)), #Simdat?
                       estimator="WLSMV"),silent=TRUE)
    if(inherits(fit1_W, "try-error")) {
      mis500[z,1] <- 1        
    } else { 
      # parameter estimates
      index <- which(fit1_W@ParTable$free != 0)
      MyAnalysisResult_WLS1est <- fit1_W@ParTable$est[index]
      MyResult_WLS_withoutC_est[z, ] <- MyAnalysisResult_WLS1est
    
      # standard errors 
      MyAnalysisResult_WLS1err <- fit1_W@ParTable$se[index]
      MyResult_WLS_withoutC_err[z, ] <- MyAnalysisResult_WLS1err
      
      ### FITINDICES
      FI_WLS_withoutC <- fitMeasures(fit1_W, c("chisq.scaled","df.scaled", 
                                         "pvalue.scaled", "cfi.scaled", "srmr"))
      MyResult_WLS_withoutC_FI[z,] <- FI_WLS_withoutC
    }
    
    ## WLS model with specified cross-loadings
    fit2_W <- try(cfa( model <- model_withC(fact), 
                       data=SimDat, std.lv=TRUE, 
                       ordered=c(colnames(SimDat)),
                       estimator="WLSMV"),silent=TRUE)
    if(inherits(fit2_W, "try-error")) {
      mis500[z,2] <- 1        
    } else {
      # parameter estimates 
      index <- which(fit2_W@ParTable$free != 0)
      MyAnalysisResult_WLS2est <- fit2_W@ParTable$est[index]
      names(MyAnalysisResult_WLS2est) <- ColnamesGeneratorEst("WLS",
                                                              "withC",
                                                              fact,
                                                              Design[RowOfDesign,3])
      MyResult_WLS_withC_est[z, ] <- MyAnalysisResult_WLS2est
      
      # standard errors
      MyAnalysisResult_WLS2err <- fit2_W@ParTable$se[index]
      names(MyAnalysisResult_WLS2err) <- ColnamesGeneratorSE("WLS", 
                                                             "withC",
                                                              fact,
                                                              Design[RowOfDesign,3])
      
      MyResult_WLS_withC_err[z, ] <- MyAnalysisResult_WLS2err
      
      ### FITINDICES
      FI_WLS_withC <-fitMeasures(fit2_W, c("chisq.scaled","df.scaled", 
                                        "pvalue.scaled", "cfi.scaled", "srmr"))
      MyResult_WLS_withC_FI[z,] <- FI_WLS_withC
    }
    
    ## PML model without specified cross-loadings
    fit1_P <- try(cfa( model <- model_withoutC(fact), 
                       data=SimDat, std.lv=TRUE, 
                       ordered=c(colnames(SimDat)), #Simdat?
                       estimator="PML"),silent=TRUE)
    if(inherits(fit1_P, "try-error")) {
      mis500[z,2] <- 1        
    } else {
      # parameter estimates
      index <- which(fit1_P@ParTable$free != 0)
      MyAnalysisResult_PML1est <- fit1_P@ParTable$est[index]
      names(MyAnalysisResult_PML1est) <- ColnamesGeneratorEst("PML", 
                                                            "withoutC", 
                                                            fact, 
                                                            Design[RowOfDesign,3])
      MyResult_PML_withoutC_est[z, ] <- MyAnalysisResult_PML1est
    
      # standard errors
      MyAnalysisResult_PML1err <- fit1_P@ParTable$se[index]
      names(MyAnalysisResult_PML1err) <- ColnamesGeneratorSE("PML", 
                                                             "withoutC",
                                                              fact, 
                                                              Design[RowOfDesign,3])
      MyResult_PML_withoutC_err[z, ] <- MyAnalysisResult_PML1err
      
      ### FITINDICES
      FI_PML_withoutC <- fitMeasures(fit1_P, 
                               c("chisq.scaled","df.scaled", 
                                 "pvalue.scaled", "cfi.scaled", "srmr"))
      MyResult_PML_withoutC_FI[z,] <- FI_PML_withoutC
    }
    
    #PML model with specified cross-loadings
    fit2_P <- try(cfa( model <- model_withC(fact), 
                       data=SimDat, std.lv=TRUE, 
                       ordered=c(colnames(SimDat)), 
                       estimator="PML"),silent=TRUE)
    if(inherits(fit2_P, "try-error")) {
      mis500[z,2] <- 1        
      z <- z+1
    } else {
      # parameter estimates
      index <- which(fit2_P@ParTable$free != 0)
      MyAnalysisResult_PML2est <- fit2_P@ParTable$est[index]
      names(MyAnalysisResult_PML2est) <- ColnamesGeneratorEst("PML", 
                                                              "withC", 
                                                              fact, 
                                                              Design[RowOfDesign,3])
      MyResult_PML_withC_est[z, ] <- MyAnalysisResult_PML2est
      
      # standard errors 
      MyAnalysisResult_PML2err <- fit2_P@ParTable$se[index]
      names(MyAnalysisResult_PML2err) <- ColnamesGeneratorSE("PML", 
                                                             "withC",
                                                             fact, 
                                                             Design[RowOfDesign,3])
      MyResult_PML_withC_err[z, ] <- MyAnalysisResult_PML2err
      
      ### FITINDICES
      FI_PML_withC <- fitMeasures(fit2_P, c("chisq.scaled","df.scaled", 
                                         "pvalue.scaled", "cfi.scaled", "srmr"))
      MyResult_PML_withC_FI[z, ] <- FI_PML_withC
    }
    
    #save the time to run the analyses of K data sets in one cell of the design.
    time <- proc.time() - tmp
    
  }
  # save all relevant results
  return(list(MyResult_WLS_withoutC_est = MyResult_WLS_withoutC_est, 
              MyResult_WLS_withoutC_err = MyResult_WLS_withoutC_err, 
              MyResult_WLS_withoutC_FI = MyResult_WLS_withoutC_FI,
              MyResult_WLS_withC_est = MyResult_WLS_withC_est,
              MyResult_WLS_withC_err = MyResult_WLS_withC_err,
              MyResult_WLS_withC_FI = MyResult_WLS_withC_FI,
              MyResult_PML_withoutC_est = MyResult_PML_withoutC_est,
              MyResult_PML_withoutC_err = MyResult_PML_withoutC_err,
              MyResult_PML_withoutC_FI = MyResult_PML_withoutC_FI,
              MyResult_PML_withC_est = MyResult_PML_withC_est,
              MyResult_PML_withC_err = MyResult_PML_withC_err,
              MyResult_PML_withC_FI = MyResult_PML_withC_FI,
              time = time,
              mis500 = mis500))
}

# collect data
MyResult_onecell <- MySimulationCell(Design, RowOfDesign = 2, Z = 2)
MyResult_onecell
str(MyResult_onecell)

################################ Simulation all cells  ###############################
TotalCells <- nrow(Design)
for (i in 1:TotalCells){
  Row <- i
  MyResult <- MySimulationCell(Design = Design, RowOfDesign = Row, Z = 2) #10!
  
  # Write output of one cell of the design
  # Save WLS results
  MyResult1 <- MyResult[1]
  write.csv(MyResult1, 
            file = paste("WLS_withoutC_est", "Row", Row,".csv" , sep = ""))
  
  MyResult2 <- MyResult[2]
  write.csv(MyResult2,
            file = paste("WLS_withoutC_err", "Row", Row,".csv" , sep =""))
  
  MyResult3 <- MyResult[3]
  write.csv(MyResult3,
            file =paste("FI_WLS_withoutC", "Row", Row, ".csv" , sep =""))
  
  MyResult4 <- MyResult[4]
  write.csv(MyResult4,
            file = paste("WLS_withC_est", "Row", Row,".csv" , sep =""))
  
  MyResult5 <- MyResult[5]
  write.csv(MyResult5,
            file = paste("WLS_withC_err", "Row", Row,".csv" , sep =""))
  
  MyResults6 <- MyResult[6]
  write.csv(MyResults6,
            file =paste("FI_WLS_withC", "Row", Row, ".csv" , sep =""))
  
  # Save PML results
  MyResult7 <- MyResult[7]
  write.csv(MyResult7,
            file = paste("PML_withoutC_est", "Row", Row,".csv" , sep =""))
  
  MyResult8 <- MyResult[8]
  write.csv(MyResult8,
            file = paste("PML_withoutC_err", "Row", Row,".csv" , sep =""))
  
  MyResults9 <- MyResult[9]
  write.csv(MyResults9,
            file =paste("FI_PML_withoutC", "Row", Row, ".csv" , sep =""))
  
  MyResult10 <- MyResult[10]
  write.csv(MyResult10,
            file = paste("PML_withC_est", "Row", Row,".csv" , sep =""))
  
  MyResult11 <- MyResult[11]
  write.csv(MyResult11,
            file = paste("PML_withC_err", "Row", Row,".csv" , sep =""))
  
  MyResults12 <- MyResult[12]
  write.csv(MyResults12,
            file =paste("FI_PML_withC", "Row", Row, ".csv" , sep =""))
  
  # Save time
  MyResults13 <- MyResult[13]
  save(MyResults13, file =paste("Time", "Row", Row, ".Rdata" , sep =""))
  
}
