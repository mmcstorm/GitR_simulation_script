# MainSimulationScript for running on my computer (following the code of the manual)

# load packages
library(lavaan)
library(usethis)

RowOfDesign <- 1 #24
Replication <- 1

############################# Simulation Design  #############################
factors <- c(2,4,6,8) 					            #number of latent variables
nobs <- c(200,400,800)                      #sample size
##Create the simulation design matrix (full factorial)
Design_mixed <- expand.grid(factors = factors,nobs = nobs)

# load functions
source("MIXED_all_functions_script.R")

################################ Simulation start (1 cell) ##########################

MySimulationCell<- function(Design_mixed = Design_mixed, RowOfDesign, K){
  
  # initialize values
  nvarp <- 6
  fact <- Design_mixed[RowOfDesign,1]
  nvar <- nvarp*fact
  
  ##### WLS ####
  # model without specified cross-loadings
  MyResult_WLS_withoutC_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("WLS", "withoutC", fact)))

  colnames(MyResult_WLS_withoutC_est) <- ColnamesGeneratorEst("WLS", "withoutC", fact)
  
  MyResult_WLS_withoutC_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("WLS", "withoutC", fact)))
  colnames(MyResult_WLS_withoutC_err) <- ColnamesGeneratorSE("WLS", "withoutC", fact)
  
  MyResult_WLS_withoutC_FI<- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_WLS_withoutC_FI) <- c("chisq","df", 
                                          "pvalue", "cfi",
                                          "srmr")
  
  # model with specified cross-loadings
  MyResult_WLS_withC_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("WLS", "withC", fact)))

  colnames(MyResult_WLS_withC_est) <- ColnamesGeneratorEst("WLS", "withC", fact)
  
  MyResult_WLS_withC_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("WLS", "withC", fact)))
  colnames(MyResult_WLS_withC_err) <- ColnamesGeneratorSE("WLS", "withC", fact)
  
  MyResult_WLS_withC_FI <- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_WLS_withC_FI) <- c("chisq","df", 
                                       "pvalue", "cfi",
                                       "srmr")
  
  ##### PML ####
  MyResult_PML_withoutC_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("PML", "withoutC", fact)))
  
  colnames(MyResult_PML_withoutC_est) <- ColnamesGeneratorEst("PML", "withoutC", fact)
  
  MyResult_PML_withoutC_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("PML", "withoutC", fact)))
  colnames(MyResult_PML_withoutC_err) <- ColnamesGeneratorSE("PML", "withoutC", fact)
  
  MyResult_PML_withoutC_FI<- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_PML_withoutC_FI) <- c("chisq.scaled", "df.scaled", 
                                          "pvalue.scaled", "cfi.scaled",
                                          "srmr")
  
  # model with specified cross-loadings
  MyResult_PML_withC_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("PML", "withC", fact)))
  
  colnames(MyResult_PML_withC_est) <- ColnamesGeneratorEst("PML", "withC", fact)
  
  MyResult_PML_withC_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("PML", "withC", fact)))
  colnames(MyResult_PML_withC_err) <- ColnamesGeneratorSE("PML", "withC", fact)
  
  MyResult_PML_withC_FI <- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_PML_withC_FI) <- c("chisq.scaled", "df.scaled", 
                                       "pvalue.scaled", "cfi.scaled",
                                       "srmr")
  
  input<- cbind(seq(1,550), rep(0,550))
  mis500<- as.matrix(input)
  nvarp <- 6
  nvar <- Design_mixed[RowOfDesign,1]*nvarp
  
  #create a loop over the replications k = 1 to K:
  tmp <- proc.time()
  for (k in 1:K){
    # Generate data
    # set a random number seed to be able to replicate the result exactly
    set.seed((k + 1000)*RowOfDesign)

    SimDat <- do.call(MyDataGeneration, Design_mixed[RowOfDesign,] )
    
    # fit model without specified cross-loadings (with silent check)

    
    #### WLS ####
    fit_WLS_withoutC <- try(cfa( model <- model_withoutC(fact), 
                                 data=SimDat, std.lv=TRUE, 
                                 ordered = colnames(SimDat[indexes_ord(fact)]),
                                 estimator="WLSMV"),silent=TRUE)
    if(inherits(fit_WLS_withoutC, "try-error")) {
      mis500[k,1] <- 1        
    } else { 
      
      # parameter estimates
      index <- which(fit_WLS_withoutC@ParTable$free != 0)
      MyAnalysisResult_WLS_withoutC_est <- fit_WLS_withoutC@ParTable$est[index]
      names(MyAnalysisResult_WLS_withoutC_est) <- ColnamesGeneratorEst("WLS", "withoutC", fact)
      MyResult_WLS_withoutC_est[k, ] <- MyAnalysisResult_WLS_withoutC_est
      
      # standard errors
      MyAnalysisResult_WLS_withoutC_err <- fit_WLS_withoutC@ParTable$se[index]
      names(MyAnalysisResult_WLS_withoutC_err) <- ColnamesGeneratorSE("WLS", "withoutC", fact)
      MyResult_WLS_withoutC_err[k, ] <- MyAnalysisResult_WLS_withoutC_err
      
      ### FITINDICES
      FI_WLS_withoutC <- fitMeasures(fit_WLS_withoutC, 
                                     c("chisq","df", 
                                       "pvalue", "cfi",
                                       "srmr"))
      MyResult_WLS_withoutC_FI[k,] <- FI_WLS_withoutC
    }

    #### Model with specified cross-loadings (with silent check)
    fit_WLS_withC <- try(cfa( model <- model_withC(fact), 
                              data=SimDat, std.lv=TRUE, 
                              ordered = colnames(SimDat[indexes_ord(fact)]),
                              estimator="WLSMV"),silent=TRUE)
    if(inherits(fit_WLS_withC, "try-error")) {
      mis500[k,2] <- 1        
    } else {
      # parameter estimates
      index <- which(fit_WLS_withC@ParTable$free != 0)
      MyAnalysisResult_WLS_withC_est <- fit_WLS_withC@ParTable$est[index]
      names(MyAnalysisResult_WLS_withC_est) <- ColnamesGeneratorEst("WLS", "withC", fact)
      MyResult_WLS_withC_est[k, ] <- MyAnalysisResult_WLS_withC_est
      
      # standard errors 
      MyAnalysisResult_WLS_withC_err <- fit_WLS_withC@ParTable$se[index]
      names(MyAnalysisResult_WLS_withC_err) <- ColnamesGeneratorSE("WLS", "withC", fact)
      MyResult_WLS_withC_err[k, ] <- MyAnalysisResult_WLS_withC_err
      
      ### FITINDICES
      FI_WLS_withC <- fitMeasures(fit_WLS_withC, c("chisq","df", 
                                                    "pvalue", "cfi",
                                                    "srmr"))
      MyResult_WLS_withC_FI[k,] <- FI_WLS_withC
    }
    
    
    
    #### PML ####
    fit_PML_withoutC <- try(cfa( model <- model_withoutC(fact), 
             data=SimDat, std.lv=TRUE, 
             ordered = colnames(SimDat[indexes_ord(fact)]),
             estimator="PML"),silent=TRUE)
    if(inherits(fit_PML_withoutC, "try-error")) {
      mis500[k,1] <- 1        
    } else { 

    # parameter estimates
    index <- which(fit_PML_withoutC@ParTable$free != 0)
    MyAnalysisResult_PML_withoutC_est <- fit_PML_withoutC@ParTable$est[index]
    names(MyAnalysisResult_PML_withoutC_est) <- ColnamesGeneratorEst("PML", "withoutC", fact)
    MyResult_PML_withoutC_est[k, ] <- MyAnalysisResult_PML_withoutC_est
    
    # standard errors
    MyAnalysisResult_PML_withoutC_err <- fit_PML_withoutC@ParTable$se[index]
    names(MyAnalysisResult_PML_withoutC_err) <- ColnamesGeneratorSE("PML", "withoutC", fact)
    MyResult_PML_withoutC_err[k, ] <- MyAnalysisResult_PML_withoutC_err
    
    ### FITINDICES
    FI_PML_withoutC <- fitMeasures(fit_PML_withoutC, 
                             c("chisq.scaled","df.scaled", 
                               "pvalue.scaled", "cfi.scaled",
                               "srmr"))
    MyResult_PML_withoutC_FI[k,] <- FI_PML_withoutC
    }
    
    #### Model with specified cross-loadings (with silent check)
    fit_PML_withC <- try(cfa( model <- model_withC(fact), 
                       data=SimDat, std.lv=TRUE, 
                       ordered = colnames(SimDat[indexes_ord(fact)]),
                       estimator="PML"),silent=TRUE)
    if(inherits(fit_PML_withC, "try-error")) {
      mis500[k,2] <- 1        
    } else {
    # parameter estimates
    index <- which(fit_PML_withC@ParTable$free != 0)
    MyAnalysisResult_PML_withC_est <- fit_PML_withC@ParTable$est[index]
    names(MyAnalysisResult_PML_withC_est) <- ColnamesGeneratorEst("PML", "withC", fact)
    MyResult_PML_withC_est[k, ] <- MyAnalysisResult_PML_withC_est
    
    # standard errors 
    MyAnalysisResult_PML_withC_err <- fit_PML_withC@ParTable$se[index]
    names(MyAnalysisResult_PML_withC_err) <- ColnamesGeneratorSE("PML", "withC", fact)
    MyResult_PML_withC_err[k, ] <- MyAnalysisResult_PML_withC_err
    
    ### FITINDICES
    FI_PML_withC <- fitMeasures(fit_PML_withC, c("chisq.scaled","df.scaled", 
                                       "pvalue.scaled", "cfi.scaled",
                                       "srmr"))
    MyResult_PML_withC_FI[k,] <- FI_PML_withC
    }
  }
    #save the time to run the analyses of K data sets in one cell of the Design_mixed.
    time <- proc.time() - tmp
  
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
              Simulated_Data = SimDat,
              time = time, mis500 = mis500))
}

# collect data
MyResult_onecell <- MySimulationCell(Design_mixed, RowOfDesign = 2, K = 2)
MyResult_onecell

################################ Simulation all cells  ###############################
TotalCells <- nrow(Design_mixed)
for (i in 1:TotalCells){
  Row <- i
  MyResult <- MySimulationCell(Design_mixed = Design_mixed, RowOfDesign = Row, K = 2) #10!

  # Write output of one cell of the design
  # Save results
  
  #### WLS ####
  MyResult_WLS_withoutC_est <- MyResult[1]
  write.csv(MyResult_WLS_withoutC_est, 
            file = paste("WLS_withoutC_est", "Row", Row,".csv" , sep = ""))
  
  MyResult_WLS_withoutC_err <- MyResult[2]
  write.csv(MyResult_WLS_withoutC_err,
            file = paste("WLS_withoutC_err", "Row", Row,".csv" , sep =""))
  
  MyResult_WLS_withoutC_FI <- MyResult[3]
  write.csv(MyResult_WLS_withoutC_FI,
            file =paste("WLS_FI_withoutC", "Row", Row, ".csv" , sep =""))
  
  # Model with specified cross-loadings
  MyResult_WLS_withC_est <- MyResult[4]
  write.csv(MyResult_WLS_withC_est,
            file = paste("WLS_withC_est", "Row", Row,".csv" , sep =""))
  
  MyResult_WLS_withC_err <- MyResult[5]
  write.csv(MyResult_WLS_withC_err,
            file = paste("WLS_withC_err", "Row", Row,".csv" , sep =""))
  
  MyResult_WLS_withC_FI <- MyResult[6]
  write.csv(MyResult_WLS_withC_FI,
            file =paste("WLS_FI_withC", "Row", Row, ".csv" , sep =""))
  
  #### PML ####
  MyResult_PML_withoutC_est <- MyResult[7]
  write.csv(MyResult_PML_withoutC_est, 
            file = paste("PML_withoutC_est", "Row", Row,".csv" , sep = ""))
  
  MyResult_PML_withoutC_err <- MyResult[8]
  write.csv(MyResult_PML_withoutC_err,
            file = paste("PML_withoutC_err", "Row", Row,".csv" , sep =""))
  
  MyResult_PML_withoutC_FI <- MyResult[9]
  write.csv(MyResult_PML_withoutC_FI,
            file =paste("PML_FI_withoutC", "Row", Row, ".csv" , sep =""))
  
  # Model with specified cross-loadings
  MyResult_PML_withC_est <- MyResult[10]
  write.csv(MyResult_PML_withC_est,
            file = paste("PML_withC_est", "Row", Row,".csv" , sep =""))
  
  MyResult_PML_withC_err <- MyResult[11]
  write.csv(MyResult_PML_withC_err,
            file = paste("PML_withC_err", "Row", Row,".csv" , sep =""))
  
  MyResult_PML_withC_FI <- MyResult[12]
  write.csv(MyResult_PML_withC_FI,
            file =paste("PML_FI_withC", "Row", Row, ".csv" , sep =""))
  
  # Save dataset
  Simulated_data <- MyResult[13]
  write.csv(Simulated_data, 
            file = paste("Simulated_Data", 'Row', Row, ".csv", sep = ""))
  
  # Save time
  time <- MyResult[14]
  save(time, file =paste("Time", "Row", Row, ".Rdata" , sep =""))
  
  # Silent Check
  mis500 <- MyResult[15]
  write.csv(mis500, 
            file = paste("mis500", 'Row', Row, ".csv", sep = ""))
}
