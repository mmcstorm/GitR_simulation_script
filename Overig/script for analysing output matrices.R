### Script for analysing output matrices

# MainSimulationScript for running on my computer (following the code of the manual)

# load packages
library(lavaan)
library(usethis)

args <- commandArgs(TRUE) #?
args <- as.numeric(args)

RowOfDesign <- args[1]
Replication <- args[2]

RowOfDesign <- 1 #24
Replication <- 1

############################# Simulation Design  #############################
factors <- c(2,4,6,8) 					           #number of latent variables
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
    
    # retrieve estimated parameters
    
    ## WLS correctly specified model
    fit1_W <- Method_old_CS(SimDat, fact = Design[RowOfDesign,1])
    
    # parameter estimates
    lamb <- data.frame(cbind(fit1_W@Model@GLIST$lambda))
    se_lambdas <- fit1_W@ParTable$se[1:nvar]
    L_SE <- cbind(lamb, se_lambdas)
    thresholds <- as.vector(fit1_W@Model@GLIST$tau)
    variances <- unlist(sapply(fit1_W@Model@GLIST$theta, 
                               function(x) x[x != 0]))
    correlations <- unique(unlist(sapply(fit1_W@Model@GLIST$psi, 
                                         function(x) x[x != 1])))
    
    errors <- NULL
    #standard errors 
    errors <- fit1_W@ParTable$se
    
    
    se_thresholds <- fit1_W@ParTable$se[(nvar+1):(nvar*2)]
    se_variances <- rep(0, nvar)
    se_correlations <- fit1_W@ParTable$se[nvar*3+fact+1]
    est_parameter <- cbind(L_SE, thresholds, se_thresholds, variances)
    correlations <- cbind(correlations, se_correlations)
    MyAnalysisResult_WLS1 <- list(est_parameter, correlations)
    
    ## WLS misspecified model
    fit2_W <- Method_old_MS(SimDat, fact = Design[RowOfDesign,1])
    
    #parameter estimates
    lamb <- data.frame(cbind(fit2_W@Model@GLIST$lambda))
    se_lambdas <- fit2_W@ParTable$se[1:nvar]
    L_SE <- cbind(lamb, se_lambdas)
    
    thresholds <- as.vector(fit2_W@Model@GLIST$tau)
    variances <- unlist(sapply(fit2_W@Model@GLIST$theta, 
                               function(x) x[x != 0]))
    correlations <- unique(unlist(sapply(fit2_W@Model@GLIST$psi, 
                                         function(x) x[x != 1])))
    
    #standard errors
    se_thresholds <- fit2_W@ParTable$se[(nvar+1):(nvar*2)]
    se_variances <- rep(0, nvar)
    se_correlations <- fit2_W@ParTable$se[nvar*3+fact+2]
    est_parameter <- cbind(L_SE, thresholds, se_thresholds, variances)
    correlations <- cbind(correlations, se_correlations)
    MyAnalysisResult_WLS2 <- list(est_parameter, correlations)
    
    ## PML correctly specified model
    fit1_P <- Method_new_CS(SimDat, fact = Design[RowOfDesign,1])
    lamb <- data.frame(cbind(fit1_P@Model@GLIST$lambda))
    se_lambdas <- fit1_P@ParTable$se[1:nvar]
    L_SE <- cbind(lamb, se_lambdas)
    
    thresholds <- as.vector(fit1_P@Model@GLIST$tau)
    variances <- unlist(sapply(fit1_P@Model@GLIST$theta, 
                               function(x) x[x != 0]))
    correlations <- unique(unlist(sapply(fit1_P@Model@GLIST$psi, 
                                         function(x) x[x != 1])))
    
    se_thresholds <- fit1_P@ParTable$se[(nvar+1):(nvar*2)]
    se_variances <- rep(0, nvar)
    se_correlations <- fit1_P@ParTable$se[nvar*3+fact+1]
    est_parameter <- cbind(L_SE, thresholds, se_thresholds, variances)
    correlations <- cbind(correlations, se_correlations)
    MyAnalysisResult_PML1 <- list(est_parameter, correlations)
    
    #PML misspecified model
    fit2_P <- Method_new_MS(SimDat, fact = Design[RowOfDesign,1])
    
    #parameter estimates
    lamb <- data.frame(cbind(fit2_P@Model@GLIST$lambda))
    se_lambdas <- fit2_P@ParTable$se[1:nvar]
    L_SE <- cbind(lamb, se_lambdas)
    
    thresholds <- as.vector(fit2_P@Model@GLIST$tau)
    variances <- unlist(sapply(fit2_P@Model@GLIST$theta, 
                               function(x) x[x != 0]))
    correlations <- unique(unlist(sapply(fit2_P@Model@GLIST$psi, 
                                         function(x) x[x != 1])))
    
    #standard errors
    se_lambdas <- fit2_P@ParTable$se[1:nvar]
    se_thresholds <- fit2_P@ParTable$se[(nvar+1):(nvar*2)]
    se_variances <- rep(0, nvar)
    se_correlations <- fit2_P@ParTable$se[nvar*3+fact+2] #adjust this
    est_parameter <- cbind(L_SE, thresholds, se_thresholds, variances)
    correlations <- cbind(correlations, se_correlations)
    MyAnalysisResult_PML2 <- list(est_parameter, correlations)
    
    # concatenate all results
    MyAnalysisResult <- list(WLS_CS = MyAnalysisResult_WLS1, 
                             WLS_MS = MyAnalysisResult_WLS2, 
                             PML_CS = MyAnalysisResult_PML1, 
                             PML_MS = MyAnalysisResult_PML2)
    #Evaluate the analysis results of Method_new (Result1) and Mehtod_old (Result2)
    #MyResult1 <- MyEvaluationPC(MyAnalysisResult1)
    #MyResult2 <- MyEvaluationPC(MyAnalysisResult2)
    #store the results in the right row k of your result matrix:
    #We only store the second result which is the evaluation criterion
    #MyResult[k, ] <- MyAnalysisResult
    #colnames(MyResult) <- colnames(MyAnalysisResult)
    #rownames(MyResult) <- rownames(MyAnalysisResult)
  }
  #save the time to run the analyses of K data sets in one cell of the design.
  time <- proc.time() - tmp
  return(MyAnalysisResult)
}

