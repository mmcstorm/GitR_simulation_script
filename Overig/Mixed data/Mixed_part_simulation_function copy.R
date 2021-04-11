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
data_type <- c("ordinal", "continuous")

Design_mixed <- expand.grid(factors = factors,nobs = nobs)
##Create the simulation design matrix (full factorial)
# Design is a data.frame with all possible combinations of the factor levels
# Each row of the design matrix represents a cell of your simulation design
# Design <- expand.grid(factors = factors, nobs = nobs, ncat = ncat)

# prepared functions
MyDataGeneration <- function(factors, nobs, nvarp = 6) {

  # model specifications
  nvar <- factors*nvarp
  BM <- matrix(c(.8,.8,.7,.7,.6,.6),nrow=1) #factor loadings
  r <- 0.3 # correlation between latent variables
  int <- c(0,0,0,0,0,0) # intercept
  s <- matrix(r, factors, factors) 
  diag(s) <- 1

  # intercepts of the items
  int2 <- rep(int,factors)

  # b is the loadingsmatrix > transposed matrix version
  b <- t(kronecker(diag(factors), BM)) # kronecker computes product of two arrays 
  
  # add crossloadings (misspecifications)
  x.mis <- sapply(1:(factors/2), function(x) return(6 + 12*(x-1)))
  y.mis <- sapply(1:(factors/2), function(x) return(2 + 2*(x-1)))
  b[x.mis, y.mis] <- .2
  
  #compute error values (theta matrix)
  ev <- diag(1 - diag(b %*% s %*% t(b)))
  
  #compute sigma (variance-covariance matrix) 
  SIGMA <- b %*% s %*% t(b) + ev  #1's on the diagonal

  #simulate data from a Multivariate Normal Distribution
  x <- data.frame(MASS:::mvrnorm(n = nobs, mu = int2, Sigma = SIGMA))
  
  item.cutpoints <- 
    matrix(c(rep(-1.20,nvar),rep(0,nvar),
             rep(1.20,nvar)), ncol=nvar, byrow=TRUE)
  
  #item cutpoints (add boundaries)
  item.cutpoints <- rbind(rep(-Inf, ncol(item.cutpoints)), 
                          item.cutpoints, rep(Inf, ncol(item.cutpoints)))

  for(i in 1:(ncol(item.cutpoints)/2)){ 
    x[,i] = cut(x[,i], br=item.cutpoints[,i], 
                labels=FALSE, include.lowest=TRUE)
  }
  
  
  #x is an nobs by nvars matrix with item scores
  return(x)
}
Design_mixed
# test function
MyDataGeneration(factors = 8, nobs = 200)

# Colnames Generator 
ColnamesGeneratorEst <- function(specif, data_type, facts){
  if(specif == "MS"){
    x.mis <- sapply(1:(facts/2), function(x) return(6 + 12*(x-1))) + seq(1, facts/2, 1)
    len <- facts*6.5
    lambs <- character(len)
    lambs[-x.mis] <- paste("L", 1:(facts*6), sep = "")
    lambs[x.mis] <- paste("L", letters[1:(facts/2)], sep = "")
  }
  else{
    lambs <- paste("L", 1:(facts*6), sep = "")
  }
  if (data_type == "ordinal"){
    thres <- paste("T", 'i', 1:(facts*6), sep = "")
    thres <- paste(rep(thres, each = 3), letters[1:3])
    ncovs <- max(cumsum(seq(1,1 +4*(facts/2-1), 4)))
    covs <- paste("C", 1:ncovs, sep = "")
    out_vec <- c(lambs, thres, covs)
  } 
  else{
    theta <- paste("Theta", 1:(facts*6), sep = "")
    err_var <- paste("Errvar", 1:(facts*6), sep = "")
    ncovs <- max(cumsum(seq(1,1 +4*(facts/2-1), 4)))
    covs <- paste("C", 1:ncovs, sep = "")
    out_vec <- c(lambs, theta, covs, err_var)
  }
  return(paste(specif, data_type, out_vec, sep = "_"))
}


ColnamesGeneratorSE <- function(specif, data_type, facts){
  if(facts %% 2 != 0){
    stop("Amount of factors should be an even number!")
  }
  if(specif == "MS"){
    x.mis <- sapply(1:(facts/2), function(x) return(6 + 12*(x-1))) + seq(1, facts/2, 1)
    len <- facts*6.5
    lambs <- character(len)
    lambs[-x.mis] <- paste("L_SE", 1:(facts*6), sep = "")
    lambs[x.mis] <- paste("L_SE", letters[1:(facts/2)], sep = "")
  }
  else{
    lambs <- paste("L_SE", 1:(facts*6), sep = "")
  }
  if (data_type == "ordinal"){
    thres <- paste("T_SE", 'i', 1:(facts*6), sep = "")
    thres <- paste(rep(thres, each = 3), letters[1:3])
    ncovs <- max(cumsum(seq(1,1 +4*(facts/2-1), 4)))
    covs <- paste("C_SE", 1:ncovs, sep = "")
    out_vec <- c(lambs, thres, covs)
  }
  else{
    theta <- paste("Theta", 1:(facts*6), sep = "")
    err_var <- paste("Errvar", 1:(facts*6), sep = "")
    ncovs <- max(cumsum(seq(1,1 +4*(facts/2-1), 4)))
    covs <- paste("C_SE", 1:ncovs, sep = "")
    out_vec <- c(lambs, theta, covs, err_var)
  }
  return(paste(data_type, specif, out_vec, sep = "_"))
}

#model building (correctly specified model)
lavaan.data.syn1 <- function(fact=1, nitems=6) {
  TXT <- ""
  for(k in 1:fact) {
    if(k==1){J <- paste0("F", k, " =~ ",
                         paste0("X", (k-1)*nitems + 1:nitems, collapse=" + "))}
    if(k>1){J <- rbind (J, K <- paste0("F", k, " =~ ",
                                       paste0("X", (k-1)*nitems + 1:nitems, collapse=" + ")))}
    TXT <- J
  }
  TXT
}

# misspecified model
lavaan.data.syn2 <- function(fact=1, nitems=6) {
  TXT <- ""
  for(k in 1:fact) {
    if(k==1){J <- paste0("F", k, " =~ ",
                         paste0("X", (k-1)*nitems + 1:nitems, collapse=" + "))}
    if(k==2|k==4|k==6|k==8){ J <- paste0("F", k, " =~ ",
                                         paste0("X", ((k-1)*nitems)-1 + 1:(nitems+1), collapse=" + "))}
    if(k==3|k==5|k==7){J <- paste0("F", k, " =~ ",
                                   paste0("X", (k-1)*nitems + 1:nitems, collapse=" + "))}
    TXT <- rbind(TXT,J)
  }
  TXT
}
################################ Simulation start (1 cell) ##########################

MySimulationCell<- function(Design_mixed = Design_mixed, RowOfDesign, K){
  
  # initialize values
  nvarp <- 6
  fact <- Design_mixed[RowOfDesign,1]
  nvar <- nvarp*fact
  
  ###### ORDINAL ######
  # correctly specified model (ordinal data)
  MyResult_CS_ord_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("CS", "ordinal", fact)))
  
  colnames(MyResult_CS_ord_est) <- ColnamesGeneratorEst("CS", "ordinal", fact)
  
  MyResult_CS_ord_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("CS", "ordinal", fact)))
  colnames(MyResult_CS_ord_err) <- ColnamesGeneratorSE("CS", "ordinal", fact)
  
  MyResult_CS_ord_FI<- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_CS_ord_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "srmr")
  
  # misspecified model (ordinal data)
  MyResult_MS_ord_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("MS", "ordinal", fact)))

  colnames(MyResult_MS_ord_est) <- ColnamesGeneratorEst("MS", "ordinal", fact)
  
  MyResult_MS_ord_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("MS", "ordinal", fact)))
  colnames(MyResult_MS_ord_err) <- ColnamesGeneratorSE("MS", "ordinal", fact)
  
  MyResult_MS_ord_FI <- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_MS_ord_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "srmr")
  
  
  ###### CONTINUOUS ######
  # correctly specified model (continuous data)
  MyResult_CS_con_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("CS", "continuous", fact)))
  colnames(MyResult_CS_con_est) <- ColnamesGeneratorEst("CS", "continuous", fact)
  
  MyResult_CS_con_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("CS", "continuous", fact)))
  
  colnames(MyResult_CS_con_err) <- ColnamesGeneratorSE("CS", "continuous", fact)
  
  MyResult_CS_con_FI<- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_CS_con_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "srmr")
  
  # specified model (continuous data)
  MyResult_MS_con_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("MS", "continuous", fact)))
  colnames(MyResult_MS_con_est) <- ColnamesGeneratorEst("MS", "continuous", fact)
  
  MyResult_MS_con_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("MS", "continuous", fact)))
  colnames(MyResult_MS_con_err) <- ColnamesGeneratorSE("MS", "continuous", fact)
  
  MyResult_MS_con_FI <- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_MS_con_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "srmr")

  #create a loop over the replications k = 1 to K:
  tmp <- proc.time()
  for (k in 1:K){
    # Generate data
    # set a random number seed to be able to replicate the result exactly
    set.seed((k + 1000)*RowOfDesign)
    
    SimDat <- do.call(MyDataGeneration, Design_mixed[RowOfDesign,] )
    
    #Method for ordinal data
    Method_CS_ord <- function(SimData, fact){
      fit_CS_ord <- lavaan:::cfa( model <- lavaan.data.syn1(fact), 
                              data=SimData, std.lv=TRUE, 
                              ordered=c(colnames(SimData)),
                              estimator="PML")
      return(fit_CS_ord)
    }
    
    Method_MS_ord <- function(SimData, fact){
      fit_MS_ord <- cfa( model <- lavaan.data.syn2(fact), 
                     data=SimData, std.lv=TRUE, 
                     ordered=c(colnames(SimData)),
                     estimator="PML")
      return(fit_MS_ord)
    }
    
    #Method for continuous data
    Method_CS_con <- function(SimData, fact){
      fit_CS_CON <- lavaan:::cfa( model <- lavaan.data.syn1(fact), 
                              data=SimData, std.lv=TRUE, 
                              estimator="PML")
      return(fit_CS_CON)
    }
    
    Method_MS_con <- function(SimData, fact){
      fit_MS_CON <- cfa( model <- lavaan.data.syn2(fact), 
                     data=SimData, std.lv=TRUE, 
                     estimator="PML")
      return(fit_MS_CON)
    }
   
    ###################################### ORDINAL #####################################
    ## correctly specified model (ordinal data)
    fit_CS_ord <- Method_CS_ord(SimDat, fact = Design_mixed[RowOfDesign,1])
    
    # parameter estimates
    index <- which(fit_CS_ord@ParTable$free != 0)
    MyAnalysisResult_CS_ord_est <- fit_CS_ord@ParTable$est[index]
    names(MyAnalysisResult_CS_ord_est) <- ColnamesGeneratorEst("CS", "ordinal", fact)
    MyResult_CS_ord_est[k, ] <- MyAnalysisResult_CS_ord_est
    
    # standard errors
    MyAnalysisResult_CS_ord_err <- fit_CS_ord@ParTable$se[index]
    names(MyAnalysisResult_CS_ord_err) <- ColnamesGeneratorSE("CS", "ordinal", fact)
    MyResult_CS_ord_err[k, ] <- MyAnalysisResult_CS_ord_err
    
    ### FITINDICES
    FI_CS_ord <- fitMeasures(fit_CS_ord, 
                             c("chisq.scaled","df.scaled", 
                               "pvalue.scaled", "cfi.scaled",
                               "srmr"))
    MyResult_CS_ord_FI[k,] <- FI_CS_ord
    
    
    #### Misspecified model (ordinal data)
    fit_MS_ord <- Method_MS_ord(SimDat, fact = Design_mixed[RowOfDesign,1])
    
    # parameter estimates
    index <- which(fit_MS_ord@ParTable$free != 0)
    MyAnalysisResult_MS_ord_est <- fit_MS_ord@ParTable$est[index]
    names(MyAnalysisResult_MS_ord_est) <- ColnamesGeneratorEst("MS", "ordinal", fact)
    MyResult_MS_ord_est[k, ] <- MyAnalysisResult_MS_ord_est
    
    # standard errors 
    MyAnalysisResult_MS_ord_err <- fit_MS_ord@ParTable$se[index]
    names(MyAnalysisResult_MS_ord_err) <- ColnamesGeneratorSE("MS", "ordinal", fact)
    MyResult_MS_ord_err[k, ] <- MyAnalysisResult_MS_ord_err
    
    ### FITINDICES
    FI_MS_ord <- fitMeasures(fit_MS_ord, c("chisq.scaled","df.scaled", 
                                       "pvalue.scaled", "cfi.scaled",
                                       "srmr"))
    MyResult_MS_ord_FI[k,] <- FI_MS_ord
    

    ###################################### CONTINUOUS #####################################
    ## correctly specified model (continuous data)
    fit_CS_con <- Method_CS_con(SimDat, fact = Design_mixed[RowOfDesign,1])
    
    # parameter estimates
    index <- which(fit_CS_con@ParTable$free != 0)
    MyAnalysisResult_CS_con_est <- fit_CS_con@ParTable$est[index]
    names(MyAnalysisResult_CS_con_est) <- ColnamesGeneratorEst("CS", "continuous", fact)
    MyResult_CS_con_est[k, ] <- MyAnalysisResult_CS_con_est
    
    # standard errors
    MyAnalysisResult_CS_con_err <- fit_CS_con@ParTable$se[index]
    names(MyAnalysisResult_CS_con_err) <- ColnamesGeneratorSE("CS", "continuous", fact)
    MyResult_CS_con_err[k, ] <- MyAnalysisResult_CS_con_err
    
    ### FITINDICES
    FI_CS_con <- fitMeasures(fit_CS_con, 
                             c("chisq.scaled","df.scaled", 
                               "pvalue.scaled", "cfi.scaled",
                               "srmr"))
    MyResult_CS_con_FI[k,] <- FI_CS_con
    
    # misspecified model (continuous data)
    fit_MS_con <- Method_MS_con(SimDat, fact = Design_mixed[RowOfDesign,1])
    # parameter estimates
    index <- which(fit_MS_con@ParTable$free != 0)
    MyAnalysisResult_MS_con_est <- fit_MS_con@ParTable$est[index]

    names(MyAnalysisResult_MS_con_est) <- ColnamesGeneratorEst("MS", "continuous", fact)
    MyResult_MS_con_est[k, ] <- MyAnalysisResult_MS_con_est
    
    # standard errors
    MyAnalysisResult_MS_con_err <- fit_MS_con@ParTable$se[index]
    names(MyAnalysisResult_MS_con_err) <- ColnamesGeneratorSE("MS", "continuous", fact)
    MyResult_MS_con_err[k, ] <- MyAnalysisResult_MS_con_err
    
    ### FITINDICES
    FI_MS_con <- fitMeasures(fit_MS_con, 
                             c("chisq.scaled", "chisq.scaling.factor", 
                               "df.scaled", "pvalue.scaled", "cfi.scaled",
                               "srmr")) 
    #########################
    MyResult_MS_con_FI[k,] <- FI_MS_con
    
    
    #save the time to run the analyses of K data sets in one cell of the Design_mixed.
    time <- proc.time() - tmp
  }
  # save all relevant results
  return(list(MyResult_CS_ord_est = MyResult_CS_ord_est, 
              MyResult_CS_ord_err = MyResult_CS_ord_err, 
              MyResult_CS_ord_FI = MyResult_CS_ord_FI,
              MyResult_MS_ord_est = MyResult_MS_ord_est,
              MyResult_MS_ord_err = MyResult_MS_ord_err,
              MyResult_MS_ord_FI = MyResult_MS_ord_FI,
              MyResult_CS_con_est = MyResult_CS_con_est,
              MyResult_CS_con_err = MyResult_CS_con_err,
              MyResult_CS_con_FI = MyResult_CS_con_FI,
              MyResult_MS_con_est = MyResult_MS_con_est,
              MyResult_MS_con_err = MyResult_MS_con_err,
              MyResult_MS_con_FI = MyResult_MS_con_FI,
              time = time))
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
  # Save ordinal data results
  MyResult_CS_ord_est <- MyResult[1]
  write.csv(MyResult_CS_ord_est, 
            file = paste("CS_est", "Row", Row,".csv" , sep = ""))
  
  MyResult_CS_ord_err <- MyResult[2]
  write.csv(MyResult_CS_ord_err,
            file = paste("CS_err", "Row", Row,".csv" , sep =""))
  
  MyResult_CS_ord_FI <- MyResult[3]
  write.csv(MyResult_CS_ord_FI,
            file =paste("FI_CS", "Row", Row, ".csv" , sep =""))
  
  # misspecified model (ordinal data)
  MyResult_MS_ord_est <- MyResult[4]
  write.csv(MyResult_MS_ord_est,
            file = paste("MS_est", "Row", Row,".csv" , sep =""))
  
  MyResult_MS_ord_err <- MyResult[5]
  write.csv(MyResult_MS_ord_err,
            file = paste("MS_err", "Row", Row,".csv" , sep =""))
  
  MyResult_MS_ord_FI <- MyResult[6]
  write.csv(MyResult_MS_ord_FI,
            file =paste("FI_MS", "Row", Row, ".csv" , sep =""))
  
  # Save continuous data results
  MyResult_CS_con_est <- MyResult[7]
  write.csv(MyResult_CS_con_est,
            file = paste("CS_est", "Row", Row,".csv" , sep =""))
  
  MyResult_CS_con_err <- MyResult[8]
  write.csv(MyResult_CS_con_err,
            file = paste("CS_err", "Row", Row,".csv" , sep =""))
  
  MyResult_CS_con_FI <- MyResult[9]
  write.csv(MyResult_CS_con_FI,
            file =paste("FI_CS", "Row", Row, ".csv" , sep =""))
  
  # misspecified  model (continuous data)
  MyResult_MS_con_est <- MyResult[10]
  write.csv(MyResult_MS_con_est,
            file = paste("MS_est", "Row", Row,".csv" , sep =""))
  
  MyResult_MS_con_err <- MyResult[11]
  write.csv(MyResult_MS_con_err,
            file = paste("MS_err", "Row", Row,".csv" , sep =""))
  
  MyResult_MS_con_FI <- MyResult[12]
  write.csv(MyResult_MS_con_FI,
            file =paste("FI_MS", "Row", Row, ".csv" , sep =""))
  
  # Save time
  time <- MyResult[13]
  save(time, file =paste("Time", "Row", Row, ".Rdata" , sep =""))
}




#---------------------------------------------------------------------------------
MyDataGeneration <- function(factors, nobs, nvarp = 6, data_type) {
  
  # model specifications
  nvar <- factors*nvarp
  BM <- matrix(c(.8,.8,.7,.7,.6,.6),nrow=1) #factor loadings
  r <- 0.3 # correlation between latent variables
  int <- c(0,0,0,0,0,0) # intercept
  s <- matrix(r, factors, factors) 
  diag(s) <- 1
  
  # intercepts of the items
  int2 <- rep(int,factors)
  
  # b is the loadingsmatrix > transposed matrix version
  b <- t(kronecker(diag(factors), BM)) # kronecker computes product of two arrays 
  
  # add crossloadings (misspecifications)
  x.mis <- sapply(1:(factors/2), function(x) return(6 + 12*(x-1)))
  y.mis <- sapply(1:(factors/2), function(x) return(2 + 2*(x-1)))
  b[x.mis, y.mis] <- .2
  
  #compute error values (theta matrix)
  ev <- diag(1 - diag(b %*% s %*% t(b)))
  
  #compute sigma (variance-covariance matrix) 
  SIGMA <- b %*% s %*% t(b) + ev  #1's on the diagonal
  
  #simulate data from a Multivariate Normal Distribution
  x <- data.frame(MASS:::mvrnorm(n = nobs, mu = int2, Sigma = SIGMA))
  
  item.cutpoints <- 
    matrix(c(rep(-1.20,nvar),rep(0,nvar),
             rep(1.20,nvar)), ncol=nvar, byrow=TRUE)
  
  #item cutpoints (add boundaries)
  item.cutpoints <- rbind(rep(-Inf, ncol(item.cutpoints)), 
                          item.cutpoints, rep(Inf, ncol(item.cutpoints)))
  if (x$X7, x$X8, x$X9, x$X10, x$X11, x$X12){
    for(i in 1:ncol(item.cutpoints)){ 
      x[,i] = cut(x[,i], br=item.cutpoints[,i], 
                  labels=FALSE, include.lowest=TRUE)
    }
  }
  #x is an nobs by nvars matrix with item scores
  return(x)
}

