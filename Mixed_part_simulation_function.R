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

# test function
MyDataGeneration(factors = 8, nobs = 200)

# Colnames Generator 
ColnamesGeneratorEst <- function(specif, facts){
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
    theta <- paste("Theta", 1:(facts*6), sep = "")
    err_var <- paste("Errvar", 1:(facts*6), sep = "")
    ncovs <- max(cumsum(seq(1,1 +4*(facts/2-1), 4)))
    covs <- paste("C", 1:ncovs, sep = "")
    out_vec <- c(lambs, theta, covs, err_var)
  return(paste(specif, out_vec, sep = "_"))
  }

#test function 
ColnamesGeneratorEst(specif = 'CS', facts = 2)
ColnamesGeneratorEst(specif = 'MS', facts = 2)

ColnamesGeneratorSE <- function(specif, facts){
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
    theta <- paste("Theta", 1:(facts*6), sep = "")
    err_var <- paste("Errvar", 1:(facts*6), sep = "")
    ncovs <- max(cumsum(seq(1,1 +4*(facts/2-1), 4)))
    covs <- paste("C_SE", 1:ncovs, sep = "")
    out_vec <- c(lambs, theta, covs, err_var)
  return(paste(specif, out_vec, sep = "_"))
}

# test function
ColnamesGeneratorSE(specif = 'CS', facts = 2)
ColnamesGeneratorSE(specif = 'MS', facts = 2)


# model building (correctly specified model)
lavaan.data.syn1 <- function(fact=1, nitems=6) {
  TXT <- ""
  for(j in 1:fact) {
    if(j==1){J <- paste0("F", j, " =~ ",
                         paste0("X", (j-1)*nitems + 1:nitems, collapse=" + "))}
    if(j>1){J <- rbind (J, K <- paste0("F", j, " =~ ",
                                       paste0("X", (j-1)*nitems + 1:nitems, collapse=" + ")))}
    TXT <- J
  }
  TXT
}

# model building (misspecified model)
lavaan.data.syn2 <- function(fact=1, nitems=6) {
  TXT <- ""
  for(j in 1:fact) {
    if(j==1){J <- paste0("F", j, " =~ ",
                         paste0("X", (j-1)*nitems + 1:nitems, collapse=" + "))}
    if(j==2|j==4|j==6|j==8){ J <- paste0("F", j, " =~ ",
                                         paste0("X", ((j-1)*nitems)-1 + 1:(nitems+1), collapse=" + "))}
    if(j==3|j==5|j==7){J <- paste0("F", j, " =~ ",
                                   paste0("X", (j-1)*nitems + 1:nitems, collapse=" + "))}
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
  
  # correctly specified model
  MyResult_CS_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("CS", fact)))
  
  colnames(MyResult_CS_est) <- ColnamesGeneratorEst("CS", fact)
  
  MyResult_CS_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("CS", fact)))
  colnames(MyResult_CS_err) <- ColnamesGeneratorSE("CS", fact)
  
  MyResult_CS_FI<- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_CS_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "srmr")
  
  # misspecified model
  MyResult_MS_est <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorEst("MS", fact)))

  colnames(MyResult_MS_est) <- ColnamesGeneratorEst("MS", fact)
  
  MyResult_MS_err <- matrix(NA, nrow = K, ncol = length(ColnamesGeneratorSE("MS", fact)))
  colnames(MyResult_MS_err) <- ColnamesGeneratorSE("MS", fact)
  
  MyResult_MS_FI <- matrix(NA, nrow = K, ncol = 5)
  colnames(MyResult_MS_FI) <- c("chisq.scaled", "df.scaled", 
                                    "pvalue.scaled", "cfi.scaled",
                                    "srmr")
  
  input<- cbind(seq(1,550), rep(0,550))
  mis500<- as.matrix(input)
  
  #create a loop over the replications k = 1 to K:
  tmp <- proc.time()
  for (k in 1:K){
    # Generate data
    # set a random number seed to be able to replicate the result exactly
    set.seed((k + 1000)*RowOfDesign)

    SimDat <- do.call(MyDataGeneration, Design_mixed[RowOfDesign,] )
    
    # fit correctly specified model (with silent check)
    
    fit_CS <- try(cfa( model <- lavaan.data.syn1(fact), 
             data=SimDat, std.lv=TRUE, 
             estimator="PML"),silent=TRUE)
    if(inherits(fit_CS, "try-error")) {
      mis500[k,1] <- 1        
    } else { 

    # parameter estimates
    index <- which(fit_CS@ParTable$free != 0)
    MyAnalysisResult_CS_est <- fit_CS@ParTable$est[index]
    names(MyAnalysisResult_CS_est) <- ColnamesGeneratorEst("CS", fact)
    MyResult_CS_est[k, ] <- MyAnalysisResult_CS_est
    
    # standard errors
    MyAnalysisResult_CS_err <- fit_CS@ParTable$se[index]
    names(MyAnalysisResult_CS_err) <- ColnamesGeneratorSE("CS", fact)
    MyResult_CS_err[k, ] <- MyAnalysisResult_CS_err
    
    ### FITINDICES
    FI_CS <- fitMeasures(fit_CS, 
                             c("chisq.scaled","df.scaled", 
                               "pvalue.scaled", "cfi.scaled",
                               "srmr"))
    MyResult_CS_FI[k,] <- FI_CS
    }
    
    #### Misspecified model
    fit_MS <- try(cfa( model <- lavaan.data.syn2(fact), 
                       data=SimDat, std.lv=TRUE, 
                       estimator="PML"),silent=TRUE)
    if(inherits(fit_MS, "try-error")) {
      mis500[k,2] <- 1        
    } else {
    # parameter estimates
    index <- which(fit_MS@ParTable$free != 0)
    MyAnalysisResult_MS_est <- fit_MS@ParTable$est[index]
    names(MyAnalysisResult_MS_est) <- ColnamesGeneratorEst("MS", fact)
    MyResult_MS_est[k, ] <- MyAnalysisResult_MS_est
    
    # standard errors 
    MyAnalysisResult_MS_err <- fit_MS@ParTable$se[index]
    names(MyAnalysisResult_MS_err) <- ColnamesGeneratorSE("MS", fact)
    MyResult_MS_err[k, ] <- MyAnalysisResult_MS_err
    
    ### FITINDICES
    FI_MS <- fitMeasures(fit_MS, c("chisq.scaled","df.scaled", 
                                       "pvalue.scaled", "cfi.scaled",
                                       "srmr"))
    MyResult_MS_FI[k,] <- FI_MS
    }
    
    #save the time to run the analyses of K data sets in one cell of the Design_mixed.
    time <- proc.time() - tmp
  }
  # save all relevant results
  return(list(MyResult_CS_est = MyResult_CS_est, 
              MyResult_CS_err = MyResult_CS_err, 
              MyResult_CS_FI = MyResult_CS_FI,
              MyResult_MS_est = MyResult_MS_est,
              MyResult_MS_err = MyResult_MS_err,
              MyResult_MS_FI = MyResult_MS_FI,
              Simulated_Data = SimDat,
              time = time, mis500 = mis500))
}

# collect data
MyResult_onecell <- MySimulationCell(Design_mixed, RowOfDesign = 1, K = 2)
MyResult_onecell

################################ Simulation all cells  ###############################
TotalCells <- nrow(Design_mixed)
for (i in 1:TotalCells){
  Row <- i
  MyResult <- MySimulationCell(Design_mixed = Design_mixed, RowOfDesign = Row, K = 2) #10!

  # Write output of one cell of the design
  # Save results
  MyResult_CS_est <- MyResult[1]
  write.csv(MyResult_CS_est, 
            file = paste("CS_est", "Row", Row,".csv" , sep = ""))
  
  MyResult_CS_err <- MyResult[2]
  write.csv(MyResult_CS_err,
            file = paste("CS_err", "Row", Row,".csv" , sep =""))
  
  MyResult_CS_FI <- MyResult[3]
  write.csv(MyResult_CS_FI,
            file =paste("FI_CS", "Row", Row, ".csv" , sep =""))
  
  # misspecified model
  MyResult_MS_est <- MyResult[4]
  write.csv(MyResult_MS_est,
            file = paste("MS_est", "Row", Row,".csv" , sep =""))
  
  MyResult_MS_err <- MyResult[5]
  write.csv(MyResult_MS_err,
            file = paste("MS_err", "Row", Row,".csv" , sep =""))
  
  MyResult_MS_FI <- MyResult[6]
  write.csv(MyResult_MS_FI,
            file =paste("FI_MS", "Row", Row, ".csv" , sep =""))
  
  # Save dataset
  Simulated_data <- MyResult[7]
  write.csv(Simulated_data, 
            file = paste("Simulated_Data", 'Row', Row, ".csv", sep = ""))
  
  # Save time
  time <- MyResult[8]
  save(time, file =paste("Time", "Row", Row, ".Rdata" , sep =""))
  
  # Silent Check
  mis500 <- MyResult[9]
  write.csv(mis500, 
            file = paste("mis500", 'Row', Row, ".csv", sep = ""))
}
