#All_code: MainSimulationScript (following the code of the manual)

# load packages
library(lavaan)
library(usethis)

args <- commandArgs(TRUE) #?
args <- as.numeric(args)

RowOfDesign <- args[1]
Replication <- args[2]

RowOfDesign <- 1
Replication <- 1

############################# Example small design matrix #############################
factors <- c(2,4)
nobs <- c(50,100)
ncat <- c(2,4)
Design_try <- expand.grid(factors = factors, nobs = nobs, ncat = ncat)

# initial values
n.replications <- 3
factors <- 2
nobs <- 50
ncat <- 2
nvarp <- 6
nvar <- 12

set.seed((Replication + 1000)*RowOfDesign)

##################################### Model building ######################################
# correctly specified model
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


############################ Data generating function ############################
MyDataGeneration <- function(factors, nobs, ncat) {
  
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
  
  # for 2 answer categories (1 threshold -> 0)
  if (ncat==2) {item.cutpoints <- 
    matrix(c(rep(0,nvar)), ncol=nvar, byrow=TRUE)}			
  
  # for 4 answer categories (3 thresholds -> -1.2, 0, 1.2)
  if (ncat==4) {item.cutpoints <- 
    matrix(c(rep(-1.20,nvar),rep(0,nvar),
             rep(1.20,nvar)), ncol=nvar, byrow=TRUE)}
  
  #item cutpoints (add boundaries)
  item.cutpoints <- rbind(rep(-Inf, ncol(item.cutpoints)), 
                          item.cutpoints, rep(Inf, ncol(item.cutpoints)))
  
  for(i in 1:ncol(item.cutpoints)){ 
    x[,i] = cut(x[,i], br=item.cutpoints[,i], 
                labels=FALSE, include.lowest=TRUE)
  }
  
  #x is an nobs by nvars matrix with item scores
  return(x)
}

############################ WLSMV method (old) ############################

# function for the correctly sepcified model 
Method_old_CS <- function(SimData, fact){
  fit1_W <- cfa( model <- lavaan.data.syn1(fact), 
                 data=SimData, std.lv=TRUE, 
                 ordered=c(colnames(SimData)),
                 estimator="WLSMV")
  return(summary(fit1_W, fit.measures = TRUE))
}
# function for the missepcified model 
Method_old_MS <- function(SimData, fact){
  fit2_W <- cfa( model <- lavaan.data.syn2(fact), 
                 data=SimData, std.lv=TRUE, 
                 ordered=c(colnames(SimData)),
                 estimator="WLSMV") 
  return(summary(fit2_W, fit.measures = TRUE))
}

############################ PML method (new) ############################
# function for the correctly sepcified model 
Method_new_CS <- function(SimData, fact){
  fit1_P <- cfa( model <- lavaan.data.syn1(fact), 
                 data=SimData, std.lv=TRUE, 
                 ordered=c(colnames(SimData)),
                 estimator="PML")
  return(summary(fit1_P, fit.measures = TRUE))
  }
  
# function for the missepcified model 
  Method_new_MS <- function(SimData, fact){
  fit2_P <- cfa( model <- lavaan.data.syn2(fact), 
                 data=SimData, std.lv=TRUE, 
                 ordered=c(colnames(SimData)),
                 estimator="PML")
  return(summary(fit2_P, fit.measures = TRUE))
}


# Generate data
SimData <- do.call(MyDataGeneration, Design_try[RowOfDesign, ] ) 

## COMPARE TWO METHODS
tmp <- proc.time()
MyAnalysisResult_WLS1 <- Method_old_CS(SimData, fact = 2)
MyAnalysisResult_WLS2 <- Method_old_MS(SimData, fact = 2)
MyAnalysisResult_PML1 <- Method_new_CS(SimData, fact = 2)
MyAnalysisResult_PML2 <- Method_new_MS(SimData, fact = 2)

time <- proc.time() - tmp

# Save all chi square statistics
MyAnalysisResult <- rbind(WLS_CS= MyAnalysisResult_WLS1$FIT[3], 
                          WLS_MS = MyAnalysisResult_WLS1$FIT[3], 
                          PML_CS = MyAnalysisResult_PML1$FIT[3], 
                          WLS_MS = MyAnalysisResult_PML2$FIT[3])


save(SimData, file = paste("Simdata","Row", 
                          RowOfDesign, "Rep", 
                          Replication ,".Rdata" , 
                          sep =""))
save(MyAnalysisResult, file =paste("Analysis","Row", 
                                   RowOfDesign, "Rep", 
                                   Replication ,".Rdata" , 
                                   sep =""))
save(time, file = paste("Time", "Row", 
                       RowOfDesign, "Rep", 
                       Replication ,".Rdata" , 
                       sep =""))
