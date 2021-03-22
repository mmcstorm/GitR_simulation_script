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
Design_small <- expand.grid(factors = factors, nobs = nobs, ncat = ncat)

############################# Official design matrix #############################
#factors <- c(4,6,8) 					          	   #number of latent variables
#nobs <- c(200,400,800)                      #sample size
#ncat <- c(2,4)                              #number of categories

##Create the simulation design matrix (full factorial)
# Design is a data.frame with all possible combinations of the factor levels
# Each row of the design matrix represents a cell of your simulation design
# Design <- expand.grid(factors = factors, nobs = nobs, ncat = ncat)

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
  return(summary(fit1_W, fit.measures = TRUE, nd = 8))
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
# function for the correctly specified model 
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
SimData <- do.call(MyDataGeneration, Design_small[RowOfDesign, ] )

## COMPARE TWO METHODS
tmp <- proc.time()
MyAnalysisResult_WLS1 <- Method_old_CS(SimData, fact = Design_small[RowOfDesign,1])
MyAnalysisResult_WLS2 <- Method_old_MS(SimData, fact = Design_small[RowOfDesign,1])
MyAnalysisResult_PML1 <- Method_new_CS(SimData, fact = Design_small[RowOfDesign,1])
MyAnalysisResult_PML2 <- Method_new_MS(SimData, fact = Design_small[RowOfDesign,1])

time <- proc.time() - tmp

# Save all chi square statistics
MyAnalysisResult <- rbind(WLS_CS=  MyAnalysisResult_WLS1$FIT[3], 
                          WLS_MS = MyAnalysisResult_WLS2$FIT[3], 
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

################################ Simulation start (1 cell) ##########################
MySimulationCell<- function(Design = Design_small, RowOfDesign = 2, K = 2){
  # Input arguments:
  # Design = designmatrix
  # RowOfDesign: number that refers to the row of the design matrix = one cell
  # K: Total number of replications = number of data sets generated in one cell
  # Create matrix or dataframe to store the results:
  MyResult <- matrix(NA, nrow = K, ncol = 4)
 
  ### Fit indices
  ### Parameter estimates
  
  #create a loop over the replications k = 1 to K:
  tmp <- proc.time()
  for (k in 1:K){
    # Generate data
    # set a random number seed to be able to replicate the result exactly
    set.seed((k + 1000)*RowOfDesign)
    SimDat <- do.call(MyDataGeneration, Design_small[RowOfDesign,] )

    # Analyze data set with Method_new
    MyAnalysisResult_WLS1 <- Method_old_CS(SimDat, fact = Design_small[RowOfDesign,1])
    MyAnalysisResult_WLS2 <- Method_old_MS(SimDat, fact = Design_small[RowOfDesign,1])
    MyAnalysisResult_PML1 <- Method_new_CS(SimDat, fact = Design_small[RowOfDesign,1])
    MyAnalysisResult_PML2 <- Method_new_MS(SimDat, fact = Design_small[RowOfDesign,1])
    
    MyAnalysisResult <- cbind(WLS_CS = MyAnalysisResult_WLS1$FIT[3], 
                              WLS_MS = MyAnalysisResult_WLS2$FIT[3], 
                              PML_CS = MyAnalysisResult_PML1$FIT[3], 
                              WLS_MS = MyAnalysisResult_PML2$FIT[3])
    #Evaluate the analysis results of Method_new (Result1) and Mehtod_old (Result2)
    #MyResult1 <- MyEvaluationPC(MyAnalysisResult1)
    #MyResult2 <- MyEvaluationPC(MyAnalysisResult2)
    #store the results in the right row k of your result matrix:
    #We only store the second result which is the evaluation criterion
    MyResult[k, ] <- MyAnalysisResult
    #colnames(MyResult) <- colnames(MyAnalysisResult)
    #rownames(MyResult) <- rownames(MyAnalysisResult)
  }
  #save the time to run the analyses of K data sets in one cell of the design.
  time <- proc.time() - tmp
  return(MyResult)
}

# collect data
#Row <- 1
MyResult_onecell <- MySimulationCell(Design_small, RowOfDesign = 2, K = 50)
colnames(MyResult_onecell) <- c('WLS_Correct', 'WLS_Misspec', 'PML_Correct', 'PML_Misspec')

plot
#MySimulationCell(Design = Design_small, RowOfDesign = 1, K = 1)
#MySimulationCell(Design = Design_small, RowOfDesign = 2, K = 2)

# ERROR: The variance-covariance matrix of the estimated parameters (vcov)
#does not appear to be positive definite! The smallest eigenvalue
#(= -2.928099e-17) is smaller than zero. This may be a symptom that
#the model is not identified.


# Write output of one cell of the design
#save(MyResult, file =paste("MyResult", "Row", Row,".Rdata" , sep =""))
#optional to save timing of analyses of K replications in 1 cell
#save(time, file =paste("Time", "Row", Row, ".Rdata" , sep =""))


################################ Simulation all cells  ###############################
TotalCells <- nrow(Design_small)
for (i in 1:TotalCells){
  Row <- i
  MyResult <- MySimulationCell(Design = Design_small, RowOfDesign = Row, K = 5 ) #10!
  # Write output of one cell of the design
  save(MyResult, file =paste("MyResult", "Row", Row,".Rdata" , sep =""))
  ##write.csv(MyResult, file =paste("MyResult", "Row", Row,".csv" , sep =""))
}
warnings()
