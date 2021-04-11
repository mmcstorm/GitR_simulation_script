#MainSimulationScript (following the code of the manual)

############################# Example small design matrix #############################
factors <- c(2,4)
nobs <- c(50,100)
ncat <- c(2,4)

##Create the simulation design matrix (full factorial)
# Design is a data.frame with all possible combinations of the factor levels
# Each row of the design matrix represents a cell of your simulation design

Design_small <- expand.grid(factors = factors, nobs = nobs, ncat = ncat)

# initial values
nvarp <- 6 

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

#load packages
library(lavaan)
library(usethis)

# prepared functions
source("MyDataGeneration.R")
source("Method_new.R")
source("Method_old.R")

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
MyResult_onecell <- MySimulationCell(Design_small, RowOfDesign = 1, K = 1)
colnames(MyResult_onecell) <- c('WLS_Correct', 'WLS_Misspec', 'PML_Correct', 'PML_Misspec')
MyResult_onecell

# Write output of one cell of the design
#save(MyResult, file =paste("MyResult", "Row", Row,".Rdata" , sep =""))
#optional to save timing of analyses of K replications in 1 cell
#save(time, file =paste("Time", "Row", Row, ".Rdata" , sep =""))


################################ Simulation all cells  ###############################
TotalCells <- nrow(Design_small)
for (i in 1:TotalCells){
  Row <- i
  MyResult <- MySimulationCell(Design = Design_small, RowOfDesign = Row, K = 1 ) #10!
  # Write output of one cell of the design
  save(MyResult, file =paste("MyResult", "Row", Row,".Rdata" , sep =""))
  save(time, file =paste("Time", "Row", Row, ".Rdata" , sep =""))
}
