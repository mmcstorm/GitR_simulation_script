# MainSimulationScriptSLURM (following the code of the manual)
args <- commandArgs(TRUE) #SLURM command
args <- as.numeric(args)

RowOfDesign <- args[1]
Replication <- args[2]

RowOfDesign <- 24
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

#load packages
library(lavaan)
library(usethis)

# prepared functions
source("MyDataGeneration.R")
source("Method_new.R")
source("Method_old.R")

################################ Simulation start (1 cell) ##########################

  
  ### Fit indices
  ### Parameter estimates
  
  #create a loop over the replications k = 1 to K:
  tmp <- proc.time()

    # Generate data
    # set a random number seed to be able to replicate the result exactly
    set.seed((Replication + 1000)*RowOfDesign)
    SimDat <- do.call(MyDataGeneration, Design[RowOfDesign,] )
    
    # Analyze data set with Method_new
    MyAnalysisResult_WLS1 <- Method_old_CS(SimDat, fact = Design[RowOfDesign,1])
    MyAnalysisResult_WLS2 <- Method_old_MS(SimDat, fact = Design[RowOfDesign,1])
    MyAnalysisResult_PML1 <- Method_new_CS(SimDat, fact = Design[RowOfDesign,1])
    MyAnalysisResult_PML2 <- Method_new_MS(SimDat, fact = Design[RowOfDesign,1])
    
    MyAnalysisResult <- cbind(WLS_CS = MyAnalysisResult_WLS1$FIT[1:3], 
                              WLS_MS = MyAnalysisResult_WLS2$FIT[3], 
                              PML_CS = MyAnalysisResult_PML1$FIT[3], 
                              PML_MS = MyAnalysisResult_PML2$FIT[3])
    MyAnalysisResult <- cbind()
    #Evaluate the analysis results of Method_new (Result1) and Mehtod_old (Result2)
    #MyResult1 <- MyEvaluationPC(MyAnalysisResult1)
    #MyResult2 <- MyEvaluationPC(MyAnalysisResult2)
    #store the results in the right row k of your result matrix:
    #We only store the second result which is the evaluation criterion

  #save the time to run the analyses of K data sets in one cell of the design.
  time <- proc.time() - tmp
  

# collect data
colnames(MyAnalysisResult) <- c('WLS_Correct', 'WLS_Misspec', 'PML_Correct', 'PML_Misspec')
MyResult_onecell

# Write output of one cell of the design
#save(MyResult, file =paste("MyResult", "Row", Row,".Rdata" , sep =""))
#optional to save timing of analyses of K replications in 1 cell
#save(time, file =paste("Time", "Row", Row, ".Rdata" , sep =""))


################################ Simulation all cells  ###############################
setwd("/exports/fsw/mmcstorm/SimData")
# create folder data
save(SimDat, file =paste("Data", "Row", RowOfDesign, "Rep", Replication ,".Rdata" , sep =""))

setwd("/exports/fsw/mmcstorm/Analysis")
save(MyAnalysisResult, file =paste("Analysis", "Row", RowOfDesign, "Rep", Replication ,".Rdata" , sep =""))

