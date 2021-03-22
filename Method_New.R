############################ PML method (new) ############################
# function for the correctly specified model 
Method_new_CS <- function(SimData, fact){
  fit1_P <- cfa( model <- lavaan.data.syn1(fact), 
                 data=SimData, std.lv=TRUE, 
                 ordered=c(colnames(SimData)),
                 estimator="PML")
  return(summary(fit1_P, fit.measures = TRUE))
}

# function for the misspecified model 
Method_new_MS <- function(SimData, fact){
  fit2_P <- cfa( model <- lavaan.data.syn2(fact), 
                 data=SimData, std.lv=TRUE, 
                 ordered=c(colnames(SimData)),
                 estimator="PML")
  return(summary(fit2_P, fit.measures = TRUE))
}