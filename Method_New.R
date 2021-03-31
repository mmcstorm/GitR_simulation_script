############################ PML method (new) ############################
# function for the correctly specified model 
Method_new_CS <- function(SimData, fact){
  fit1_P <- lavaan:::cfa( model <- lavaan.data.syn1(fact), 
                 data=SimData, std.lv=TRUE, 
                 ordered=c(colnames(SimData)),
                 estimator="PML")
  return(fit1_P)
}

# function for the misspecified model 
Method_new_MS <- function(SimData, fact){
  fit2_P <- cfa( model <- lavaan.data.syn2(fact), 
                 data=SimData, std.lv=TRUE, 
                 ordered=c(colnames(SimData)),
                 estimator="PML")
  return(fit2_P)
}