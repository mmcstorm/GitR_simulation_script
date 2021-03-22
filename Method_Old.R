############################ WLSMV method (old) ############################

# function for the correctly specified model 
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