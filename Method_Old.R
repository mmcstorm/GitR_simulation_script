############################ WLSMV method (old) ############################
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