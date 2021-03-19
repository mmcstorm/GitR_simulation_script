######### Example small design matrix ########
factors <- c(2,4)
nobs <- c(50,100)
ncat <- c(2,4)
Design_try <- expand.grid(factors = factors, nobs = nobs, ncat = ncat)

n.replications <- 2
factors <- 2
nobs <- 50
ncat <- 2
nvarp <- 6
nvar <- 12
# Adding a comment

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
