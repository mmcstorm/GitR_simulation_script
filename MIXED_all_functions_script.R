# All functions for simulation (mixed part)

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

  # determine coordinates of cross-loadings (misspecifications)
  x.mis <- sapply(1:(factors/2), function(x) return(6 + 12*(x-1)))
  y.mis <- sapply(1:(factors/2), function(x) return(2 + 2*(x-1)))
  
  xy.mis <- cbind(x.mis, y.mis)
  b[xy.mis] <- .2
  
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
  
  #create indexes for ordinal
  .indexes <- function(fact){
    reps <- fact / 2
    base.vec <- rep(seq(7,12,1), times = reps)
    add.vec <- rep(1:reps - 1, each = 6) * 12
    return(base.vec + add.vec)
  }
  
  for(i in .indexes(factors)){ 
    x[,i] = cut(x[,i], br=item.cutpoints[,i], 
                labels=FALSE, include.lowest=TRUE)
  }
  
  
  #x is an nobs by nvars matrix with item scores
  return(x)
}

# test function


#create indexes for ordinal and continuous
indexes_con <- function(fact){
  reps <- fact / 2
  base.vec <- rep(seq(1,6,1), times = reps)
  add.vec <- rep(1:reps - 1, each = 6) * 12
  return(base.vec + add.vec)
}

indexes_ord <- function(fact){
  reps <- fact / 2
  base.vec <- rep(seq(7,12,1), times = reps)
  add.vec <- rep(1:reps - 1, each = 6) * 12
  return(base.vec + add.vec)
}

ColnamesGeneratorEst <- function(method, specif, facts){
  if(specif == "withC"){
    x.mis <- sapply(1:(facts/2), function(x) return(6 + 12*(x-1))) + seq(1, facts/2, 1)
    len <- facts*6.5
    lambs <- character(len)
    lambs[-x.mis] <- paste("L", 1:(facts*6), sep = "")
    lambs[x.mis] <- paste("L", letters[1:(facts/2)], sep = "")
  }
  else{
    lambs <- paste("L", 1:(facts*6), sep = "")
  }
  intercept <- paste("Int", indexes_con(facts), sep = "")
  thres <- paste("T", 'i', indexes_ord(facts), sep = "")
  thres <- paste(rep(thres, each = 3), letters[1:3])
  err_var <- paste("Errvar", indexes_con(facts), sep = "")
  ncovs <- max(cumsum(seq(1,1 +4*(facts/2-1), 4)))
  covs <- paste("C", 1:ncovs, sep = "")
  out_vec <- c(lambs, thres, err_var, covs, intercept)
  return(paste(method, specif, out_vec, sep = "_"))
}

ColnamesGeneratorSE <- function(method, specif, facts){
  if(specif == "withC"){
    x.mis <- sapply(1:(facts/2), function(x) return(6 + 12*(x-1))) + seq(1, facts/2, 1)
    len <- facts*6.5
    lambs <- character(len)
    lambs[-x.mis] <- paste("L", 1:(facts*6), sep = "")
    lambs[x.mis] <- paste("L", letters[1:(facts/2)], sep = "")
  }
  else{
    lambs <- paste("L_SE", 1:(facts*6), sep = "")
  }
  intercept <- paste("Int_SE", indexes_con(facts), sep = "")
  thres <- paste("T_SE", 'i', indexes_ord(facts), sep = "")
  thres <- paste(rep(thres, each = 3), letters[1:3])
  err_var <- paste("Errvar_SE", indexes_con(facts), sep = "")
  ncovs <- max(cumsum(seq(1,1 +4*(facts/2-1), 4)))
  covs <- paste("C_SE", 1:ncovs, sep = "")
  out_vec <- c(lambs, thres, err_var, covs, intercept)
  return(paste(method, specif, out_vec, sep = "_"))
}

# test function
#ColnamesGeneratorSE(method = "PML", specif = 'withoutC', facts = 4)
#ColnamesGeneratorSE(method = "WLS", specif = 'withC', facts = 2)

# model building (model without specified cross-loadings)
model_withoutC <- function(fact=1, nitems=6) {
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

# model building (model with specified cross-loadings)
model_withC <- function(fact=1, nitems=6) {
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

