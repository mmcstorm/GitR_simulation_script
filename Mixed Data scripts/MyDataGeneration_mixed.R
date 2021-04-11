# My Data Generation function 

MyDataGeneration <- function(factors, nobs, ncat, nvarp = 6) {
  
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