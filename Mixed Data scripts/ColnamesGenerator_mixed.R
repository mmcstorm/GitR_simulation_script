# Colnames Generator 
ColnamesGeneratorEst <- function(method, specif, facts, ncat){
  if(facts %% 2 != 0){stop("Amount of factors should be an even number!")}
  if(specif == "MS"){
    x.mis <- sapply(1:(facts/2), function(x) return(6 + 12*(x-1))) + seq(1, facts/2, 1)
    len <- facts*6.5
    lambs <- character(len)
    lambs[-x.mis] <- paste("L", 1:(facts*6), sep = "")
    lambs[x.mis] <- paste("L", letters[1:(facts/2)], sep = "")
  }
  else{
    lambs <- paste("L", 1:(facts*6), sep = "")
  }
  if (ncat == 4){
    thres <- paste("T", 'i', 1:(facts*6), sep = "")
    thres <- paste(rep(thres, each = 3), letters[1:3])
  }
  else{
    thres <- paste("T", 'i', 1:(facts*6), sep = "")
  }
  ncovs <- max(cumsum(seq(1,1 +4*(facts/2-1), 4)))
  covs <- paste("C", 1:ncovs, sep = "")
  out_vec <- c(lambs, thres, covs)
  return(paste(method, specif, out_vec, sep = "_"))
}

ColnamesGeneratorSE <- function(method, specif, facts, ncat){
  if(facts %% 2 != 0){stop("Amount of factors should be an even number!")}
  if(specif == "MS"){
    x.mis <- sapply(1:(facts/2), function(x) return(6 + 12*(x-1))) + seq(1, facts/2, 1)
    len <- facts*6.5
    lambs <- character(len)
    lambs[-x.mis] <- paste("L_SE", 1:(facts*6), sep = "")
    lambs[x.mis] <- paste("L_SE", letters[1:(facts/2)], sep = "")
  }
  else{
    lambs <- paste("L_SE", 1:(facts*6), sep = "")
  }
  if (ncat == 4){
    thres <- paste("T_SE", 'i', 1:(facts*6), sep = "")
    thres <- paste(rep(thres, each = 3), letters[1:3])
  }
  else{
    thres <- thres <- paste("T_SE", 'i', 1:(facts*6), sep = "")
  }
  ncovs <- max(cumsum(seq(1,1 +4*(facts/2-1), 4)))
  covs <- paste("C_SE", 1:ncovs, sep = "")
  out_vec <- c(lambs, thres, covs)
  return(paste(method, specif, out_vec, sep = "_"))
}