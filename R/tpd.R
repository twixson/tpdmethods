TPDF <- function(ts, thresh = 0.975, maxlag = 50){
  tpdf <- rep(0, maxlag+1)
  if(class(ts) == "list"){
    tpdf[1] <- 1
    for(i in 2:(length(ts)+1)){
      temp_mat <- ts[[i-1]]
      temp_r   <- rowSums(temp_mat^2)
      tpdf[i]  <- 2*mean(temp_mat[, 1]*temp_mat[, 2]/temp_r)
    }
  } else {
    rt0 <- sqrt(ts^2 + ts^2)
    r00 <- sort(rt0)[floor(thresh*(length(rt0)+1))]
    tpdf <- rep(0, maxlag+1)
    tpdf[1] <- 2/length(which(rt0>r00)) * sum((ts^2)[which(rt0>r00)] /
                                                (rt0[which(rt0>r00)]^2))
    n <- length(ts)
    rt <- matrix(NA, nrow = n-1, ncol = maxlag)
    r0 <- rep(NA, maxlag)
    for(i in 1:maxlag){
      temp_means <- c(mean(ts[1:(n-i)]), mean(ts[(1+i):n]))
      temp_mat <- matrix(
        c(pmax(rep(0, (n-i)), ts[1:(n-i)] - temp_means[1]),
          pmax(rep(0, (n-i)), ts[(1+i):n] - temp_means[2])),
        ncol = 2, byrow = F)
      rt[,i] <- c(sqrt(temp_mat[,1]^2 + temp_mat[,2]^2), rep(NA, i-1))
      temp <- na.omit(rt[,i])
      if(length(temp)){ # ignore zero-length columns
        k <- floor(thresh*(length(temp)+1))
        r0[i] <- sort(temp)[k]
      } else {r0[i] <- 0}
      indices <- which(rt[,i]>r0[i])
      tpdf[i+1] <- 2/length(indices) * sum((temp_mat[,1]*temp_mat[,2])[indices] /
                                             (rt[indices, i]^2))
    }
  }
  return(tpdf)
}
