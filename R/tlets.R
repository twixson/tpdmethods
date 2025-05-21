innovations <- function(tpdf, max_q =50){
  #initialize variables nu and theta_hat
  nu <- rep(NA, max_q + 1)
  theta_hat <- matrix(0, nrow = max_q, ncol = max_q)
  #set nu_0 equal to tpdf(0)
  nu[1] <- tpdf[1]
  #compute n theta_hat's and nu_n up to n=max_q iterations
  for(n in 1:max_q){
    for(k in 0:(n-1)){
      if(k==0){
        temp <- tpdf[n+1]
      } else {temp <- 0
      for(j in 0:(k-1)){
        temp <- temp + theta_hat[k, k-j]*theta_hat[n, n-j]*nu[j+1]}
      temp <- tpdf[n-k+1] - temp}
      theta_hat[n, n-k] <- nu[k+1]^(-1)*temp}
    nu[n+1] <- tpdf[1]
    for(l in 0:(n-1)){nu[n+1] <- nu[n+1] - theta_hat[n, n-l]^2*nu[l+1]}}
  results <- list()
  results[[1]] <- theta_hat
  results[[2]] <- nu
  return(results)
}


transform_marginal <- function(x, q = 0.975){
  quant_q <- unname(quantile(prob = q, x))
  large_x <- which(x >= quant_q)
  params  <- fpot(x, threshold = quant_q, shape = 1/2, std.err = F)$estimate
  x_new   <- ecdf(x)(x)
  for(i in 1:length(x)){
    if(i %in% large_x){
      x_new[i] <- q + (1-q)*pgpd(x[i], shape = 1/2, loc = quant_q, scale = params)
    }
  }

  return(qfrechet(x_new, shape = 2))
}


gen_ar1 <- function(n, phi){
  RVnoise   <- rfrechet(1000+n, shape = 2)
  ar1_ts    <- numeric(1000 + n)
  ar1_ts[1] <- RVnoise[1]
  for(i in 2:(1000+n)){
    ar1_ts[i] <- tadd(tmult(phi, ar1_ts[i-1]), RVnoise[i])
  }
  ar1_ts <- ar1_ts[1001:(1000+n)]
  transform_marginal(ar1_ts)
}


gen_arma11 <- function(n, phi, theta){
  RVnoise   <- rfrechet(1000+n, shape = 2)
  arma11_ts    <- numeric(1000 + n)
  arma11_ts[1] <- RVnoise[1]
  for(i in 2:(1000+n)){
    arma11_ts[i] <- f(phi * finv(arma11_ts[i-1]) +
                        finv(RVnoise[i]) +
                        theta * finv(RVnoise[i-1]))
  }
  arma11_ts <- arma11_ts[1001:(1000+n)]
  a <- transform_marginal(arma11_ts)
}


gen_maq <- function(n, thetas){
  q         <- length(thetas)
  RVnoise   <- rfrechet(q+n+1, shape = 2)
  maq_ts    <- numeric(n)
  temp_vals <- numeric(q+1)
  for(i in (q+1):(n+q+1)){
    temp_vals[1] <- RVnoise[i]
    for(j in 1:q){
      temp_vals[j+1] <- tmult(thetas[j], RVnoise[i-j])
    }
    temp_val <- temp_vals[1]
    for(j in 1:q){
      temp_val <- tadd(temp_val, temp_vals[j+1])
    }
    maq_ts[i-q] <- temp_val
  }
  transform_marginal(maq_ts)
}
