
#' Extremal analogue to the innovations algorithm
#'
#' This function estimates TL-MA(q) coefficients using the extremes-analogue
#'    to the innovations algorithm for model of order q = 1, ..., max_q. The
#'    user inputs the estimated (or known) TPDF values and chooses a maximum
#'    order to fit with the recursive algorithm.
#'
#' @param tpdf a `vector` of tpdf values where the first value in the vector is
#'    the lag-0 tpd value.
#' @param max_q the largest order MA(q) that you want to fit. This could be your
#'    desired order or some other large value to see when the estimates
#'    stabilize.
#'
#' @returns a `list` with two components. The first component `coefs`  is a
#'    `matrix` of innovations-estimated coefficients. The `n`th row contains
#'    coefficients for the TL-MA(`n`). The second component `nus` is a
#'    vector of `nu` values from the innovations algorithm.
#' @export
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{mhatre2023innov}{tpdmethods}
#'
#' @examples
#' myTPD <- c(1, 0.4, 0.34, 0.2, 0.11, 0.05, 0.01)
#' out <- innovations(myTPD, max_q = 5)
innovations <- function(tpdf, max_q =50){
  if(is.tpd(tpdf)){
    stopifnot(attr(tpdf, "alpha") == 2)
  } else {
    print("The input is not a `tpd` object.")
    print("... The following assumes alpha=2.")
  }

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
  results$coefs <- theta_hat
  results$nus <- nu
  return(results)
}



#' Transform the marginal distribution to be Fr\'echet(2)
#'
#' This function transforms the marginal distribution of a dataset. It is used on
#'    one margin at a time. It uses a GPD above a quantile `q` and the ECDF below.
#'
#' @param x univariate data (`vector`) to be transformed
#' @param q the quantile to use as a cutoff between the ECDF and GPD components
#' @param use_gpd change to `FALSE` if you want to use a rank-transformation
#'    only. When `TRUE` a GPD is used to estimate the upper tail. (default is
#'    `TRUE`)
#' @param gpd_scale fix the scale at your desired scale value. If no value is
#'    input then it will be estimated. (default is -99)
#' @param gpd_shape fix the shape at your desired shape value. If no value is
#'    input then it will be estimated. Note that gpd shape is 1/alpha so if you
#'    have generated data with alpha=2 you should enter 1/2. (default is -99)
#'
#' @returns a `vector` of data on Fr\'echet(2) margins
#' @export
#'
#' @importFrom stats ecdf quantile
#'
#' @examples
#' myData <- rnorm(1000)
#' out <- transform_marginal(myData, q = 0.95)
transform_marginal <- function(x,
                               q = 0.975,
                               use_gpd = TRUE,
                               fix_gpd_params = FALSE,
                               gpd_scale = -99,
                               gpd_shape = -99){
  quant_q <- unname(quantile(prob = q, x))
  large_x <- which(x >= quant_q)
  x_new   <- ecdf(x)(x)
  x_new   <- x_new - 0.5 * min(x_new)
  if(use_gpd == TRUE){
    if(gpd_scale != -99 && gpd_shape != -99){
      for(i in 1:length(x)){
        if(i %in% large_x){
          x_new[i] <- q + (1-q)*evd::pgpd(x[i],
                                          loc = quant_q,
                                          scale = gpd_scale,
                                          shape = gpd_shape)
        }
      }
    } else if(gpd_scale != -99){
      params <- evd::fpot(x,
                          threshold = quant_q,
                          scale = gpd_scale,
                          std.err = F)$estimate
      for(i in 1:length(x)){
        if(i %in% large_x){
          x_new[i] <- q + (1-q)*evd::pgpd(x[i],
                                          loc = quant_q,
                                          scale = gpd_scale,
                                          shape = params[1])
        }
      }
    } else if(gpd_shape != -99){
      params <- evd::fpot(x,
                          threshold = quant_q,
                          shape = gpd_shape,
                          std.err = F)$estimate
      for(i in 1:length(x)){
        if(i %in% large_x){
          x_new[i] <- q + (1-q)*evd::pgpd(x[i],
                                          loc = quant_q,
                                          scale = params[1],
                                          shape = gpd_shape)
        }
      }
    } else {
      params <- evd::fpot(x,
                          threshold = quant_q,
                          std.err = F)$estimate
      for(i in 1:length(x)){
        if(i %in% large_x){
      x_new[i] <- q + (1-q)*evd::pgpd(x[i],
                                      loc = quant_q,
                                      scale = params[1],
                                      shape = params[2])
        }
      }
    }
  }

  return(evd::qfrechet(x_new, shape = 2))
}


#' Generate a TL-AR(1) time series
#'
#' This function generates a TL-AR(1) time series of length `n`.
#'
#' @param n the length of the desired time series
#' @param phi the AR(1) parameter
#'
#' @returns a `vector`. The length-`n` times series with AR(1) parameter `phi`.
#' @export
#'
#' @references
#' \insertRef{mhatre2024arma}{tpdmethods}
#'
#' @examples
#' out <- gen_ar1(1000, phi = 0.2)
gen_ar1 <- function(n, phi){
  RVnoise   <- evd::rfrechet(1000+n, shape = 2)
  ar1_ts    <- numeric(1000 + n)
  ar1_ts[1] <- RVnoise[1]
  for(i in 2:(1000+n)){
    ar1_ts[i] <- tadd(tmult(phi, ar1_ts[i-1]), RVnoise[i])
  }
  ar1_ts <- ar1_ts[1001:(1000+n)]
  transform_marginal(ar1_ts)
}



#' Generate a TL-ARMA(1, 1) time series
#'
#' This function generates a TL-ARMA(1,1) time series of length `n`.
#'
#' @param n the length of the desired time series
#' @param phi the AR(1) parameter
#' @param theta the MA(1) parameter
#'
#' @returns a `vector`. The length-`n` times series with AR(1) parameter `phi`
#'    and MA(1) parameter `theta`.
#' @export
#'
#'@references
#' \insertRef{mhatre2024arma}{tpdmethods}
#'
#' @examples
#' out <- gen_arma11(1000, phi = 0.3, theta = -0.1)
gen_arma11 <- function(n, phi, theta){
  RVnoise   <- evd::rfrechet(1000+n, shape = 2)
  arma11_ts    <- numeric(1000 + n)
  arma11_ts[1] <- RVnoise[1]
  for(i in 2:(1000+n)){
    arma11_ts[i] <- f(phi * finv(arma11_ts[i-1]) +
                        finv(RVnoise[i]) +
                        theta * finv(RVnoise[i-1]))
  }
  arma11_ts <- arma11_ts[1001:(1000+n)]
  transform_marginal(arma11_ts)
}


#' Generate a TL-MA(`q`) time series
#'
#' This function generates a TL-MA(`q`) time series of length `n` where `q`
#'    is determined by `length(thetas)`.
#'
#' @param n the length of the desired time series
#' @param thetas the `vector` of MA(`q`) parameters
#'
#' @returns a `vector`. The length-`n` times series with MA(`q`) parameters
#'    `thetas`.
#' @export
#'
#'@references
#' \insertRef{mhatre2024arma}{tpdmethods}
#'
#' @examples
#' out <- gen_maq(n = 1200, thetas = c(0.8, 0.2, 0.3, 0.1))
gen_maq <- function(n, thetas){
  q         <- length(thetas)
  RVnoise   <- evd::rfrechet(q+n+1, shape = 2)
  maq_ts    <- numeric(n)
  temp_vals <- numeric(q+1)
  for(i in (q+1):(n+q)){
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


#' Generate a TL-ARMA(p, q) time series
#'
#' This function generates a TL-ARMA(p, q) time series of length `n`.
#'
#' @param n the length of the desired time series
#' @param phi a p-`vector` of the auto-regressive (AR) parameters where p is the
#'    order of the AR component of the TL-ARMA(p,q) process.
#' @param theta a q-`vector` of the moving average (MA) parameters where q is the
#'    order of the MA component of the TL-ARMA(p,q) process.
#'
#' @returns a `vector`. The length-`n` times series with AR parameters `phi`
#'    and MA parameters `theta`.
#' @export
#'
#' @references
#' \insertRef{mhatre2024arma}{tpdmethods}
#'
#' @examples
#' out <- gen_arma(1000, phi = c(0.3), theta = c(-0.1, 0.4, 0.2))
gen_arma <- function(n, phi, theta){
  p <- length(phi)
  q <- length(theta)
  if(p != 0){ # check if invertible
    if(root_check(phi) < 1){
      stop("#!#!# The time series is not generated because the AR polynomial
            has roots in the unit circle (the process is not invertible).#!#!#")
    }
  }
  if(p == 0){ # MA(q)
    out <- gen_maq(n, thetas = theta)
  } else if(p == 1 && q == 0){ # AR(1)
    out <- gen_ar1(n, phi = phi)
  } else if(p == 1 && q == 1){ # ARMA(1,1)
    out <- gen_arma11(n, phi = phi, theta = theta)
  } else if(q == 0){ # AR(p)
    RVnoise <- evd::rfrechet(1000*p+n, shape = 2)
    out     <- numeric(1000*p + n)
    out[1:(p+q)]  <- RVnoise[1:(p+q)]
    for(i in (p+q+1):(1000*p+n)){
      temp_ar <- 0
      for(a in 1:p){
        temp_ar <- temp_ar + phi[a] * finv(out[i-a])
      }
      out[i] <- f(temp_ar + finv(RVnoise[i]))
    }
  } else { # general p and q
    RVnoise <- evd::rfrechet(1000*p+n, shape = 2)
    out     <- numeric(1000*p + n)
    out[1:(p+q)]  <- RVnoise[1:(p+q)]
    for(i in (p+q+1):(1000*p+n)){
      temp_ar <- 0
      for(a in 1:p){
        temp_ar <- temp_ar + phi[a] * finv(out[i-a])
      }
      temp_ma <- finv(RVnoise[i])
      for(m in 1:q){
        temp_ma <- temp_ma + theta[m] * finv(RVnoise[i-m])
      }
      out[i] <- f(temp_ar + temp_ma)
    }
  }
  out <- out[(1000*p+1):(1000*p+n)]
  transform_marginal(out)
}




#' Returns the modulus of the root of a polynomial that is closest to the origin.
#'
#' @param phi the parameter vector of the AR polynomial to be checked.
#'
#' @returns the scalar modulus of the smallest root
#' @export
#'
#' @examples
#' root_check(phi = c(0.8, 0.11, 0.1))
#' root_check(phi = c(0.8, 0.01, 0.1))
root_check <- function(phi){
  min(Mod(polyroot(c(1, -rev(phi)))))
}




#' Compute the TPD function from TL-MA(q) parameters
#'
#' This function takes a parameter vector `thetas` as input and returns the
#'    model TPDF for lags 0 through `max_lag`
#'
#' @param thetas the q-`vector` of parameters from a TL-MA(q) model
#' @param max_lag the TPDF is computed for lags 1 through `max_lag`
#'    (default is 20)
#'
#' @returns a `vector` of model TPDF values
#' @export
#'
#' @references
#' \insertRef{mhatre2024arma}{tpdmethods}
#'
#' @examples
#' out <- maq_tpdf(c(0.8, 0.4, 0.01), max_lag = 5)
maq_tpdf <- function(thetas, max_lag = 20){
  q <- length(thetas)
  sigmas <- rep(0, max_lag)
  thetas <- c(1, thetas, rep(0, max_lag))
  for(i in 1:q){
    sigmas[i] <- sum(pmax(thetas[1:(q+1)], rep(0, q+1)) *
                       pmax(thetas[(i+1):(i+q+1)], rep(0, q+1))) /
      sum(thetas^2) # this is sigma_0 and thus forces the scale to be 1
  }
  return(validate_tpd(new_tpd(sigmas, alpha = 2)))
}


#' Compute the TPD function from a TL-AR(1) parameter
#'
#' This function takes a parameter scalar `phi` as input and returns the
#'    model TPDF for lags 1 through `max_lag`
#'
#' @param phi the TL-AR(1) parameter
#' @param max_lag the TPDF is computed for lags 0 through`max_lag`
#'
#' @returns a `vector` of model TPDF values
#' @export
#'
#' @references
#' \insertRef{mhatre2024arma}{tpdmethods}
#'
#' @examples
#' out <- ar1_tpdf(0.3)
ar1_tpdf <- function(phi, max_lag = 20){
  sigmas <- pmax(phi^(1:max_lag), rep(0, max_lag)) # scale 1 implies just phi^h
  return(validate_tpd(new_tpd(sigmas, alpha = 2)))
}


#' Compute the TPD value for the lag-h pair from a TL-ARMA(1,1) model
#'
#' This function takes a `vector` pair of parameter values as input and outputs
#'    the model TPDF for a single lag `h`. This function is called by
#'    `arma11_tpdf` for each lag.
#'
#'
#' @param params a vector of the form (`theta`, `phi`) where `theta` is the
#'    TL-MA coefficient and `phi` is the TL-AR coefficient
#' @param h the lag for the desired TPD value
#'
#' @returns a scalar TPD value
#' @export
#'
#' @references
#' \insertRef{mhatre2024arma}{tpdmethods}
#'
#' @examples
#' out <- get_arma11_h(c(0.4, 0.1), h = 3)
get_arma11_h <- function(params, h){
  theta <- params[1]
  phi <- params[2]
  out <- 0
  if(phi > 0){
    if(phi + theta > 0){
      num <- (phi + theta) * phi^(h-1) * (1 + phi*theta)
      den <- 1 + 2 * phi * theta + theta^2
      out <- num / den
    }
  } else if(phi + theta > 0){
    if(h %% 2 == 0){
      num <- (phi + theta)^2 * phi ^ h
      den <- 1 - phi^4 + (phi + theta)^2
      out <- num / den
    } else {
      num <- (phi + theta) * phi^(h - 1) * (1 - phi^4)
      den <- 1 - phi^4 + (phi + theta)^2
      out <- num / den
    }
  } else if(phi + theta < 0){
    if(h %% 2 == 0){
      num <- (phi + theta) * phi^(h - 1) * (1 + theta * phi^3)
      den <- 1 + phi^2 * theta^2 + 2 * phi^3 * theta
      out <- num / den
    }
  }
  return(out)
}

#' Compute the TPD function from a TL-ARMA(1,1) parameter `vector`
#'
#' This function takes a parameter `vector` as input and returns the
#'    model TPDF for lags 1 through `max_lag`
#'
#' @param params a vector of the form (`theta`, `phi`) where `theta` is the
#'    TL-MA coefficient and `phi` is the TL-AR coefficient
#' @param max_lag the TPDF is computed for lags 0 through`max_lag`
#'
#' @returns a `vector` of model TPDF values
#' @export
#'
#' @references
#' \insertRef{mhatre2024arma}{tpdmethods}
#'
#' @examples
#' out <- arma11_tpdf(c(0.9, 0.1), max_lag = 15)
arma11_tpdf <- function(params, max_lag = 20){
  sigmas <- rep(0, max_lag)
  for(h in 1:max_lag){
    sigmas[h] <- get_arma11_h(params, h)
  }
  return(validate_tpd(new_tpd(sigmas, alpha = 2)))
}
