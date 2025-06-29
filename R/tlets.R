
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
#'    one margin at a time. It uses a GPD above a quantile q and the ECDF below.
#'
#' @param x univariate data (`vector`) to be transformed
#' @param q the quantile to use as a cutoff between the ECDF and GPD components
#' @param use_gpd change to `FALSE` if you want to use a rank-transformation
#'    only. When `TRUE` a GPD is used to estimate the upper tail. (default is
#'    `TRUE`)
#' @param fix_gpd_params change to `TRUE` if you want to fix either the scale
#'    or the shape of the GPD used to fit the upper tail of the marginal
#'    distribution. (default is `FALSE`)
#' @param gpd_scale set at your desired scale value. Only used if
#'    `fix_gpd_params = TRUE`. (default is -99)
#' @param gpd_shape set at your desired shape value. Only used if
#'    `fix_gpd_params = TRUE`. Note that gpd shape is 1/alpha. (default is -99)
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
    if(fix_gpd_params == TRUE){
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
  a <- transform_marginal(arma11_ts)
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
  return(sigmas)
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
  pmax(phi^(1:max_lag), rep(0, max_lag)) # scale 1 implies just phi^h
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
  temp_out <- rep(0, max_lag)
  for(h in 1:max_lag){
    temp_out[h] <- get_arma11_h(params, h)
  }
  return(temp_out)
}
