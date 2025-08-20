
#' Compute the Tail Pairwise Dependence (TPD) of an `R` object
#'
#' This function puts the users input into a form recognizable by simpler
#'    functions which transform the marginal distribution, extract the large
#'    points, and compute the TPD value. It currently recognizes `matrices`,
#'    `vectors` as time series, and `lists` of bivariate data. `Matrices` can
#'    be bivariate data, multivariate data, or time series data split into
#'    seasons (use `matrix_as_seasons` = `TRUE`).
#'
#' @param data as a `matrix`, `vector`, or `list`
#' @param radial_quantile the quantile cutoff for computing the TPD. Points with
#'    radial quantiles above this value will be used to estimate the TPD.
#'    (default is 0.975)
#' @param max_lag the maximum lag-value to compute the TPD for. Only used with
#'    `vector` inputs. (default is 50)
#' @param radial_thresh instead of specifying the `radial_quantile` users can
#'    specify a threshold such that points with radii larger than the
#'    `radial_thresh` are included in the TPD estimation. `radial_quantile` is
#'    ignored if this is greater than zero. (default is 0)
#' @param get_large change to `FALSE` if you have already removed all small
#'    points. If `true` the function will remove points below the with radii
#'    below the `radial_quantile` so that only large points are used in the
#'    TPD computation. (default is `TRUE`)
#' @param trans_marginal change to `FALSE` if the margins are already
#'    Fr\'echet(2). If `TRUE` the function will transform marginal distributions
#'    to be Fr\'echet(2). (default is `TRUE`)
#' @param gpd_scale fix the scale at your desired scale value during the
#'    marginal transformation. This is only used if `trans_marginal == TRUE`.
#'    If no value is input then it will be estimated. (default is -99)
#' @param gpd_shape fix the shape at your desired shape value. If no value is
#'    input then it will be estimated. Note that gpd shape is 1/alpha so if you
#'    have generated data with alpha=2 you should enter 1/2. This is only used
#'    if `trans_marginal == TRUE`. If no value is input then it will be
#'    estimated. (default is -99)
#' @param fix_alpha used to fix the value of the tail index `alpha` at a user
#'    specified value so that PCA can be done with different `alpha` values.
#'    This overrides the estimation of `alpha` for the `tpd` object.
#' @param marginal_thresh the quantile to use as a cutoff between the ECDF and
#'    GPD components in the marginal transformation. (default is 0.975)
#' @param matrix_as_seasons change to `TRUE` if your time series data are stored
#'    as a matrix with columns representing plausibly iid replicates (e.g.,
#'    seasons). This protects against lagging across seasons. (default is
#'    `FALSE`)
#' @param vector_norm change to `TRUE` to perform the pseudo-polar decomposition
#'    using the full vector-norm for the radial components rather than the
#'    pairwise norm. This functionality is limited at present. (default
#'    is `TRUE`)
#' @param verbose change to `FALSE` if you do not want to see comments.
#'    (default is `TRUE`)
#'
#' @returns a `tpd` S3 object which the estimated tpd as the value and an
#'    estimate of `alpha` (the tail index). When the user supplies a value to
#'    `fix_alpha` then `alpha` is set to that value. Otherwise, `alpha` is set
#'    to 2 when `transform_marginal == TRUE` and is estimated by `evd::fpot()`
#'    for all other cases.
#' @export
#'
#' @references
#' \insertRef{cooley_thibaud_2019}{tpdmethods}
#'
#' @examples
#' myData <- matrix(rnorm(10000), ncol = 2)
#' out1 <- tpd(myData, radial_quantile = 0.95)
#'
#' myMatrix <- matrix(rnorm(10000), ncol = 5)
#' out2 <- tpd(myMatrix, radial_thresh = 2)
#'
#' myVector <- gen_ar1(1000, phi = 0.6)
#' out3 <- tpd(myVector, trans_marginal = FALSE)
tpd <- function(data,
                radial_quantile = 0.975,
                max_lag = 50,
                radial_thresh = 0,
                get_large = TRUE,
                trans_marginal = TRUE,
                gpd_scale = -99,
                gpd_shape = -99,
                fix_alpha = NULL,
                marginal_thresh = 0.975,
                matrix_as_seasons = FALSE,
                vector_norm = FALSE,
                verbose = TRUE){

  if(vector_norm == TRUE){
    if(trans_marginal == TRUE){
      stop("trans_marginal == TRUE and vector_norm == TRUE: we're working on it")
    }
    if(is.matrix(data)){
      if(min(dim(data)) < 2){
        stop("!#!#! Matrix dimension too small, try using a vector !#!#!")
      }
      if(dim(data)[1] == 2){
        data <- t(data)
      }
      if(dim(data)[2] == 2){
        tpds <- tpd.once(data,
                         radial_quantile = radial_quantile,
                         radial_thresh = radial_thresh,
                         get_large = get_large,
                         trans_marginal = trans_marginal,
                         gpd_scale = gpd_scale,
                         gpd_shape = gpd_shape,
                         marginal_thresh = marginal_thresh)
      } else {
        if(verbose == TRUE){
          print("Matrix input, we assume rows represent replicates")
          print(". . . and columns represent variables.")}
        d <- dim(data)[2] # dimension of the problem (e.g., number of stations)
        tpds <- matrix(1, nrow = d, ncol = d)
        radii <- apply(data, 1, function(x){sqrt(sum(x^2))})
        radial_cutoff <- max(radial_thresh,
                             unname(stats::quantile(radii,
                                                    probs = radial_quantile)))
        large_obs <- which(radii > radial_cutoff)
        data      <- data[large_obs, ] / radii[large_obs]
        tictoc::tic(paste0("finished 100 of ", d, " rows of the TPDM"))
        for(j in 1:d){
          for(k in j:d){
            tpds[j, k] <- sum(data[,j] * data[,k]) / length(large_obs)
            tpds[k, j] <- tpds[j, k]
          }
          if(j %% 100 == 0){
            tictoc::toc()
            tictoc::tic(paste0("finished ", j+100, " of ", d, " rows of the TPDM"))
          }
        }
      }
    } else {
      stop("vector_norm = TRUE but not a matrix: we're working on it")
    }
  } else { # Normalize pairwise instead
    if(matrix_as_seasons == TRUE){
      if(is.matrix(data)){
        n_seasons <- dim(data)[2]
        if(verbose == TRUE){
          print(paste0("Assuming this is a time series with ", n_seasons,
                       " seasons."))
          }
        tpds <- numeric(max_lag + 1)
        tpds[1] <- 1
        n_each <- dim(data)[1] # length of each season
        for(i in 1:max_lag){
          lagged_data <- matrix(NA, ncol = 2, nrow = n_seasons * (n_each - i))
          for(j in 1:n_seasons){ # lag within seasons
            lagged_data[((j-1)*(n_each-i) + 1):(j*(n_each-i)), 1] <-
              data[1:(n_each - i), j]
            lagged_data[((j-1)*(n_each-i) + 1):(j*(n_each-i)), 2] <-
              data[(1 + i):n_each, j]
          }
          tpds[i + 1] <- tpd.once(lagged_data, # tpd on bivariate matrix
                                  radial_quantile = radial_quantile,
                                  radial_thresh = radial_thresh,
                                  get_large = get_large,
                                  trans_marginal = trans_marginal,
                                  gpd_scale = gpd_scale,
                                  gpd_shape = gpd_shape,
                                  marginal_thresh = marginal_thresh)
        }
      } else {
        stop("!#!#! matrix_as_seasons is TRUE but object is not a matrix !#!#!")
      }
    } else if(is.matrix(data)){
      if(min(dim(data)) < 2){
        stop("!#!#! Matrix dimension too small, try using a vector !#!#!")
      }
      if(dim(data)[1] == 2){
        data <- t(data)
      }
      if(dim(data)[2] == 2){
        tpds <- tpd.once(data,
                         radial_quantile = radial_quantile,
                         radial_thresh = radial_thresh,
                         get_large = get_large,
                         trans_marginal = trans_marginal,
                         gpd_scale = gpd_scale,
                         gpd_shape = gpd_shape,
                         marginal_thresh = marginal_thresh)
      } else {
        if(verbose == TRUE){
          print("Matrix input, we assume rows represent replicates")
          print(". . . and columns represent variables.")}
        d <- dim(data)[2] # dimension of the problem (e.g., number of stations)
        tpds <- matrix(1, nrow = d, ncol = d)
        tictoc::tic(paste0("finished 100 of ", d, " rows of the TPDM"))
        for(j in 1:(d - 1)){
          for(k in (j+1):d){
            temp_data <- data[, c(j,k)]
            tpds[j, k] <- tpd.once(temp_data,
                                   radial_quantile = radial_quantile,
                                   radial_thresh = radial_thresh,
                                   get_large = get_large,
                                   trans_marginal = trans_marginal,
                                   gpd_scale = gpd_scale,
                                   gpd_shape = gpd_shape,
                                   marginal_thresh = marginal_thresh)
            tpds[k, j] <- tpds[j, k]
          }
          if(j %% 100 == 0){
            tictoc::toc()
            tictoc::tic(paste0("finished ", j+100, " of ", d, " rows of the TPDM"))
          }
        }
      }
    } else if(is.vector(data)){
      if(verbose == TRUE){print("Vector input, we assume this is a time series")}
      tpds <- numeric(max_lag + 1)
      tpds[1] <- 1
      n <- length(data)
      for(i in 1:max_lag){
        temp_data <- matrix(NA, nrow = n - i, ncol = 2)
        temp_data[,1] <- data[1:(n - i)]
        temp_data[,2] <- data[(1 + i):n]
        tpds[i + 1] <- tpd.once(temp_data,
                                radial_quantile = radial_quantile,
                                radial_thresh = radial_thresh,
                                get_large = get_large,
                                trans_marginal = trans_marginal,
                                gpd_scale = gpd_scale,
                                gpd_shape = gpd_shape,
                                marginal_thresh = marginal_thresh)
      }
    } else if(is.list(data)){
     all_same_class <- all(sapply(data,
                                  function(x){class(x) == class(data[[1]])}))
     all_same_dimension <- all(sapply(data,
                                      function(x){dim(x) == dim(data[[1]])}))
     if(all_same_class & all_same_dimension){
       if(dim(data[[1]])[1] == 2){
         data <- lapply(data, function(x){t(x)}) # transpose list objects
       }
       if(dim(data[[1]])[2] == 2){
         tpds <- sapply(data,
                        FUN = tpd.once,
                        radial_quantile = radial_quantile,
                        radial_thresh = radial_thresh,
                        get_large = get_large,
                        trans_marginal = trans_marginal,
                        gpd_scale = gpd_scale,
                        gpd_shape = gpd_shape,
                        marginal_thresh = marginal_thresh)
       } else {
         stop("!#!#! Objects in list are not bivariate !#!#!")
       }
     } else {
       stop("!#!#! Objects in list are different types or sizes !#!#!")
     }
    }
  }

  # get alpha for the `tpd` object - needs work
  temp_alpha <- 2
  if(!is.null(fix_alpha)){
    temp_alpha <- fix_alpha
  } else if(trans_marginal != TRUE){
    temp_fit   <- evd::fpot(c(data), threshold = quantile(data, probs = 0.975))
    temp_alpha <- 1/unname(temp_fit$estimate[2])
  }

  tpds <- validate_tpd(new_tpd(x = tpds, alpha = temp_alpha))

  return(tpds)
}



#' Get the Tail Pairwise Dependence of one pair of variables
#'
#' This function takes a bivariate matrix as input and estimates the TPD. If
#'    necessary it transforms the marginal and gets the large points first.
#'
#' @param data a bivariate `matrix`
#' @param radial_quantile the quantile cutoff for computing the TPD. Points with
#'    radial quantiles above this value will be used to estimate the TPD.
#'    (default is 0.975)
#' @param radial_thresh instead of specifying the `radial_quantile` users can
#'    specify a threshold such that points with radii larger than the
#'    `radial_thresh` are included in the TPD estimation. `radial_quantile` is
#'    ignored if this is greater than zero. (default is 0)
#' @param get_large change to `FALSE` if you have already removed all small
#'    points. If `true` the function will remove points below the with radii
#'    below the `radial_quantile` so that only large points are used in the
#'    TPD computation. (default is `TRUE`)
#' @param trans_marginal change to `FALSE` if the margins are already
#'    Fr\'echet(2). If `TRUE` the function will transform marginal distributions
#'    to be Fr\'echet(2). (default is `TRUE`)
#' @param gpd_scale fix the scale at your desired scale value during the
#'    marginal transformation. This is only used if `trans_marginal == TRUE`.
#'    If no value is input then it will be estimated. (default is -99)
#' @param gpd_shape fix the shape at your desired shape value. If no value is
#'    input then it will be estimated. Note that gpd shape is 1/alpha so if you
#'    have generated data with alpha=2 you should enter 1/2. This is only used
#'    if `trans_marginal == TRUE`. If no value is input then it will be
#'    estimated. (default is -99)
#' @param marginal_thresh the quantile to use as a cutoff between the ECDF and
#'    GPD components in the marginal transformation. (default is 0.975)
#'
#' @returns the estimated TPD value
#' @export
#'
#' @examples
#' myData <- matrix(rnorm(10000), ncol = 2)
#' out1 <- tpd(myData, radial_quantile = 0.95)
tpd.once <- function(data,
                     radial_quantile = 0.975,
                     radial_thresh = 0,
                     get_large = TRUE,
                     trans_marginal = TRUE,
                     gpd_scale = -99,
                     gpd_shape = -99,
                     marginal_thresh = 0.975){

  if(trans_marginal == TRUE){
    data <- apply(data, 2,
                  transform_marginal,
                  q = marginal_thresh,
                  gpd_scale = gpd_scale,
                  gpd_shape = gpd_shape)
  }
  if(get_large == TRUE){
    data <- tpd.get_large(data, thresh = radial_quantile, r0 = radial_thresh)
  }
  return(tpd.est(data))
}


# This function estimates the TPD from a bivariate matrix of large RV(2) points.
#' Estimate the Tail Pairwise Dependence from large points
#'
#' This function estimates the TPD from a bivariate `matrix` of large regularly
#'    varying (alpha = 2) points.
#'
#' @param data a bivariate `matrix` of large RV(2) points.
#'
#' @returns the estimated TPD value
#' @export
#'
#' @references
#' \insertRef{cooley_thibaud_2019}{tpdmethods}
#'
#' @examples
#' myData <- matrix(evd::rfrechet(1000, shape = 2), ncol = 2)
#' radii <- sqrt(rowSums(myData^2))
#' largeData <- myData[order(radii)[450:500], ]
#' out <- tpd.est(largeData)
tpd.est <- function(data){
  temp_r <- rowSums(data^2)
  return(2*mean(data[,1]*data[,2]/temp_r))
}


#
#' Get large points from a bivariate `matrix`
#'
#' This function extracts the large points from a bivariate `matrix`
#'
#' @param data a bivariate `matrix` of points
#' @param thresh the radial quantile used to distinguish large points from
#'    points in the bulk of the distribution. (default is 0.975)
#' @param r0 instead of setting a quantile the user may set a radial value. If
#'    `r0` is greater than zero the `thresh` is ignored. (default is 0).
#'
#' @returns a bivariate `matrix` of large points.
#' @export
#'
#' @examples
#' myData <- matrix(evd::rfrechet(1000, shape = 2), ncol = 2)
#' out <- tpd.get_large(myData)
tpd.get_large <- function(data, thresh = 0.975, r0 = 0){
  radii <- sqrt(rowSums(data^2))
  if(r0 > 0){
    r00 <- r0 # user defined threshold
  } else {
    r00 <- sort(radii)[floor(thresh*(length(radii)+1))]
  }
  return(data[which(radii > r00), ]) # return large points
}



#' Constructor for the S3 class `tpd`.
#'
#' @param x the `numeric` values of the tpd object. A `vector` for a tpdf and a
#'    `matrix` for a tpdm.
#' @param alpha the tail index of the underlying data.
#'
#' @returns a `tpd` S3 object
#' @export
#'
#' @examples
#' my_tpd <- new_tpd(x = matrix(c(0.7, 1, 1, 0.7), nrow = 2), alpha = 2)
new_tpd <- function(x = double(), alpha = numeric(1)){
  stopifnot(is.double(x))
  stopifnot(is.numeric(alpha))
  stopifnot(NROW(alpha) == 1 && NCOL(alpha) == 1)

  structure(x,
            class = c("tpd", class(x)),
            alpha = alpha
  )
}

#' Validator for the S3 class `tpd`
#'
#' @param x the `tpd` object to be validated
#'
#' @returns either an ERROR or the `tpd` object
#' @export
#'
#' @examples
#' my_tpd <- new_tpd(x = matrix(c(0.7, 1, 1, 0.7), nrow = 2), alpha = 2)
#' my_validated_tpd <- validate_tpd(my_tpd)
validate_tpd <- function(x){
  values     <- unclass(x)
  tail_index <- attr(x, "alpha")

  if(!all(!is.na(values)) || any(values < 0)){
    stop(
      "All `tpd` values must be non-missing and non-negative",
      call. = FALSE
    )
  }

  if(tail_index <= 0){
    stop(
      "The `tpd` requires heavy-tailed data but the `alpha` is non-positive",
      call. = FALSE
    )
  }

  x
}

#' Check whether `x` is a `tpd` S3 object
#'
#' @param x an object
#'
#' @returns `TRUE` if `x` is a `numeric` `tpd` object, `FALSE` otherwise
#' @export
#'
#' @examples
#' my_tpd <- new_tpd(x = matrix(c(0.7, 1, 1, 0.7), nrow = 2), alpha = 2)
#' is.tpd(my_tpd)
#'
#' is.tpd(1:5)
is.tpd <- function(x){
  is.numeric(x) && inherits(x, "tpd")
}


#' Create a plot of confidence intervals for alpha from each variable (station)
#'
#' The alpha (index of regular variation) for each variable is estimated with
#'    the Hill estimator with a user specified `k`. In addition, the `ci`%
#'    confidence intervals are computed. These are plotted with the a horizontal
#'    line which is the estimate for a single joint alpha. The joint alpha is
#'    computed as the average of the individual estimates because `k` is the
#'    same for each variable. If a user wants `k` to differ then this tool is
#'    inappropriate.
#'
#' @param data the `matrix` of observations with variables (e.g., stations) in
#'    the columns.
#' @param k the number of largest order statistics to be used with the Hill
#'    estimator
#' @param ci the confidence level (as a proportion) for the confidence intervals
#'    that will be plotted (the default is 0.95 which makes 95% CIs).
#'
#' @returns a ggplot2 `plot` object which shows 100*`ci`% confidence intervals
#'    for the tail index (alpha) of each variable computed using the Hill
#'    estimator as well as an estimate of the shared (joint) alpha under the
#'    assumption that the alpha is the same for each variable.
#' @export
#'
#' @references
#' \insertRef{hill1975}{tpdmethods}
#'
#' @examples
#' alpha_plot(financial_data[, -1], k = 300)
alpha_plot <- function(data, k, ci = 0.95){
  if(!is.null(names(data))){
    index <- names(data)
  } else {
    index  <- 1:NCOL(data)
  }
  df            <- data.frame(index = index)
  df$alpha_hats <- apply(data, 2, hill_est, k = k)
  df$ses        <- df$alpha_hats / sqrt(k)
  crit_val      <- qnorm(1 - (1 - ci)/2)
  df$upper      <- df$alpha_hats + crit_val * df$ses
  df$lower      <- df$alpha_hats - crit_val * df$ses

  temp_plot <- ggplot2::ggplot(df,
                               ggplot2::aes(x = index,
                                            y = alpha_hats,
                                            ymin = lower,
                                            ymax = upper)) +
    ggplot2::geom_pointrange() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = mean(alpha_hats)), color = 2) +
    ggplot2::annotate(geom = "text",
                      x = 1.1,
                      y = mean(df$alpha_hats) + 0.05,
                      label = paste(round(mean(df$alpha_hats), 3)),
                      color = 2) +
    ggplot2::labs(title = paste0(100*ci, "% confidence intervals for alpha"),
                  subtitle = paste0("(Hill estimator with k=", k, ")"),
                  x = "",
                  y = expression(hat(alpha))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30))
  return(temp_plot)
}



#' The Hill estimator for alpha
#'
#' Estimates the index of regular variation (alpha) given a vector of
#'    observations `x` and a number of large order statistics `k` to be used.
#'
#' @param x the `vector` of observations to be used to estimate alpha
#' @param k the number of large order statistics to be used
#'
#' @returns a scalar value which is the point estimate for alpha
#' @export
#'
#' @references
#' \insertRef{hill1975}{tpdmethods}
#'
#' @examples
#' hill_est(financial_data[,2], k = 300)
hill_est <- function(x, k){
  ordered <- rev(sort(x))
  ordered <- ordered[ordered > 0][1:k]
  return(1 / (mean(log(ordered)) - log(ordered[k])))
}


#' Estimates the joint alpha for an entire data matrix
#'
#' This function computes the hill estimate for each column in the `data` and
#'    then takes the `mean()` because `k` is shared. This is the estimate for
#'    alpha under the assumption that all of the variable have the same alpha.
#'
#' @param data the `matrix` of observations to be used to estimate alpha
#' @param k the number of large order statistics to be used for each variable in
#'    `data` (this value must be a scalar).
#'
#' @returns a scalar estimate of the joint alpha rounded to four decimal places.
#' @export
#'
#' @examples
#' alpha_hat(financial_data[,-1], k = 300)
alpha_hat <- function(data, k){
  alpha_hats <- apply(data, 2, hill_est, k = k)
  round(mean(alpha_hats), 4)
}


# # The following function is what was used in the generate_tlets.R script for
# #   the proxy-likelihood work. It is here as a reference.
# TPDF <- function(ts, thresh = 0.975, maxlag = 50){
#   tpdf <- rep(0, maxlag+1)
#   if(is.list(ts)){
#     tpdf[1] <- 1
#     for(i in 2:(length(ts)+1)){
#       temp_mat <- ts[[i-1]]
#       temp_r   <- rowSums(temp_mat^2)
#       tpdf[i]  <- 2*mean(temp_mat[, 1]*temp_mat[, 2]/temp_r)
#     }
#   } else {
#     rt0 <- sqrt(ts^2 + ts^2)
#     r00 <- sort(rt0)[floor(thresh*(length(rt0)+1))]
#     tpdf <- rep(0, maxlag+1)
#     tpdf[1] <- 2/length(which(rt0>r00)) * sum((ts^2)[which(rt0>r00)] /
#                                                 (rt0[which(rt0>r00)]^2))
#     n <- length(ts)
#     rt <- matrix(NA, nrow = n-1, ncol = maxlag)
#     r0 <- rep(NA, maxlag)
#     for(i in 1:maxlag){
#       temp_means <- c(mean(ts[1:(n-i)]), mean(ts[(1+i):n]))
#       temp_mat <- matrix(
#         c(pmax(rep(0, (n-i)), ts[1:(n-i)] - temp_means[1]),
#           pmax(rep(0, (n-i)), ts[(1+i):n] - temp_means[2])),
#         ncol = 2, byrow = F)
#       rt[,i] <- c(sqrt(temp_mat[,1]^2 + temp_mat[,2]^2), rep(NA, i-1))
#       temp <- na.omit(rt[,i])
#       if(length(temp)){ # ignore zero-length columns
#         k <- floor(thresh*(length(temp)+1))
#         r0[i] <- sort(temp)[k]
#       } else {r0[i] <- 0}
#       indices <- which(rt[,i]>r0[i])
#       tpdf[i+1] <- 2/length(indices) * sum((temp_mat[,1]*temp_mat[,2])[indices] /
#                                              (rt[indices, i]^2))
#     }
#   }
#   return(tpdf)
# }
