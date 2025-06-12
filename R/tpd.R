
# this function puts the users input into a form recognizable by simpler
#   functions which transform the marginal distribution, extract the large
#   points, and compute the TPD value. It currently recognizes matrices,
#   vectors as time series, and lists as
tpd <- function(data,
                radial_quantile = 0.975,
                max_lag = 50,
                radial_thresh = 0,
                get_large = TRUE,
                trans_marginal = TRUE,
                marginal_thresh = 0.975,
                verbose = TRUE){
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
                       marginal_thresh = marginal_thresh)
    } else {
      if(verbose == TRUE){
        print("Matrix input, we assume rows represent replicates.")}
      d <- dim(data)[2] # dimension of the problem (e.g., number of stations)
      tpds <- matrix(1, nrow = d, ncol = d)
      for(j in 1:(d - 1)){
        for(k in (j+1):d){
          temp_data <- data[, c(j,k)]
          tpds[j, k] <- tpd.once(temp_data,
                                 radial_quantile = radial_quantile,
                                 radial_thresh = radial_thresh,
                                 get_large = get_large,
                                 trans_marginal = trans_marginal,
                                 marginal_thresh = marginal_thresh)
          tpds[k, j] <- tpds[j, k]
        }
      }
    }
  }

  if(is.vector(data)){
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
                              marginal_thresh = marginal_thresh)
    }
  }

  if(is.list(data)){
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
                      marginal_thresh = marginal_thresh)
     } else {
       stop("!#!#! Objects in list are not bivariate !#!#!")
     }
   } else {
     stop("!#!#! Objects in list are different types or sizes !#!#!")
   }
 }
}

# This function takes a bivariate matrix as input and estimates the TPD. If
#   necessary it transforms the marginal and gets the large points first.
tpd.once <- function(data,
                     radial_quantile = 0.975,
                     radial_thresh = 0,
                     get_large = TRUE,
                     trans_marginal = TRUE,
                     marginal_thresh = 0.975){

  if(trans_marginal == TRUE){
    data <- transform_marginal(data, q = marginal_thresh)
  }
  if(get_large == TRUE){
    data <- tpd.get_large(data, thresh = radial_quantile, r0 = radial_thresh)
  }
  return(tpd.est(data))
}


# This function estimates the TPD from a bivariate matrix of large RV(2) points.
tdp.est <- function(data){
  temp_r <- rowSums(data^2)
  return(2*mean(data[,1]*data[.2]/temp_r))
}


# This function extracts the large points from a bivariate matrix
tpd.get_large <- function(data, thresh = 0.975, r0 = 0){
  radii <- sqrt(rowSums(data^2))
  if(r0 > 0){
    r00 <- r0 # user defined threshold
  } else {
    r00 <- sort(radii)[floor(thresh*(length(radii)+1))]
  }
  return(data[which(rt0 > r00), ]) # return large points
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
