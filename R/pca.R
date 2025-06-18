
#' Get the eigenvalues and eigenvectors from a `matrix` of TPD values
#'
#' This is a wrapper of the `eigen()` function that includes a step to find
#'    the closest PD matrix to the input TPDM using the `nearPD()` function.
#'
#' @param tpdm a square `matrix` of TPD values
#' @param make_pd change to `FALSE` if you do not want to find a PD matrix
#'    that is close to the input matrix.
#'
#' @returns a `list` with two components. The `values` are the eigenvalues of
#'    the TPDM. The `vectors` are the eigenvectors of the TPDM.
#' @export
#'
#' @examples
#' set.seed(22)
#' myData <- matrix(rnorm(625), ncol = 5)
#' myTPDM <- tpd(myData)
#' out <- get_eigen(myTPDM)
get_eigen <- function(tpdm, make_pd = TRUE){
  if(make_pd == TRUE){
    tpdm <- make_pd(tpdm)$mat
  }
  temp_out <- eigen(tpdm)

  for(i in 1:length(temp_out$values)){
    if(max(temp_out$vectors[,i]) < 0){
      temp_out$vectors[,i] <- -temp_out$vectors[,i]
    }
  }
  return(temp_out)
}


#' Find a positive definite `matrix` that is close to the input `matrix`.
#'
#' This is just a wrapper for the `Matrix` function `nearPD()`
#'
#' @param x square `matrix` that may or may not be PD
#'
#' @returns a `matrix` that is PD
#' @export
#'
#' @examples
#' myData <- matrix(rnorm(25), ncol = 5)
#' out <- make_pd(myData)
make_pd <- function(x){
  Matrix::nearPD(x)
}


#' Plot a single eigenvector.
#'
#' This is a wrapper fuction to plot a map showing the value of a single
#'    eigenvector using latitude and longitude.
#' Note that eigenvalues and eigenvectors are separated for easier
#'    plotting of transformed-linear reconstructions.
#'
#' @param eigen_num which eigenvector you want to plot
#' @param eigen_vecs the `matrix` of eigenvectors
#' @param eigen_vals the `vector` of eigenvalues
#' @param lat a `vector` of latitudes
#' @param lon a `vector` of longitudes
#' @param same_scale change to `TRUE` if you want to force several plots to
#'    have the same scale. (default is `FALSE`)
#' @param scale_lims a bivariate `vector` of the form (scale minimum, scale
#'    maximum). (default is c(-99, -99))
#'
#' @returns a ggplot2 `plot` object
#' @export
#'
#' @examples
#' myData <- matrix(rnorm(625), ncol = 5)
#' myTPDM <- tpd(myData)
#' myEigen <- get_eigen(myTPDM)
#' plot1 <- plot_eigen(1,
#'                     eigen_vecs = myEigen$vectors,
#'                     eigen_vals = myEigen$values,
#'                     lat = 1:25,
#'                     lon = 1:25)
#' plot1
plot_eigen <- function(eigen_num,
                       eigen_vecs,
                       eigen_vals,
                       lat,
                       lon,
                       same_scale = FALSE,
                       scale_lims = c(-99, -99)){
  temp_df <- data.frame("lat" = lat,
                        "lon" = lon,
                        "eigen" = eigen_vecs[,eigen_num])
  my_sf <- sf::st_as_sf(temp_df, coords = c("lon", "lat"), crs = 4326)

  if(same_scale == TRUE){
    if(max(scale_lims) > -99){
      plot1 <- ggplot2::ggplot(my_sf) +
        ggplot2::geom_sf(ggplot2::aes(col = eigen), shape = 15, size = 2) +
        ggplot2::scale_color_viridis_c(limits = scale_lims) +
        ggplot2::labs(title = paste0("Map of Eigenvector ", eigen_num),
             subtitle = paste0("Eigenvalue = ", round(eigen_vals[eigen_num], 2))) +
        ggplot2::theme(axis.text = ggplot2::element_text(angle = 45))
    }
    plot1 <- ggplot2::ggplot(my_sf) +
      ggplot2::geom_sf(ggplot2::aes(col = eigen), shape = 15, size = 2) +
      ggplot2::scale_color_viridis_c("limits = c(-0.037, 0.037)") +
      ggplot2::labs(title = paste0("Map of Eigenvector ", eigen_num),
           subtitle = paste0("Eigenvalue = ", round(eigen_vals[eigen_num], 2))) +
      ggplot2::theme(axis.text = ggplot2::element_text(angle = 45))
  } else {
    plot1 <- ggplot2::ggplot(my_sf) +
      ggplot2::geom_sf(ggplot2::aes(col = eigen), shape = 15, size = 2) +
      ggplot2::scale_color_viridis_c() +
      ggplot2::labs(title = paste0("Map of Eigenvector ", eigen_num),
           subtitle = paste0("Eigenvalue = ", round(eigen_vals[eigen_num], 2))) +
      ggplot2::theme(axis.text = ggplot2::element_text(angle = 45))
  }


  return(plot1)
}

#' Get the principal components of a data vector
#'
#' @param data_vec the data vector
#' @param eigenvecs the eigenvectors of the TPDM of the data vector
#'
#' @returns a `vector` of principal components
#' @export
#'
#' @references
#' \insertRef{cooley_thibaud_2019}{tpdmethods}
#'
#' @examples
#' myData <- matrix(rnorm(625), ncol = 5)
#' myTPDM <- tpd(myData)
#' myEigen <- get_eigen(myTPDM)
#' myDatavec <- abs(rnorm(5))
#' out <- get_pcs(myDatavec, eigenvecs = myEigen$vectors)
get_pcs <- function(data_vec, eigenvecs){
  data_vec[which(data_vec == 0)] <- 0.00001
  t(eigenvecs) %*% finv(data_vec)
}



#' Reconstruct a data vector with a subset of the eigenvectors.
#'
#' This function allows for the partial reconstruction of a data vector (for
#'    example, of a big storm) using the `num_eigen` largest eigenvectors.
#'
#' @param pcs the principal components of the data vector to be reconstructed.
#'    Often this is the output of `get_pcs()`.
#' @param eigenvecs the eigenvectors of the TPDM of the data vector to be
#'    reconstructed.
#' @param num_eigen the number of eigenvectors to be used in the partial
#'    reconstruction
#'
#' @returns a `vector` that contains the values of the partially reconstructed
#'    data vector.
#' @export
#'
#' @references
#' \insertRef{cooley_thibaud_2019}{tpdmethods}
#'
#' @examples
#' myData <- matrix(rnorm(625), ncol = 5)
#' myTPDM <- tpd(myData)
#' myEigen <- get_eigen(myTPDM)
#' myDatavec <- abs(rnorm(5))
#' myPCs <- get_pcs(myDatavec, eigenvecs = myEigen$vectors)
#' out <- reconstruct_data(myPCs, eigenvecs = myEigen$vectors, num_eigen = 3)
reconstruct_data <- function(pcs, eigenvecs, num_eigen){
  temp_val <- pcs[1] * eigenvecs[,1]
  if(num_eigen > 1){
    for(i in 2:num_eigen){
      temp_val <- temp_val + (pcs[i] * eigenvecs[,i])
    }
  }
  return(f(temp_val))
}



#' Plot a `vector` of data using latitude and longitude
#'
#' This is a wrapper function for creating and plotting an `sf` object. We use
#'    it to plot partial reconstructions created with `reconstruct_data()`.
#'
#' @param data the `vector` of data to be plotted
#' @param plot_title the name of the plot
#' @param lat a `vector` of latitudes
#' @param lon a `vector` of longitudes
#'
#' @returns a ggplot2 `plot` object
#' @export
#'
#' @examples
#' myData <- rnorm(25)
#' out <- plot_sf(myData, plot_title = "Example Plot", lat = 1:5, lon = 1:5)
plot_sf <- function(data, plot_title, lat, lon){
  temp_df <- data.frame("lat" = lat,
                        "lon" = lon,
                        "value" = data)
  my_sf <- sf::st_as_sf(temp_df, coords = c("lon", "lat"), crs = 4326)

  plot1 <- ggplot2::ggplot(my_sf) +
    ggplot2::geom_sf(ggplot2::aes(col = value), shape = 15, size = 2) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(title = plot_title) +
    ggplot2::theme(axis.text = ggplot2::element_text(angle = 45))
}


#' Compute the time series of a single principle component.
#'
#' @param eigen_num the principal component to compute
#' @param eigen_vecs the `matrix` of eigenvectors.
#' @param data the full `matrix` of data.
#'
#' @returns a `vector` of the principal components for each observation.
#' @export
#'
#' @references
#' \insertRef{cooley_thibaud_2019}{tpdmethods}
#'
#' @examples
#' myData <- matrix(rnorm(625), ncol = 5)
#' myTPDM <- tpd(myData)
#' myEigen <- get_eigen(myTPDM)
#' out <- get_pc_ts(1, eigen_vecs = myEigen$vectors, data = myData)
get_pc_ts <- function(eigen_num, eigen_vecs, data){
  return(data %*% eigen_vecs[,eigen_num])
}

#' Plot a time series of a single principle component
#'
#' This function is a wrapper to compute and plot a time series for a single
#'    principle component.
#'
#' @param eigen_num the principal component to plot.
#' @param eigen_vecs the `matrix` of eigenvectors.
#' @param data the full `matrix` of data.
#' @param dates the dates to put on the x-axis. If `dates == FALSE` no a default
#'    x-axis will be used.
#'
#' @returns a ggplot2 `plot` object
#' @export
#'
#' @references
#' \insertRef{cooley_thibaud_2019}{tpdmethods}
#'
#' @examples
#' myData <- matrix(rnorm(625), ncol = 5)
#' myTPDM <- tpd(myData)
#' myEigen <- get_eigen(myTPDM)
#' out <- get_ts_plot(eigen_num = 1,
#'                    eigen_vecs = myEigen$vectors,
#'                    data = myData,
#'                    dates = FALSE)
get_ts_plot <- function(eigen_num,
                        eigen_vecs,
                        data,
                        dates){
  ts_vec <- data %*% eigen_vecs[,eigen_num]
  temp_df <- data.frame("Date" = dates,
                        "PC_Score" = ts_vec,
                        "x_index" = 1:length(dates))
  x_ticks <- grep("0331", dates)

  if(dates[1] == FALSE){
    temp_plot <- ggplot2::ggplot(temp_df,
                                 ggplot2::aes(x = x_index, y = PC_Score)) +
      ggplot2::geom_line() +
      ggplot2::ylab("PC Score") +
      ggplot2::xlab("Date") +
      ggplot2::labs(title = paste0("Principle Component ", eigen_num)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30))
  } else {
    temp_plot <- ggplot2::ggplot(temp_df,
                                 ggplot2::aes(x = x_index, y = PC_Score)) +
      ggplot2::geom_line() +
      ggplot2::ylab("PC Score") +
      ggplot2::xlab("Date") +
      ggplot2::scale_x_continuous(breaks = x_ticks[seq(5, 40, length.out = 8)],
                         label = substr(temp_df$Date[x_ticks], start = 1,
                                        stop = 4)[seq(5, 40, length.out = 8)]) +
      ggplot2::labs(title = paste0("Principle Component ", eigen_num)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30))
  }

  return(temp_plot)
}
