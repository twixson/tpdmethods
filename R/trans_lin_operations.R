
#
#
#' Transformed-linear transformation function
#'
#' This function is the softplus function used to map the reals to the positive
#'    half-line.
#'
#' @param y real number to be transformed
#'
#' @returns a scalar on the positive half-line
#' @export
#'
#' @examples
#' f(2.5)
#' f(-1.2)
f <- function(y){
  if(min(y) < -36){
    print("*!*!* mass at machine zero *!*!*")
  }
  y <- ifelse(y < -36, -36, y)
  return(ifelse(y > 100, y, log(1+exp(y))))
}

#
#
#' Transformed-linear inverse transformation
#'
#' This function is the inverse of the softplus (`f`). It maps the positive
#'    half-line to the reals.
#'
#' @param x a non-negative number
#'
#' @returns a real number
#' @export
#'
#' @examples
#' finv(2.5)
#' # finv(-0.1) # this should return an error.
finv <- function(x){
  if(min(x) < 0){
    stop("*!*!* negative value *!*!*")
  } else if(min(x) < 2e-16){
    print("*!*!* mass at machine zero *!*!*")
  }
  x <- ifelse(x < 2e-16, 2e-16, x)
  return(ifelse(x > 100, x, log(exp(x)-1)))
}


#' Transformed-linear addition
#'
#' This function performes transformed-linear addition.
#'
#' @param a a non-negative number
#' @param b a non-negative number
#'
#' @returns a non-negative number
#' @export
#'
#' @references
#' \insertRef{cooley_thibaud_2019}{tpdmethods}
#'
#' @examples
#' tadd(1.2, 0.05)
#' tadd(125, 23.1)
tadd <- function(a,b){
  return(f(finv(a)+finv(b)))
}


#' Transformed-linear scalar multiplication
#'
#' This function performs transformed-linear scalar multiplication between for
#'    scalar `a` and vector `b`
#'
#' @param a scalar
#' @param b `vector` of non-negative numbers
#'
#' @returns a `vector`
#' @export
#'
#' @references
#' \insertRef{cooley_thibaud_2019}{tpdmethods}
#'
#' @examples
#' tmult(0.5, c(1.2, 2.5, 0.3, 25.1))
tmult <- function(a,b){
  return(f(a*finv(b)))
}
