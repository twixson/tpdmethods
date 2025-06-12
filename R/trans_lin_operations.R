
# This function is the softplus function used to map the reals to the positive
#   half-line.
f <- function(y){
  if(min(y) < -36){
    print("*!*!* mass at machine zero *!*!*")
  }
  y <- ifelse(y < -36, -36, y)
  return(ifelse(y > 100, y, log(1+exp(y))))
}

# This function is the inverse of the softplus. It maps the positive half-line
#   to the reals.
finv <- function(x){
  if(min(x) < 0){
    stop("*!*!* negative value *!*!*")
  } else if(min(x) < 2e-16){
    print("*!*!* mass at machine zero *!*!*")
  }
  x <- ifelse(x < 2e-16, 2e-16, x)
  return(ifelse(x > 100, x, log(exp(x)-1)))
}

# This function performes transformed-linear addition.
tadd <- function(a,b){
  return(f(finv(a)+finv(b)))
}

# This function performs transformed-linear scalar multiplication between for
#   scalar "a" and vector "b"
tmult <- function(a,b){
  return(f(a*finv(b)))
}
