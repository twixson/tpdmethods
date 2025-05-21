f <- function(y){
  if(min(y) < -36){
    print("*!*!* mass at machine zero *!*!*")
  }
  y <- ifelse(y < -36, -36, y)
  return(ifelse(y > 100, y, log(1+exp(y))))
}

finv <- function(x){
  if(min(x) < 0){
    stop("*!*!* negative value *!*!*")
  } else if(min(x) < 2e-16){
    print("*!*!* mass at machine zero *!*!*")
  }
  x <- ifelse(x < 2e-16, 2e-16, x)
  return(ifelse(x > 100, x, log(exp(x)-1)))
}

tadd <- function(a,b){
  return(f(finv(a)+finv(b)))
}

tmult <- function(a,b){
  return(f(a*finv(b)))
}
