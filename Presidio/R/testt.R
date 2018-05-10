f<-function(lower,uper,log = TRUE){
  if (length(lower) == 2 && missing(upper)) {
    upper <- lower[2]
    lower <- lower[1]
  }
  n <- length(lower)
  if (length(upper) != n)
    stop("'lower' and 'upper' both be the same length")
  p.in <- 1/(upper - lower)
  p.out <- 0
  if (log) {
    p.in <- log(p.in)
    p.out <- -Inf
  }
  res<-function(x) {
    ret <- rep(p.in, length.out = length(x))
    ret[x < lower | x > upper] <- p.out
    if (log)
      sum(ret)
    else prod(ret)
  }
  return(res)
}
f(lower,upper)

t<-function (lower, upper, log = TRUE)
{
  if (length(lower) == 2 && missing(upper)) {
    upper <- lower[2]
    lower <- lower[1]
  }
  n <- length(lower)
  if (length(upper) != n)
    stop("'lower' and 'upper' both be the same length")
  p.in <- 1/(upper - lower)
  p.out <- 0
  if (log) {
    p.in <- log(p.in)
    p.out <- -Inf
  }
  function(x) {
    ret <- rep(p.in, length.out = length(x))
    ret[x < lower | x > upper] <- p.out
    if (log)
      sum(ret)
    else prod(ret)
  }
}
t(lower,upper)
