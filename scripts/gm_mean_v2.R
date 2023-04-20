gm_mean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  }
  if(length(x) == 0){
    return(0)
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x[x>0 & !is.na(x)]))
  }
}
