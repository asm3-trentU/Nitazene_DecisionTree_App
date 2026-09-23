dt2PM <- function(MS1){
  return(round(MS1$mz[which.max(MS1$int)]))
}