csv2dt <- function(exampleSpec_csv){
  tempdata <- read.csv(exampleSpec_csv, header = FALSE)
  all_mzs  <- tempdata[[1]]
  all_ints <- as.numeric(tempdata[[2]])
  all_ints = 100*all_ints/max(all_ints)
  min_signal = 1
  counter = 0
  mz = NULL
  int = NULL
  for(j in 1:length(all_ints)){
    if(as.numeric(all_ints[j])>=min_signal){
      counter = counter + 1
      mz[counter] = as.numeric(all_mzs[j])
      int[counter] = as.numeric(all_ints[j])
    }
  }
  return(data.frame(mz,int))
}