asm_DT <- function(qSpec,emma_tree){

  nsteps = dim(emma_tree)[1]
  i = 1
  while(i <= nsteps){
    a = which(round(qSpec$mz)==emma_tree$mz[i])
    b = which(qSpec$int>emma_tree$int[i])
    d = intersect(a,b)
    if(length(d)>0){ # this means the answer at the decision point i is yes
      test1 = is.na(suppressWarnings(as.numeric(emma_tree$yes[i])))
      if(test1==FALSE){
        i = as.numeric(emma_tree$yes[i]) # this means we go to a new i
      } else {
        return(data.frame(mz = emma_tree$mz[i],
                          ion = emma_tree$yes[i],
                          R = emma_tree$R_yes[i]))
      }
    } else { # this means the answer at decision point i is no
      test2 = is.na(suppressWarnings(as.numeric(emma_tree$no[i])))
      if(test2==FALSE){
        i = as.numeric(emma_tree$no[i]) # this means we go to a new i
      } else {
        return(data.frame(mz = emma_tree$mz[i],
                          ion = emma_tree$no[i],
                          R = emma_tree$R_no[i]))
      }
    }
  }

  
  
}