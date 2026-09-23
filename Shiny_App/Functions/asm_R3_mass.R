asm_Subtraction <- function(PM,mI1,mI2){
  if(PM == "unknown"){
    return("unknown")
  }
  
  if(mI1=="unknown"){
    return("unknown")
  }
  
  if(mI2=="unknown"){
    return("unknown")
  }
  
  R3 = PM-mI1-mI2
  
  return(R3)
}