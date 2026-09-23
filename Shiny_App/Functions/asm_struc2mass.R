asm_struc2mass <- function(struc){
  
  if(struc=="unknown"){
    return(NA)  
  }
  
  atomic_weights <- c(
    H  = 1,
    C  = 12,
    N  = 14,
    O  = 16,
    F  = 19,
    Cl = 35
  )
  
  # Match element symbols and counts
  matches <- gregexpr("(Cl|C|H|O|N|F)(\\d*)", struc, perl = TRUE)
  tokens <- regmatches(struc, matches)[[1]]
  
  mw <- 0
  
  for (token in tokens) {
    
    parts <- regmatches(
      token,
      regexec("(Cl|C|H|O|N|F)(\\d*)", token, perl = TRUE)
    )[[1]]
    
    element <- parts[2]
    count <- parts[3]
    
    if (count == "") {
      count <- 1
    } else {
      count <- as.numeric(count)
    }
    
    mw <- mw + atomic_weights[element] * count
  }
  
  return(mw)
  
}