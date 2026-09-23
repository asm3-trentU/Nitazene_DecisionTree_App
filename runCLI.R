# Emma's ESI-MS/MS decision trees for Nitazene analog classification
#
# Example decision trees encoded by Emma Hardwick; software implemented by Arun Moorthy
# ==============================================================================

## SETUP =======================================================================
rm(list=ls())
setwd("Shiny_App")
source("setup.R")
setwd("..")

consistencyFlag = 0
exceptionFlag = 0
## =============================================================================

## USER PROVIDED DATA ==========================================================
PM = 411
MSMS = csv2dt("ExampleData/S38.csv")
## =============================================================================

## STEPS =======================================================================
## USE QUERY MASS SPECTRA AND DECISION TREES TO PREDICT:
## 1. Protonated Molecule (from MS1)
## 2. R1 Modification (from MS2 and decision tree)
## 3. R2 Modification (from MS2 and decision tree)
## 4. R3 Mass (from results of 1, 2, and 3)
## 5. R3 Modification (from look up table of common R3 masses)
## =============================================================================


## MAKE PREDICTIONS

## using Tree 1
P1 = asm_DT(MSMS,Tree1)
M1 = P1$mz
I1 = P1$ion
R1 = P1$R

## using Tree 2
P2 = asm_DT(MSMS,Tree2)
M2 = P2$mz
I2 = P2$ion
R2 = P2$R

if(I1!="unknown"){
  molI1 = parse.smiles(I1)[[1]]
  mI1 = round(get.exact.mass(molI1),0)
  fI1 = get.mol2formula(molI1)@string
  iI1 = view.image.2d(molI1)
} else {
  fI1 = "unknown"
  mI1 = "unknown"
  iI1 = NULL
}

if(I2!="unknown"){
  molI2 = parse.smiles(I2)[[1]]
  mI2 = round(get.exact.mass(molI2),0)
  fI2 = get.mol2formula(molI2)@string
  iI2 = view.image.2d(molI2)
} else {
  fI2 = "unknown"
  mI2 = "unknown"
  iI2 = NULL
}

Mass3 = asm_Subtraction(PM,mI1,mI2)

if(Mass3=="unknown"){
  fR3 = "unknown"
  iR3 = "unknown"
  consistencyFlag = 2
} else {
  a = which(expectedR3_mass==(Mass3+1))
  if(length(a)==1){
    R3 = expectedR3_SMILES[a]
    message_R3 = expectedR3_message[a]
    consistencyFlag = 0
  } else {
    consistencyFlag = 1
  }
}

if(consistencyFlag==0){
  if(message_R3==""){
    Full_Mol = asm_custom_structure_joiner(R3,R2,R1)  
  } else {
    Full_Mol = asm_custom_structure_joiner(R3,"",R1)
  }
  
  
  generate.2d.coordinates(Full_Mol)
  # writeLines(asm_Mol2SVG(Full_Mol),"CLI_generated.svg")
  # write.molecules(Full_Mol,"FullMolecule.sdf")
  Full_mf = get.mol2formula(Full_Mol)@string
  Full_smiles = get.smiles(Full_Mol)
}

Message = "DETAILS:\n\n"
Message = paste0(Message,"The provided protonated molecule m/z values is ",PM,".\n\n")

if(I1=="unknown" && I2=="unknown"){
  Message = paste0(Message,"Using the MS2 spectra, the supplied decision trees do not make any predictions about the structure of the unknown.")
} else if(I1!="unknown" && I2=="unknonwn"){
  Message = paste0(Message,"Using the MS2 spectra and the supplied decision trees, the ion at m/z ", mI1," likely corresponds to a product ion with an elemental composition of ", fI1,".")
} else if (I1=="unknown" && I2!="unknonwn") { 
  Message = paste0(Message,"Using the MS2 spectra and the supplied decision trees, the ion at m/z ", mI2," likely corresponds to a product ion with an elemental composition of ", fI2,".")
} else {
  Message = paste0(Message, "Using the MS2 spectrum and supplied decision trees, the ion at m/z ", mI1," likely corresponds to a product ion with an elemental composition of ", fI1)
  Message = paste0(Message, ", and the ion at m/z ",mI2," likely corresponds to a product ion with an elemental composition of ",fI2,".\n\n")
  if(consistencyFlag==0){
  if(message_R3 ==""){
    Message = paste0(Message, "Based on our previously observed data, a possible molecular formula for your unknown is ", Full_mf,", and a potential SMILES is", Full_smiles,".")
  } else {
    Message = paste0(Message, message_R3)
  }
  }
}

Message = paste0(Message,"\n\nPlease note that all predictions are based on models developed with emerging data; the authors do not make any guarantees.\n\n")


cat(Message)


