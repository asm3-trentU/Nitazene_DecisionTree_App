functions = list.files("Functions",full.names = TRUE)
if(length(functions)>0){
  for(i in functions){
    source(i)
  }
}

expectedR3 = read.csv("Data/Expected Cores_09172026.csv")
expectedR3_SMILES = expectedR3$Expected.Cores
expectedR3_mass = expectedR3$Mass
expectedR3_message = expectedR3$Message

Tree1 = read.csv("Data/DECISION_TREES/ESI-MSMS Decision Tree R1_09222026.csv")
Tree2 = read.csv("Data/DECISION_TREES/ESI-MSMS Decision Tree R2_08272026.csv")


