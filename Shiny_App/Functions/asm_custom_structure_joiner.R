asm_custom_structure_joiner<-function(R3,R2,R1){
 
  if(R2==""){
    molR3 = generate.2d.coordinates(parse.smiles(R3)[[1]])
    molR1 = generate.2d.coordinates(parse.smiles(R1)[[1]])
    
    filename = "StrucData.sdf"
    write.molecules(c(molR3,molR1),filename,together=TRUE)
    
    strucData = readLines(filename)
    unlink(filename)
    
    breaks = which(strucData == "$$$$")
    
    ## Dealing with R3
    struc_data = strucData[1:breaks[1]]
    
    temp = strsplit(struc_data[4],' ')[[1]]
    temp = temp[temp!=""]
    natoms = as.numeric(temp[1])
    nconnections = as.numeric(temp[2])
    s1 = 5;
    s2 = s1+natoms-1
    atomInfo = seq(s1,s2)
    
    s3 = s2+1
    s4 = s2+nconnections
    connInfo = seq(s3,s4)
    
    chargeInfo = grep("CHG",struc_data)
    
    DeuteratedInfoLine = grep("ISO",struc_data)
    if(length(DeuteratedInfoLine)>0){
      DeuteratedInfo = strsplit(struc_data[DeuteratedInfoLine]," ")[[1]]
      spaces = which(DeuteratedInfo=="")
      DeuteratedInfo = DeuteratedInfo[-spaces]
      DeuteratedSet = DeuteratedInfo[seq(4,length(DeuteratedInfo),2)]
    } else {
      DeuteratedSet = NULL
    }
    
    x3 = numeric(natoms)
    y3 = numeric(natoms)
    c3 = character(natoms)
    ch3 = numeric(natoms)
    
    p3 = numeric(nconnections)
    q3 = numeric(nconnections)
    t3 = character(nconnections)
    
    for(i in 1:natoms){
      temp2 = strsplit(struc_data[atomInfo[i]]," ")[[1]]
      spaces = which(temp2=="");
      temp2 = temp2[-spaces]
      x3[i] = as.numeric(temp2[1])
      y3[i] = as.numeric(temp2[2])
      ch3[i] = as.numeric(temp2[6])
      if(i %in% DeuteratedSet){
        if(temp2[4]=="H"){
          c3[i] = "D"
        }
      } else {
        c3[i] = temp2[4]  
      }
    }
    
    for(i in 1:nconnections){
      temp3 = strsplit(struc_data[connInfo[i]]," ")[[1]]
      spaces = which(temp3=="");
      temp3 = temp3[-spaces]
      p3[i] = as.numeric(temp3[1])
      q3[i] = as.numeric(temp3[2])
      t3[i] = temp3[3]
    }
    
    if(length(chargeInfo)>0){
      for(i in 1:length(chargeInfo)){
        temp4 = strsplit(struc_data[chargeInfo[i]]," ")[[1]]
        spaces = which(temp4=="");
        temp4 = temp4[-spaces]
        a = as.numeric(temp4[4])
        if(as.numeric(temp4[5])==1){
          symbol = "+"
        } else if (as.numeric(temp4[5])==-1) {
          symbol = "-"
        }
        c3[a] = paste0(c3[a],"[",symbol,"]")
      }
    }
    y3 = -y3
    nconnections3 = nconnections
    
    
    ## Dealing with R1
    struc_data = strucData[(breaks[1]+1):breaks[2]]
    
    temp = strsplit(struc_data[4],' ')[[1]]
    temp = temp[temp!=""]
    natoms = as.numeric(temp[1])
    nconnections = as.numeric(temp[2])
    s1 = 5;
    s2 = s1+natoms-1
    atomInfo = seq(s1,s2)
    
    s3 = s2+1
    s4 = s2+nconnections
    connInfo = seq(s3,s4)
    
    chargeInfo = grep("CHG",struc_data)
    
    DeuteratedInfoLine = grep("ISO",struc_data)
    if(length(DeuteratedInfoLine)>0){
      DeuteratedInfo = strsplit(struc_data[DeuteratedInfoLine]," ")[[1]]
      spaces = which(DeuteratedInfo=="")
      DeuteratedInfo = DeuteratedInfo[-spaces]
      DeuteratedSet = DeuteratedInfo[seq(4,length(DeuteratedInfo),2)]
    } else {
      DeuteratedSet = NULL
    }
    
    x2 = numeric(natoms)
    y2 = numeric(natoms)
    c2 = character(natoms)
    ch2 = numeric(natoms)
    
    p2 = numeric(nconnections)
    q2 = numeric(nconnections)
    t2 = character(nconnections)
    
    for(i in 1:natoms){
      temp2 = strsplit(struc_data[atomInfo[i]]," ")[[1]]
      spaces = which(temp2=="");
      temp2 = temp2[-spaces]
      x2[i] = as.numeric(temp2[1])
      y2[i] = as.numeric(temp2[2])
      ch2[i] = as.numeric(temp2[6])
      if(i %in% DeuteratedSet){
        if(temp2[4]=="H"){
          c2[i] = "D"
        }
      } else {
        c2[i] = temp2[4]  
      }
    }
    
    for(i in 1:nconnections){
      temp3 = strsplit(struc_data[connInfo[i]]," ")[[1]]
      spaces = which(temp3=="");
      temp3 = temp3[-spaces]
      p2[i] = as.numeric(temp3[1])
      q2[i] = as.numeric(temp3[2])
      t2[i] = temp3[3]
    }
    
    if(length(chargeInfo)>0){
      for(i in 1:length(chargeInfo)){
        temp4 = strsplit(struc_data[chargeInfo[i]]," ")[[1]]
        spaces = which(temp4=="");
        temp4 = temp4[-spaces]
        a = as.numeric(temp4[4])
        if(as.numeric(temp4[5])==1){
          symbol = "+"
        } else if (as.numeric(temp4[5])==-1) {
          symbol = "-"
        }
        c2[a] = paste0(c2[a],"[",symbol,"]")
      }
    }
    x2 = -x2
    y2 = y2
    p2 = p2 + max(p3,q3)
    q2 = q2 + max(p3,q3)
    nconnections2 = nconnections
    
    ## shift R1 so that the question marks are aligned
    deltaX = x3[which(c3=="?")]-x2[which(c2=="?")]
    deltaY = y3[which(c3=="?")]-y2[which(c2=="?")]
    x2 = x2+deltaX
    y2 = y2+deltaY
    
    
    x = c(x3,x2)
    y = c(y3,y2)
    c = c(c3,c2)
    ch = c(ch3,ch2)
    p = c(p3,p2)
    q = c(q3,q2)
    t = c(t3,t2)
    
    q_index = which(c=="?")
    a = which(p==q_index[2])
    if(length(a)>0){
      p[a] = q_index[1]
    }
    a = which(q==q_index[2])
    if(length(a)>0){
      q[a] = q_index[1]
    }
    
    b = q_index[2]
    for(i in b){
      a = which(p>i)
      p[a] = p[a]-1
      a = which(q>i)
      q[a] = q[a]-1
    }
    
    c = c[-b]
    x = x[-b]
    y = y[-b]
    ch = ch[-b]
    
    # c = gsub("!","C",c)
    c = gsub("\\?","N",c)
    
    if(FALSE){
      minX = min(x); maxX = max(x); midX = 0.5*(minX+maxX)
      minY = min(y); maxY = max(y); midY = 0.5*(minY+maxY)
      
      x = x-midX;
      y = y-midY;
      
      plot(x,y,
           xlim=c(min(x),max(x)),
           ylim=c(min(y),max(y)),
           axes=TRUE,
           cex=0.1,
           pch=18,
           col = "white",
           xlab="",ylab="")
      
      offset=0.006*(maxX-minX)
      offset3 = 1.1*offset
      
      for(i in 1:(nconnections3+nconnections2+nconnections1)){
        if(t[i]=="1"){
          segments(x[p[i]],y[p[i]],x[q[i]],y[q[i]])
        } else if (t[i] =="2"){
          if(x[p[i]] == x[q[i]]){ # vertical
            segments((x[p[i]]-offset3),y[p[i]],(x[q[i]]-offset3),y[q[i]],col="blue")
            segments((x[p[i]]+offset3),y[p[i]],(x[q[i]]+offset3),y[q[i]],col="blue")
          } else {
            segments(x[p[i]],(y[p[i]]-offset),x[q[i]],(y[q[i]]-offset),col="blue")
            segments(x[p[i]],(y[p[i]]+offset),x[q[i]],(y[q[i]]+offset),col="blue")
          }
        } else if (t[i] == "3"){
          if(x[p[i]] == x[q[i]]){
            segments((x[p[i]]-offset3),y[p[i]],(x[q[i]]-offset3),y[q[i]],col="red")
            segments(x[p[i]],y[p[i]],x[q[i]],y[q[i]],col="red")
            segments((x[p[i]]+offset3),y[p[i]],(x[q[i]]+offset3),y[q[i]],col="red")
          } else {
            segments(x[p[i]],(y[p[i]]-offset3),x[q[i]],(y[q[i]]-offset3),col="red")
            segments(x[p[i]],y[p[i]],x[q[i]],y[q[i]],col="red")
            segments(x[p[i]],(y[p[i]]+offset3),x[q[i]],(y[q[i]]+offset3),col="red")
          }
        }
        
      }
      nonC = which(c!="C")
      matplot(x[nonC],y[nonC],col = "white",pch=15,
              cex = 1.5,add=TRUE)
      
      text(x[nonC],y[nonC],c[nonC],cex=0.7)
    }
    
    
    chargeMols = grep("\\[",c)
    if(length(chargeMols)>0){
      chargeInfo = numeric(length(chargeMols))
      for(i in chargeMols){
        temp = strsplit(c[i],"\\[")[[1]]
        c[i] = temp[1]
        temp2 = strsplit(temp[2],"\\]")[[1]][1]
        chargeInfo[i] = as.numeric(paste0(temp2,1))
      }  
    }
    
    filename = "custom.sdf"
    sink(filename)
    cat("\n")
    cat(" aruns nitazine construction\n")
    cat("\n")
    cat(paste0(" ",length(x)," ",length(p)," 0 0 0 0 0 0 0 0999 V2000\n"))
    for(i in 1:length(x)){
      cat(paste0("    0.0000    0.0000    0.0000 ",c[i],"   0 ",sprintf("%3d",ch[i]),"  0  0  0  0  0  0  0  0  0  0\n"))
    }
    for(i in 1:length(p)){
      cat(paste0(sprintf("%3d",p[i]),sprintf("%3d",q[i]),sprintf("%3s",t[i]),"  0  0  0  0\n"))
    }
    if(length(chargeInfo)>0){
      for(i in 1:length(chargeInfo)){
        cat(paste0("M  CHG  1",sprintf("%4s",chargeMols[i]),sprintf("%4s",chargeInfo[i]),"\n"))
      }
    }
    cat("M  END\n")
    cat("$$$$\n")
    sink()
    
    MOL = load.molecules(filename)[[1]]
    unlink(filename)
    return(MOL)  
  }
  

    molR3 = generate.2d.coordinates(parse.smiles(R3)[[1]])
    molR2 = generate.2d.coordinates(parse.smiles(R2)[[1]])
    molR1 = generate.2d.coordinates(parse.smiles(R1)[[1]])
    
    filename = "StrucData.sdf"
    write.molecules(c(molR3,molR2,molR1),filename,together=TRUE)
    
    strucData = readLines(filename)
    unlink(filename)
    
    breaks = which(strucData == "$$$$")
    
    ## Dealing with R3
    struc_data = strucData[1:breaks[1]]
    
    temp = strsplit(struc_data[4],' ')[[1]]
    temp = temp[temp!=""]
    natoms = as.numeric(temp[1])
    nconnections = as.numeric(temp[2])
    s1 = 5;
    s2 = s1+natoms-1
    atomInfo = seq(s1,s2)
    
    s3 = s2+1
    s4 = s2+nconnections
    connInfo = seq(s3,s4)
    
    chargeInfo = grep("CHG",struc_data)
    
    DeuteratedInfoLine = grep("ISO",struc_data)
    if(length(DeuteratedInfoLine)>0){
      DeuteratedInfo = strsplit(struc_data[DeuteratedInfoLine]," ")[[1]]
      spaces = which(DeuteratedInfo=="")
      DeuteratedInfo = DeuteratedInfo[-spaces]
      DeuteratedSet = DeuteratedInfo[seq(4,length(DeuteratedInfo),2)]
    } else {
      DeuteratedSet = NULL
    }
    
    x3 = numeric(natoms)
    y3 = numeric(natoms)
    c3 = character(natoms)
    ch3 = numeric(natoms)
    
    p3 = numeric(nconnections)
    q3 = numeric(nconnections)
    t3 = character(nconnections)
    
    for(i in 1:natoms){
      temp2 = strsplit(struc_data[atomInfo[i]]," ")[[1]]
      spaces = which(temp2=="");
      temp2 = temp2[-spaces]
      x3[i] = as.numeric(temp2[1])
      y3[i] = as.numeric(temp2[2])
      ch3[i] = as.numeric(temp2[6])
      if(i %in% DeuteratedSet){
        if(temp2[4]=="H"){
          c3[i] = "D"
        }
      } else {
        c3[i] = temp2[4]  
      }
    }
    
    for(i in 1:nconnections){
      temp3 = strsplit(struc_data[connInfo[i]]," ")[[1]]
      spaces = which(temp3=="");
      temp3 = temp3[-spaces]
      p3[i] = as.numeric(temp3[1])
      q3[i] = as.numeric(temp3[2])
      t3[i] = temp3[3]
    }
    
    if(length(chargeInfo)>0){
      for(i in 1:length(chargeInfo)){
        temp4 = strsplit(struc_data[chargeInfo[i]]," ")[[1]]
        spaces = which(temp4=="");
        temp4 = temp4[-spaces]
        a = as.numeric(temp4[4])
        if(as.numeric(temp4[5])==1){
          symbol = "+"
        } else if (as.numeric(temp4[5])==-1) {
          symbol = "-"
        }
        c3[a] = paste0(c3[a],"[",symbol,"]")
      }
    }
    y3 = -y3
    nconnections3 = nconnections
    
    
    ## Dealing with R2
    struc_data = strucData[(breaks[1]+1):breaks[2]]
    
    temp = strsplit(struc_data[4],' ')[[1]]
    temp = temp[temp!=""]
    natoms = as.numeric(temp[1])
    nconnections = as.numeric(temp[2])
    s1 = 5;
    s2 = s1+natoms-1
    atomInfo = seq(s1,s2)
    
    s3 = s2+1
    s4 = s2+nconnections
    connInfo = seq(s3,s4)
    
    chargeInfo = grep("CHG",struc_data)
    
    DeuteratedInfoLine = grep("ISO",struc_data)
    if(length(DeuteratedInfoLine)>0){
      DeuteratedInfo = strsplit(struc_data[DeuteratedInfoLine]," ")[[1]]
      spaces = which(DeuteratedInfo=="")
      DeuteratedInfo = DeuteratedInfo[-spaces]
      DeuteratedSet = DeuteratedInfo[seq(4,length(DeuteratedInfo),2)]
    } else {
      DeuteratedSet = NULL
    }
    
    x2 = numeric(natoms)
    y2 = numeric(natoms)
    c2 = character(natoms)
    ch2 = numeric(natoms)
    
    p2 = numeric(nconnections)
    q2 = numeric(nconnections)
    t2 = character(nconnections)
    
    for(i in 1:natoms){
      temp2 = strsplit(struc_data[atomInfo[i]]," ")[[1]]
      spaces = which(temp2=="");
      temp2 = temp2[-spaces]
      x2[i] = as.numeric(temp2[1])
      y2[i] = as.numeric(temp2[2])
      ch2[i] = as.numeric(temp2[6])
      if(i %in% DeuteratedSet){
        if(temp2[4]=="H"){
          c2[i] = "D"
        }
      } else {
        c2[i] = temp2[4]  
      }
    }
    
    for(i in 1:nconnections){
      temp3 = strsplit(struc_data[connInfo[i]]," ")[[1]]
      spaces = which(temp3=="");
      temp3 = temp3[-spaces]
      p2[i] = as.numeric(temp3[1])
      q2[i] = as.numeric(temp3[2])
      t2[i] = temp3[3]
    }
    
    if(length(chargeInfo)>0){
      for(i in 1:length(chargeInfo)){
        temp4 = strsplit(struc_data[chargeInfo[i]]," ")[[1]]
        spaces = which(temp4=="");
        temp4 = temp4[-spaces]
        a = as.numeric(temp4[4])
        if(as.numeric(temp4[5])==1){
          symbol = "+"
        } else if (as.numeric(temp4[5])==-1) {
          symbol = "-"
        }
        c2[a] = paste0(c2[a],"[",symbol,"]")
      }
    }
    x2 = -x2
    y2 = y2
    p2 = p2 + max(p3,q3)
    q2 = q2 + max(p3,q3)
    nconnections2 = nconnections
    
    ## shift R2 so that the exclamation marks are aligned
    deltaX = x3[which(c3=="!")]-x2[which(c2=="!")]
    deltaY = y3[which(c3=="!")]-y2[which(c2=="!")]
    x2 = x2+deltaX
    y2 = y2+deltaY
    
    ## Dealing with R1
    struc_data = strucData[(breaks[2]+1):breaks[3]]
    
    temp = strsplit(struc_data[4],' ')[[1]]
    temp = temp[temp!=""]
    natoms = as.numeric(temp[1])
    nconnections = as.numeric(temp[2])
    s1 = 5;
    s2 = s1+natoms-1
    atomInfo = seq(s1,s2)
    
    s3 = s2+1
    s4 = s2+nconnections
    connInfo = seq(s3,s4)
    
    chargeInfo = grep("CHG",struc_data)
    
    DeuteratedInfoLine = grep("ISO",struc_data)
    if(length(DeuteratedInfoLine)>0){
      DeuteratedInfo = strsplit(struc_data[DeuteratedInfoLine]," ")[[1]]
      spaces = which(DeuteratedInfo=="")
      DeuteratedInfo = DeuteratedInfo[-spaces]
      DeuteratedSet = DeuteratedInfo[seq(4,length(DeuteratedInfo),2)]
    } else {
      DeuteratedSet = NULL
    }
    
    x1 = numeric(natoms)
    y1 = numeric(natoms)
    c1 = character(natoms)
    ch1 = numeric(natoms)
    
    p1 = numeric(nconnections)
    q1 = numeric(nconnections)
    t1 = character(nconnections)
    
    for(i in 1:natoms){
      temp2 = strsplit(struc_data[atomInfo[i]]," ")[[1]]
      spaces = which(temp2=="");
      temp2 = temp2[-spaces]
      x1[i] = as.numeric(temp2[1])
      y1[i] = as.numeric(temp2[2])
      ch1[i] = as.numeric(temp2[6])
      if(i %in% DeuteratedSet){
        if(temp2[4]=="H"){
          c1[i] = "D"
        }
      } else {
        c1[i] = temp2[4]  
      }
    }
    
    for(i in 1:nconnections){
      temp3 = strsplit(struc_data[connInfo[i]]," ")[[1]]
      spaces = which(temp3=="");
      temp3 = temp3[-spaces]
      p1[i] = as.numeric(temp3[1])
      q1[i] = as.numeric(temp3[2])
      t1[i] = temp3[3]
    }
    
    if(length(chargeInfo)>0){
      for(i in 1:length(chargeInfo)){
        temp4 = strsplit(struc_data[chargeInfo[i]]," ")[[1]]
        spaces = which(temp4=="");
        temp4 = temp4[-spaces]
        a = as.numeric(temp4[4])
        if(as.numeric(temp4[5])==1){
          symbol = "+"
        } else if (as.numeric(temp4[5])==-1) {
          symbol = "-"
        }
        c1[a] = paste0(c1[a],"[",symbol,"]")
      }
    }
    x1 = -x1
    y1 = -y1
    p1 = p1 + max(p2,q2)
    q1 = q1 + max(p2,q2)
    nconnections1 = nconnections
    
    ## shift R1 so that the exclamation marks are aligned
    deltaX = x3[which(c3=="?")]-x1[which(c1=="?")]
    deltaY = y3[which(c3=="?")]-y1[which(c1=="?")]
    x1 = x1+deltaX
    y1 = y1+deltaY
    
    
    x = c(x3,x2,x1)
    y = c(y3,y2,y1)
    c = c(c3,c2,c1)
    ch = c(ch3,ch2,ch1)
    p = c(p3,p2,p1)
    q = c(q3,q2,q1)
    t = c(t3,t2,t1)
    
    e_index = which(c=="!")
    a = which(p==e_index[2])
    if(length(a)>0){
      p[a] = e_index[1]
    }
    a = which(q==e_index[2])
    if(length(a)>0){
      q[a] = e_index[1]
    }
    
    
    q_index = which(c=="?")
    a = which(p==q_index[2])
    if(length(a)>0){
      p[a] = q_index[1]
    }
    a = which(q==q_index[2])
    if(length(a)>0){
      q[a] = q_index[1]
    }
    
    b = sort(c(e_index[2],q_index[2]),decreasing = TRUE)
    for(i in b){
      a = which(p>i)
      p[a] = p[a]-1
      a = which(q>i)
      q[a] = q[a]-1
    }
    
    c = c[-b]
    x = x[-b]
    y = y[-b]
    ch = ch[-b]
    
    c = gsub("!","C",c)
    c = gsub("\\?","N",c)
    
    if(FALSE){
      minX = min(x); maxX = max(x); midX = 0.5*(minX+maxX)
      minY = min(y); maxY = max(y); midY = 0.5*(minY+maxY)
      
      x = x-midX;
      y = y-midY;
      
      plot(x,y,
           xlim=c(min(x),max(x)),
           ylim=c(min(y),max(y)),
           axes=TRUE,
           cex=0.1,
           pch=18,
           col = "white",
           xlab="",ylab="")
      
      offset=0.006*(maxX-minX)
      offset3 = 1.1*offset
      
      for(i in 1:(nconnections3+nconnections2+nconnections1)){
        if(t[i]=="1"){
          segments(x[p[i]],y[p[i]],x[q[i]],y[q[i]])
        } else if (t[i] =="2"){
          if(x[p[i]] == x[q[i]]){ # vertical
            segments((x[p[i]]-offset3),y[p[i]],(x[q[i]]-offset3),y[q[i]],col="blue")
            segments((x[p[i]]+offset3),y[p[i]],(x[q[i]]+offset3),y[q[i]],col="blue")
          } else {
            segments(x[p[i]],(y[p[i]]-offset),x[q[i]],(y[q[i]]-offset),col="blue")
            segments(x[p[i]],(y[p[i]]+offset),x[q[i]],(y[q[i]]+offset),col="blue")
          }
        } else if (t[i] == "3"){
          if(x[p[i]] == x[q[i]]){
            segments((x[p[i]]-offset3),y[p[i]],(x[q[i]]-offset3),y[q[i]],col="red")
            segments(x[p[i]],y[p[i]],x[q[i]],y[q[i]],col="red")
            segments((x[p[i]]+offset3),y[p[i]],(x[q[i]]+offset3),y[q[i]],col="red")
          } else {
            segments(x[p[i]],(y[p[i]]-offset3),x[q[i]],(y[q[i]]-offset3),col="red")
            segments(x[p[i]],y[p[i]],x[q[i]],y[q[i]],col="red")
            segments(x[p[i]],(y[p[i]]+offset3),x[q[i]],(y[q[i]]+offset3),col="red")
          }
        }
        
      }
      nonC = which(c!="C")
      matplot(x[nonC],y[nonC],col = "white",pch=15,
              cex = 1.5,add=TRUE)
      
      text(x[nonC],y[nonC],c[nonC],cex=0.7)
    }
    
    
    chargeMols = grep("\\[",c)
    if(length(chargeMols)>0){
      chargeInfo = numeric(length(chargeMols))
      for(i in chargeMols){
        temp = strsplit(c[i],"\\[")[[1]]
        c[i] = temp[1]
        temp2 = strsplit(temp[2],"\\]")[[1]][1]
        chargeInfo[i] = as.numeric(paste0(temp2,1))
      }  
    }
    
    filename = "custom.sdf"
    sink(filename)
    cat("\n")
    cat(" aruns nitazine construction\n")
    cat("\n")
    cat(paste0(" ",length(x)," ",length(p)," 0 0 0 0 0 0 0 0999 V2000\n"))
    for(i in 1:length(x)){
      cat(paste0("    0.0000    0.0000    0.0000 ",c[i],"   0 ",sprintf("%3d",ch[i]),"  0  0  0  0  0  0  0  0  0  0\n"))
    }
    for(i in 1:length(p)){
      cat(paste0(sprintf("%3d",p[i]),sprintf("%3d",q[i]),sprintf("%3s",t[i]),"  0  0  0  0\n"))
    }
    if(length(chargeInfo)>0){
      for(i in 1:length(chargeInfo)){
        cat(paste0("M  CHG  1",sprintf("%4s",chargeMols[i]),sprintf("%4s",chargeInfo[i]),"\n"))
      }
    }
    cat("M  END\n")
    cat("$$$$\n")
    sink()
    
    MOL = load.molecules(filename)[[1]]
    unlink(filename)
    return(MOL)  

}