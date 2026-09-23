function(input, output, session) {
  
    session$onSessionEnded(stopApp)
  
    consistencyFlag <- reactiveVal(0) # Number (1 - no core, 2 - too many core, 0 - one core found)
    exceptionFlag <- reactiveVal(0) # Number (1 - potential core but not from data, 0 - no core found )
    
    PM <- reactiveVal(2)     # Number
    
    
    M1 <-reactiveVal(0)         # m/z of product ion associated with R1
    I1 <- reactiveVal("unknown") # ion SMILES associated with R1
    R1 <- reactiveVal("unknown")  # Neutral SMILES associated with R1
  
    molI1 <-reactiveVal(NULL)    # CDK molecule of I1
    fI1 <-reactiveVal("unknown") # FORMULA associated with ion I1
    mI1 <- reactiveVal(0)        # mass associated with ion I1 
    
    
    M2 <-reactiveVal(0)         # m/z of product ion associated with R2
    I2 <- reactiveVal("unknown") # ion SMILES associated with R2
    R2 <- reactiveVal("unknown")  # Neutral SMILES associated with R2
    
    molI2 <-reactiveVal(NULL)    # CDK molecule of I2
    fI2 <-reactiveVal("unknown") # FORMULA associated with ion I2
    mI2 <- reactiveVal(0)        # mass associated with ion I2
    
    
    M3 <-reactiveVal(0)      # mass computed by subtracting mass of ionR1, ionR2 from protonated molecule
    R3 <- reactiveVal("unknown")  # Neutral SMILES associated with core + R3
    message_R3 <- reactiveVal("")
    
    Full_Mol <- reactiveVal(NULL)     # CDK molecule of paste0(R2,R3,R1)
    Full_mf <- reactiveVal(NULL)      # molecular formula of full molecule
    Full_smiles <- reactiveVal(NULL)  # smiles for full molecule
    
    MSMS_loaded <- reactiveVal(0)
    MSMS <- reactiveVal(NA)
    
  observeEvent(input$PreMZ,{
    PM(round(input$PreMZ))
  })
  
  observeEvent(input$MSMS,{
    
    consistencyFlag(0) # Number (1 - no core, 2 - too many core, 0 - one core found)
    exceptionFlag (0) # Number (1 - potential core but not from data, 0 - no core found )
    
    M1(0)         # m/z of product ion associated with R1
    M2(0)         # m/z of product ion associated with R2
    M3(0)         # mass of core + R3
    
    I1("unknown") # ion SMILES associated with R1
    I2("unknown") # ion SMILES associated with R2
    
    R1("unknown")  # Neutral SMILES associated with R1
    R2("unknown")  # Neutral SMILES associated with R2
    R3("unknown")  # Neutral SMILES associated with core + R3
    
    fI1("unknown") # FORMULA associated with ion I1
    fI2("unknown") # FORMULA associated with ion I2
    
    mI1(0)        # mass associated with ion I1 
    mI2(0)        # mass associated with ion I2
    
    molI1(NULL)    # CDK molecule of I1
    molI2(NULL)    # CDK molecule of I2
    
    
    message_R3("")
    Full_Mol(NULL)     # CDK molecule of paste0(R2,R3,R1)
    Full_mf(NULL)      # molecular formula of full molecule
    Full_smiles(NULL)  # smiles for full molecule
    
    MSMS_loaded(0)
    MSMS(NA)
    
    MSMS(csv2dt(input$MSMS$datapath))
    MSMS_loaded(1)
  })
  
  output$report_out <- renderUI({
    
    output$Generic <- renderUI({
      if(PM()==2){
        message = paste0("<p>Please enter a nominal protonated molecule <i>m/z</i> value and upload an MS<sup>2</sup> mass spectrum to get started.</p><br>")
      } else {
        message = paste0("<p>Please upload an MS<sup>2</sup> mass spectrum to continue.</p><br>")
      }
      HTML(message)
    })
    
    if(MSMS_loaded()==1){
      P1=asm_DT(MSMS(),Tree1)
      M1(P1$mz)
      I1(P1$ion)
      R1(P1$R)
      
      if(I1()!="unknown"){
        molI1(parse.smiles(I1())[[1]])
        mI1(round(get.exact.mass(molI1()),0))
        fI1(get.mol2formula(molI1())@string)
      } else {
        fI1("unknown")
        mI1("unknown")
      }

      P2 = asm_DT(MSMS(),Tree2)
      M2(P2$mz)
      I2(P2$ion)
      R2(P2$R); 
      
      if(I2()!="unknown"){
        molI2(parse.smiles(I2())[[1]])
        mI2(round(get.exact.mass(molI2()),0))
        fI2(get.mol2formula(molI2())@string)
      } else {
        fI2("unknown")
        mI2("unknown")
      }
      
    }
    
    M3(asm_Subtraction(PM(),mI1(),mI2()))

    if(M3()=="unknown"){
      consistencyFlag(2)
    } else {
      a = which(expectedR3_mass==(M3()+1))
      if(length(a)==1){
        R3(expectedR3_SMILES[a])
        message_R3(expectedR3_message[a])
        consistencyFlag(0)
      } else {
        consistencyFlag(1)
      }
    }
    
    
    if(consistencyFlag()==0){
      if(message_R3()==""){
        temp = asm_custom_structure_joiner(R3(),R2(),R1())  
      } else {
        temp = asm_custom_structure_joiner(R3(),"",R1())
      }
      Full_Mol(generate.2d.coordinates(temp))
      Full_mf(get.mol2formula(Full_Mol())@string)
      Full_smiles(get.smiles(Full_Mol()))
    }
    
    
    output$status <- renderPrint({
        rvalues <- c(
          MSMS_loaded()
        )

        rvalues
    })
    
    output$Narrative <- renderUI({
      
      output$intro <- renderUI({
        HTML("<p>Using the supplied protonated molecule <i>m/z</i> and MS<sup>2</sup> spectrum, the tree-based predictions are as follows.</p>")
      })
      output$detailsI1 <- renderUI({
        if(I1()!= "unknown"){
          Message = paste0("<p>The ion at <i>m/z</i> ", mI1()," likely corresponds to a product ion with elemental composition ",formula2html(fI1()),". ")
          Message = paste0(Message,"A proposed  structure for this product ion is shown below.</p>")
          HTML(Message)
        } 
      })
      output$imageI1 <- renderUI({
        if(I1()!= "unknown"){
          svg <-asm_Mol2SVG(molI1(),width=50,height=50,zoom=10)
          div(
            style="width:100%,text-align: center;",
            HTML(svg)
          )
        } 
      })
      output$detailsI2 <- renderUI({
        if(I2()!= "unknown"){
          Message = paste0("<p>The ion at <i>m/z</i> ", mI2()," likely corresponds to a product ion with elemental composition ",formula2html(fI2()),". ")
          Message = paste0(Message,"A proposed  structure for this product ion is shown below.</p>")
          HTML(Message)
        } 
      })
      output$imageI2 <- renderUI({
        if(I2()!= "unknown"){
          svg <-asm_Mol2SVG(molI2(),width=50,height=50,zoom=10)
          div(
            style="width:100%,text-align: center;",
            HTML(svg)
          )
        } 
      })
      
      output$fail <- renderUI({
        return(HTML("<p>Based on the provided protonated molecule and mass spectra, the decision trees were unable to predict any information about the unknown.</p>"))
      })
      
      output$detailsR3 <- renderUI({
        if(I1()== "unknown" | I2()== "unknown") {
          Message = "<p>No further information about the unknown was predicted.</p>"
          return(HTML(Message))
        }
        Message = "<p>"
        if(consistencyFlag()==2){
            Message = paste0(Message, " In our dataset, we have not observed a nitazene analog that would be consistent with your supplied protonated molecule and the other predictions. It is possible that either our predicted product ions for the R<sub>1</sub> or R<sub>2</sub> group are incorrect; historically, the R<sub>1</sub> group product ion predictions are more accurate. ")
        } else if (consistencyFlag()==1){
          Message = paste0(Message,"<p>No further information about the unknown was predicted.</p>")
        } else if (consistencyFlag()==0){
          if(message_R3() ==""){
            Message = paste0(Message, Message = paste0(Message, "Based on our previously observed data, a possible molecular formula for your unknown is ", formula2html(Full_mf()),", with a potential structure shown below."))
          } else {
            Message = paste0(Message, message_R3())
          }
          
        }

        
          Message = paste0(Message,"</p>")

        return(HTML(Message))
      })
      
      output$imageFull <- renderUI({
        if(consistencyFlag()==0){
        svg <-asm_Mol2SVG(Full_Mol(),width=100,height=100,zoom=10)
        div(
          style="width:100%,text-align: center;",
          HTML(svg)
        )
        }
      })
      
      output$imageSMILES <-renderUI({
        HTML(paste0("<p>The SMILES associated with this predicted structure is:<br> ",Full_smiles(),"</p>"))
      })
      
      output$outro <- renderUI({
        Message = paste0("<p><b>Please note that all predictions are based on models developed with emerging data; the authors do not make any guarantees.</b><br><br></p>")
        HTML(Message)
      })
      
      
      if(I1()!= "unknown" && I2()!= "unknown"){
          if(consistencyFlag()==0){
            return(tagList(
              htmlOutput("intro"),
              htmlOutput("detailsI1"),
              htmlOutput("imageI1"),
              htmlOutput("detailsI2"),
              htmlOutput("imageI2"),
              htmlOutput("detailsR3"),
              htmlOutput("imageFull"),
              htmlOutput("imageSMILES"),
              htmlOutput("outro")
            ))  
          } else {
            return(tagList(
              htmlOutput("intro"),
              htmlOutput("detailsI1"),
              htmlOutput("imageI1"),
              htmlOutput("detailsI2"),
              htmlOutput("imageI2"),
              htmlOutput("detailsR3"),
              htmlOutput("outro")
            ))
          }
      } else if (I1()!= "unknown" && I2()=="unknown") {
            return(tagList(
              htmlOutput("intro"),
              htmlOutput("detailsI1"),
              htmlOutput("imageI1"),
              htmlOutput("detailsR3"),
              htmlOutput("outro")
            )) 
      } else if (I1()== "unknown" && I2()!="unknown") {
            return(tagList(
              htmlOutput("intro"),
              htmlOutput("detailsI2"),
              htmlOutput("imageI2"),
              htmlOutput("detailsR3"),
              htmlOutput("outro")
            )) 
      } else {
         return(tagList(
              htmlOutput("fail"),
              htmlOutput("outro")
         ))
      } 
      
    })
    
    
    if(MSMS_loaded()==1){
        return(tagList(
          uiOutput("Narrative")
          #verbatimTextOutput("status")
        ))
    } else {
      return(tagList(
        uiOutput("Generic")
        #verbatimTextOutput("status")
      ))
    }

  })
  
  
  

}
