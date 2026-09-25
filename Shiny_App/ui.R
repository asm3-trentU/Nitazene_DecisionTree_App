source("setup.R")
page_sidebar(
    
    tags$head(
    tags$style(HTML("

    .sidebar .form-group {
      margin-bottom: 0.5rem;
    }

    .sidebar .btn {
      padding: 0.25rem 0.5rem;
    }
    
    .shiny-input-container {
      margin-bottom: 6px;
    }
    
    /* Tab text */ 
    .nav-tabs .nav-link { 
    color: #555555; 
    } 
    
    /* Active tab */ 
    .nav-tabs .nav-link.active {
    color: #2C5F7C; 
    font-weight: 600; 
    border-bottom: 
    3px solid #2C5F7C; 
    }
    
    "))
    ),
    
    # Application title
    title="ESI-MS/MS Decision Tree Application (Nitazene)",

    
    sidebar=sidebar(
        
      width=350,
        
      h3("Analyst Input"),
      
      numericInput(
        inputId = "PreMZ",
        label = HTML("Nominal protonated molecule <i>m/z</i>"),
        value = 369,
        min = 100,
        max = 2000,
        step = 1
      ),
      
      fileInput(
          "MSMS",
          HTML("Upload MS<sup>2</sup> spectrum"),
          accept = c(".csv")
      )

    ),
    
    navset_tab(
      nav_panel(
        "Predicted Results",
        card(
          uiOutput("report_out")
        )
      ),
      nav_panel(
        "About the app",
        card(
          p(HTML("The Nitazene Decision Tree app employs encoded decision trees and a product ion spectrum to predict the potential structure of an unknown. The product ion mass spectrum should be provided as a two-column (<i>m/z</i>, intensity) .csv file  without header information. You can download an example mass spectrum "),tags$a(href="exampleMS2.csv",download="exampleMS2.csv","here"),HTML(".")),
          p("For reference, the structure predictions are based on the following scaffold:"),
          img(src="Scaffold.svg",width="50%"),
          p( HTML("<b>Note on core + R<sub>3</sub> values:</b>")),
          p(HTML("The full structure created as a part of the results is made by combining the R<sub>1</sub> substitution, the R<sub>2</sub> substitution and the core group + R<sub>3</sub> substitution. The core + R<sub>3</sub> value is calculated by subtracting the mass of the R<sub>1</sub> and R<sub>2</sub> product ions from the mass of the protonated molecule. Assuming the R<sub>1</sub> and R<sub>2</sub> predictions are correct, the core + R<sub>3</sub> value should be one of five potential outcomes based on our dataset.")),
          p(HTML("However, several compounds that we tested consistently produced incorrect R<sub>2</sub> predictions, leading to an incorrect core + R<sub>3</sub> prediction. This happened for isopropoxy/propoxy substitutions (e.g., isotonitazene, protonitazene, etc.) and butoxy substitutions (e.g., butonitazene, <i>sec</i>-butonitazene, and <i>iso</i>-butonitazene.).")),
          p(HTML("The incorrect R<sub>2</sub> prediction is attributed to the diagnostic ion appearing at very low abundances, meaning that it is occasionally not present in the centroid data, and cannot be identified by the decision tree. For example, isotonitazene should produce a product ion at <i>m/z</i> 149, but this ion was sometimes at a relative abundance below 1% in our data. That said, the incorrect core + R<sub>3</sub> values were consistent for these compounds, creating several exceptions to the five standard core + R<sub>3</sub> values.")),
          p(HTML("Therefore, a resulting full structure may include an isopropoxy/propoxy or butoxy group rather than the R<sub>2</sub> that was predicted based on the ions present in the product ion spectrum. When this is the case, the full structure was identified by one of the core + R<sub>3</sub> exceptions.")),
          p(HTML("For more information, please see the associated manuscript." ))
        )
      )
    ),
    div(
      style = "width: 100%; margin: 0 auto",
      p( "This software tool was co-developed by researchers from the Davidson Research Group at Sam Houston State University (Huntsville, TX, USA) and the CRAFTS Lab at Trent University (Peterborough, ON, Canada).", style = "font-size: 80%;"),
      img(src="Logos.png",width="30%")
    )

)

