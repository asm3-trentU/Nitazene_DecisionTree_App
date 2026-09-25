This repository contains source code for the decision-tree web app for Nitazene Classification developed by Hardwick, Moorthy and Davidson (2026). 

A manuscript describing the research behind this application is forthcoming. 
A video showcasing use of the application is forthcoming.

The source code is maintained by Dr. Arun Moorthy; contact him at arunmoorthy@trentu.ca if you have any questions.

### How to work with the app
In order to use either the command line or shiny app version of the program, you need to have R installed on your computer. You also need to have a system appropriate version of Java Development Kit (JDK) installed. You can download your appropriate JDK from <a href="https://www.oracle.com/java/technologies/downloads/" target="_blank" rel="noopener noreferrer">here</a> and follow the installation instructions.

To use the command line version of the program, just run "runCLI.R" in your R console or RStudio--you can specify your protonated molecule *m/z* and the path to your two-column mass spectral data (as a .csv) on lines 17 and 18.

To use the program as a shiny app, run the "runApp.R"  code; this will automatically open the application in your default web browser. In the app, you just enter your protonated molecule *m/z* and upload your 2-column mass spec data (as a .csv). 
Example mass spectral data is available in the "ExampleData/" folder.

### How to work with the source code

There are several ways you can modify the utility of the application. Most conveniently will be updating your decision trees. To do this, you will need to click into the "Shiny_App>Data>DECISION_TREES" directory and update the .csv files. 
Alternatively, you might want to update behaviour of the actual source code. A walk-through video describing each part of the code is forthcoming. 


### Copyright and permissions

This software is distributed under an MIT license.
