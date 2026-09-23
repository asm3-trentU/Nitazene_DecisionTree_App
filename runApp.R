# The purpose of this code is to load any packages the source code depends on
# Also to clear working variables and run the shiny app

rm(list=ls())
shiny::runApp('./Shiny_App',port=7777,launch.browser="TRUE")


