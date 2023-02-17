
#'This file contains the script for the interactive online application (ShinyApp)
#'which generates results using the AHHME-B model


# loading packages --------------------------------------------------------

library("shiny")
library("rsconnect")
library("tidyr")
library("shinydashboard")
library("PKI")
library("packrat")
library("data.table")
library("readxl")
library("stargazer")
library("tidyverse")
library("tseries") 
library("forecast") 
library("dynlm") 
library("seastests")
library("forecast")
library("TSA")
library("epiR")
library("extraDistr")
library("MonoInc")
library("pksensi")
library("sensitivity")
library("xlsx")
library("gridExtra")
library("ggplot2")
library("reshape2")
library("here")
library("devtools")
library("multisensi")
library("rsq")
library("forcats")
library("patchwork")
library("scales")
library("ggpubr")



# The UI ------------------------------------------------------------------


#creating the tabs in the dashboard

ui <- dashboardPage(
  dashboardHeader(title = "Interactive AHHME-B Model Tool"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("How to Use This Tool", tabName = "howto", icon = icon("question")),
      #will have a guide on how to use the app
      menuItem("Input Spreadsheet", tabName = "inputs", icon = icon("file-csv")), 
      #will have the original spreadsheet, and a place to upload one
      menuItem("Main Model Outputs", tabName = "outputs-main", icon = icon("file-contract")),
      #will show the twin bar graph (Main Outputs.jpg) and the numeric results (Main Outputs.xlsx)
      menuItem("Sensitivity to Parameters (tornado plot)", tabName = "tornado", icon = icon("wind"))
      #will show the tornados (Tornado Plots.jpeg)
    )
  ),
  dashboardBody(
    tabItems(
      

# 'how to' tab  -----------------------------------------------------------

      tabItem("howto",
              box(
                title = "How to Use This Tool",
                background = "green",
                width = 12,
                textOutput("introtext")
              ),
      ),
      

# inputs tab --------------------------------------------------------------

      tabItem("inputs",
              box(
                title = "Download Link for Input Spreadsheet (default parameters)",
                background = "green",
                width = 6,
                textOutput("inputslink")
              ),
              
              box(
                title = "Upload Your Own Inputs",
                background = "green",
                width = 6,
                fileInput("inputsheet", "Upload CSV File",
                          multiple = F,
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv"))
              ),
              
              box(
                title = "Uploaded Inputs Sheet",
                background = "green",
                width = 6,
                tableOutput("inputsheettable")
              ),
              
              box(
                title = "Inputs Sheet Filepath",
                background = "green",
                width = 6,
                tableOutput("inputsheetpath")
              ),
              
              box(
                title = "Inputs Sheet",
                background = "green",
                width = 12,
                tableOutput("table")
              )
      ),
      
      

# main outputs tab --------------------------------------------------------

      tabItem("outputs-main",
              
              box(
                title = "Bar Chart Showing the Burden of AMR by Sector",
                width = 12,
                plotOutput("chart", width = "100%", height = "400px"))
              ),
      
      

# tornado plot tab --------------------------------------------------------
      
      tabItem("tornado",
              box(
                title = "Tornado Plots",
                background = "green",
                width = 12,
                plotOutput("tornado", width = "100%", height = "400px"))
      )
    )
  )
)



# server function ---------------------------------------------------------

server <- function(input, output){

  #load the AHHME-B functions
  source("Functions.R")
  
# intro text --------------------------------------------------------------


  introtext <- 
    
    "Hello and welcome to the AHHME-B interactive modelling tool. AHHME-B stands for
  'Agricultrue Human Health Micro-Economic Burden' and is a model which estimates the 
  burden of antimicrobial resistance (AMR) at the national level. It does this using an
  agent-based state transition model which allows both people and livestock animals to
  transition between health states. The burden consists of: the effct of AMR on labour
  productivity, the effect of AMR on livestock productivity, the treatment cost of 
  antimicrobial-resistant infections to the healthcare system, and the loss of human 
  (quality-adjusted) life years (QALYs) to resistant infections via morbidity and mortality.
  The model can use two approaches to calculating productivity loss from illness,
  namely the human capital and friction cost approaches, and estimates both the attributable
  and associated burdens of resistance, giving a total of four possible interpretations
  of the AMR burden.
  
  In order to use this model, navigate to the 'Input Spreadsheet' tab and download 
  the input parameter spreadsheet. You can alter the input parameters to reflect the 
  economic and epidemiological context of your own country, as well as your own methodological
  assumptions (e.g. the discount rate, the willingness to pay per QALY, or the time 
  horizon). After having done this, you can upload your input parameter spreadsheet to
  the 'Input Parameter' tab, and the app will calculate the AMR burden accordingly. 
  From here, navigate to the 'Main Model Outputs' tab, where you will find a bar
  chart,  displaying the the AMR burden decomposed by sector. The chart 
  will respectively display the attributable and associated burden under both the 
  human capital and friction cost approaches to estimating productivity loss from
  disease. From here, you can navigate to the 'Sensitivity to Parameters (tornado plot)'
  tab, which displays (under each of the four sets of assumptions) the sensitivity of
  the AMR burden to one-way variation in each of a set of key parameters, allowing you
  to gain insight into the parameters which have the greatest influence.
  
  For access to the code used in this model, please visit the GitHub repository at
  https://github.com/Trescovia/national_AMR_burden_calculator/"
  
  
  
  #create the link to the inputs spreadsheet
  inputslink <- "https://drive.google.com/file/d/1z5sUTgVDlSPMGdAiSxpwhlk24UVm2luW/view?usp=share_link"
  

# create a meta-function which produce the output plot --------------------

  shinyplotfunction <- function(inputsheet){
    
    inputsheet <- read.csv(file = input$inputsheet$datapath)
    inputsheet <- as.data.table(inputsheet)
    colnames(inputsheet) <- c("parameter", "description", "Value", "country min", "country max",
                              "Distribution", "country param 1", "country param 2")
    
    outputplot <- plot_function(inputsheet)
    
    return(outputplot)
    
  }

# create a meta-function which produces tornado plots ---------------------

  shinytornadofunction <- function(inputsheet){
    
    inputsheet <- read.csv(file = input$inputsheet$datapath)
    inputsheet <- as.data.table(inputsheet)
    colnames(inputsheet) <- c("parameter", "description", "Value", "country min", "country max",
                              "Distribution", "country param 1", "country param 2")
    
    tornadoplot <- tornado_function(inputsheet)
    
    return(tornadoplot)
    
  }

# render outputs for the UI -----------------------------------------------

  output$introtext <- renderText(introtext)
  
  output$inputslink <- renderText(inputslink)
  
  #output$chart <- renderPlot(plotfunction(inputsheet))
  
  #output$tornado <- renderPlot(tornadofunction(inputsheet))
  
  output$chart <- renderPlot({
    if(is.null(input$inputsheet)){return()}
    shinyplotfunction(inputsheet)
  })
  
  output$tornado <- renderPlot({
    if(is.null(input$inputsheet)){return()}
    shinytornadofunction(inputsheet)
  })
  
  output$inputsheettable <- renderTable({
    if(is.null(input$inputsheet)){return()}
    input$inputsheet
  })
  
  output$inputsheetpath <- renderTable({
    if(is.null(input$inputsheet)){return()}
    input$inputsheet$datapath
  })
  
  output$table <- renderTable({
    if(is.null(input$inputsheet)){return()}
    read.csv(file = input$inputsheet$datapath)
  })
  
}      
      
shinyApp(ui, server)      
      
      

    
      
      

      
      
      
