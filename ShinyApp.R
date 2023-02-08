
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
      #will have a guideo on how to use the app
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
                textOutput(inputslink)
              ),
              box(
                title = "Upload Your Own Inputs",
                background = "green",
                width = 6,
                fileInput("inputs", "Upload CSV File",
                          multiple = F,
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv"))
              )
      ),
      
      

# main outputs tab --------------------------------------------------------

      tabItem("outputs-main",
              
              box(
                title = "Bar Chart Showing the Burden of AMR by Sector
                using the Human Captial Approach to Estimating Productivity Losses",
                width = 12,
                plotOutput(chart, width = "100%", height = "400px"))
                
              # ),
              # 
              # box(
              #   title = "Bar Chart Showing the Associated Burden of AMR (total and by sector)
              #   using the Friction Cost Approach to Estimating Productivity Losses",
              #   width = 6,
              #   plotOutput(chartAssFCA, width = "100%", height = "400px")
              # ),
              # 
              # box(
              #   title = "Bar Chart Showing the Attributable Burden of AMR (total and by sector)
              #   using the Human Captial Approach to Estimating Productivity Losses",
              #   width = 6,
              #   plotOutput(chartAttHCA, width = "100%", height = "400px")
              #   
              # ),
              # 
              # box(
              #   title = "Bar Chart Showing the Attributable Burden of AMR (total and by sector)
              #   using the Friction Cost Approach to Estimating Productivity Losses",
              #   width = 6,
              #   plotOutput(chartAttFCA, width = "100%", height = "400px")
              #   )
              ),
      
      

# tornado plot tab --------------------------------------------------------
      
      tabItem("tornado",
              box(
                title = "Tornado Plot (associated burden, human capital approach)",
                background = "green",
                width = 12,
                plotOutput(tonadoAssHCA, width = "100%", height = "400px"))
              # ),
              # 
              # box(
              #   title = "Tornado Plot (associated burden, friction cost approach)",
              #   background = "green",
              #   width = 6,
              #   plotOutput(tonadoAssFCA, width = "100%", height = "400px")
              # 
              # ),
              # 
              # box(
              #   title = "Tornado Plot (attributable burden, human capital approach)",
              #   background = "green",
              #   width = 6,
              #   plotOutput(tonadoAttHCA, width = "100%", height = "400px")
              # 
              # ),
              # 
              # box(
              #   title = "Tornado Plot (attributable burden, friction cost approach)",
              #   background = "green",
              #   width = 6,
              #   plotOutput(tonadottFCA, width = "100%", height = "400px")
              # )
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
  inputslink <- "https://drive.google.com/file/d/1VHPIm9QozNLQ9n2KC_C-NRdTtqCDyucD/view?usp=sharing"
  

# create a meta-function which produce the output plot --------------------

  plotfunction <- function(inputs){
    inputs <- read.csv("inputs")
    inputs <- as.data.table(inputs)
    colnames(inputs) <- c("parameter", "description", "Value", "country min", "country max",
                          "Distribution", "country param 1", "country param 2")
    
    
    #table (not displayed)
    main_outputs <- matrix(rep(0), ncol = 5, nrow = 11)
    colnames(main_outputs) = c("Output", "Associated Burden (HCA)", "Associated Burden (FCA)",
                               "Attributable Burden (HCA)", "Attributable Burden (FCA)")
    main_outputs[,1] <- Model(inputs, "HCA")[,1]
    main_outputs[,2] <- as.numeric(Model(inputs, "HCA")[,2])
    main_outputs[,3] <- as.numeric(Model(inputs, "FCA")[,2])
    main_outputs[,4] <- as.numeric(Model(inputs, "HCA")[,3])
    main_outputs[,5] <- as.numeric(Model(inputs, "FCA")[,3])
    
    #graphical representation
    main_outputs_df <- matrix(as.numeric(unlist(main_outputs)),nrow=nrow(main_outputs))[,2:5] #numeric
    main_outputs_df <- main_outputs_df[c(1,3,4,5,10),]
    
    ##HCA
    
    scenario_prod = "HCA"
    
    main_outputs_df_rearranged_HCA <- as.data.frame(matrix(rep(0), nrow = 10, ncol = 3))
    colnames(main_outputs_df_rearranged_HCA) = c("sector", "burden_estimation_method", "value")
    main_outputs_df_rearranged_HCA[,1] <- rep(main_outputs[c(1,3,4,5,10),1])
    main_outputs_df_rearranged_HCA$burden_estimation_method <- c(rep("Associated", 5), rep("Attributable",5))
    main_outputs_df_rearranged_HCA$value <- c(main_outputs_df[,1], main_outputs_df[,3])
    
    main_outputs_plot_HCA <- ggplot(main_outputs_df_rearranged_HCA, aes(fill = burden_estimation_method,
                                                                        y = value, x = sector)) +
      geom_bar(position = "dodge", stat = "identity") +
      ggtitle("Distribution of AMR Burden by Sector over 20 Years, Human Capital Method") +
      labs(y = "Burden (2019 $USD)", x = "Sector")
    
    #main_outputs_plot_HCA
    
    ##FCA
    main_outputs_df_rearranged_FCA <- as.data.frame(matrix(rep(0), nrow = 10, ncol = 3))
    colnames(main_outputs_df_rearranged_FCA) = c("sector", "burden_estimation_method", "value")
    main_outputs_df_rearranged_FCA[,1] <- rep(main_outputs[c(1,3,4,5,10),1])
    main_outputs_df_rearranged_FCA$burden_estimation_method <- c(rep("Associated", 5), rep("Attributable",5))
    main_outputs_df_rearranged_FCA$value <- c(main_outputs_df[,2], main_outputs_df[,4])
    
    main_outputs_plot_FCA <- ggplot(main_outputs_df_rearranged_FCA, aes(fill = burden_estimation_method,
                                                                        y = value, x = sector)) +
      geom_bar(position = "dodge", stat = "identity") +
      ggtitle("Distribution of AMR Burden by Sector over 20 Years, Friction Cost Method") +
      labs(y = "Burden (2019 $USD)", x = "Sector")
    
    #main_outputs_plot_FCA
    
    outputplot <- ggarrange(main_outputs_plot_HCA, main_outputs_plot_FCA,
                            labels = c("HCA", "FCA"),
                            ncol = 1, nrow = 2)
    
    return(outputplot)
    
  }
  


# create a meta-function which produces tornado plots ---------------------


  tornadofunction <- function(inputs){
    
    #loading files
    inputs <- read.csv("inputs")
    inputs <- as.data.table(inputs)
    colnames(inputs) <- c("parameter", "description", "Value", "country min", "country max",
                          "Distribution", "country param 1", "country param 2")
    
    #producing tornados
    HCA_base_associated <-  as.numeric(Model(inputs, "HCA")[1,2])
    HCA_base_attributable <-  as.numeric(Model(inputs, "HCA")[1,3])
    FCA_base_associated <-  as.numeric(Model(inputs, "FCA")[1,2])
    FCA_base_attributable <-  as.numeric(Model(inputs, "FCA")[1,3])
    
    inputs_tornado <- inputs
    
    #create tornado results matrix
    tornado_values <- as.data.frame(matrix(ncol = 19, nrow = 4))
    rownames(tornado_values) = c("HCA - Associated Burden", "HCA - Attributable Burden",
                                 "FCA - Associated Burden", "FCA - Attributable Burden")
    colnames(tornado_values) = c("baseline",
                                 "discount rate - low", "discount rate - high",
                                 "population growth - low", "population growth - high",
                                 "WTP threshold - low", "WTP threshold - high",
                                 "productivity growth - low", "productivity growth - high",
                                 "LFPR - low", "LFPR - high",
                                 "infection incidence - low", "infection incidence - high",
                                 "AMR incidence - low", "AMR incidence - high",
                                 "AMR growth - low", "AMR growth - high",
                                 "livestock industrialisation - low", "livestock industrialisation - high")
    
    #generate tornado values
    
    tornado_values["baseline"] <- c(HCA_base_associated, HCA_base_attributable,
                                    FCA_base_associated, FCA_base_attributable)
    
    j = 1
    
    for(i in c(3, 5:7, 10, 27:28, 48)){
      
      j = j + 1 #j lets us choose the column of tornado_values which we put data in!!!!!! 
      
      inputs_tornado <- inputs
      inputs_tornado[i,3] <- inputs_tornado[i,4]
      tornado_values[1,j] <- as.numeric(Model(inputs_tornado, "HCA")[1,2]) - tornado_values[1,1]
      tornado_values[2,j] <- as.numeric(Model(inputs_tornado, "HCA")[1,3]) - tornado_values[2,1]
      tornado_values[3,j] <- as.numeric(Model(inputs_tornado, "FCA")[1,2]) - tornado_values[3,1]
      tornado_values[4,j] <- as.numeric(Model(inputs_tornado, "FCA")[1,3]) - tornado_values[4,1]
      
      j = j+1
      
      inputs_tornado <- inputs
      inputs_tornado[i,3] <- inputs_tornado[i,5]
      tornado_values[1,j] <- as.numeric(Model(inputs_tornado, "HCA")[1,2]) - tornado_values[1,1]
      tornado_values[2,j] <- as.numeric(Model(inputs_tornado, "HCA")[1,3]) - tornado_values[2,1]
      tornado_values[3,j] <- as.numeric(Model(inputs_tornado, "FCA")[1,2]) - tornado_values[3,1]
      tornado_values[4,j] <- as.numeric(Model(inputs_tornado, "FCA")[1,3]) - tornado_values[4,1]
      
      
    }
    
    #fill in data frames for the plots
    
    tornado_HCA_associated <- data.frame("variable" = c("Discount Rate",
                                                        "Population Growth Rate",
                                                        "Willingness-to-Pay Threshold",
                                                        "Productivity Growth Rate",
                                                        "Labour Force Participation Rate",
                                                        "Incidence of Bacterial Infection",
                                                        "Rate of AMR",
                                                        "AMR Growth Rate",
                                                        "Livestock Industrialisation"),
                                         
                                         "min" = c(tornado_values[1,2],
                                                   tornado_values[1,4],
                                                   tornado_values[1,6],
                                                   tornado_values[1,8],
                                                   tornado_values[1,10],
                                                   tornado_values[1,12],
                                                   tornado_values[1,14],
                                                   tornado_values[1,16],
                                                   tornado_values[1,18]),
                                         
                                         "max" = c(tornado_values[1,3],
                                                   tornado_values[1,5],
                                                   tornado_values[1,7],
                                                   tornado_values[1,9],
                                                   tornado_values[1,11],
                                                   tornado_values[1,13],
                                                   tornado_values[1,15],
                                                   tornado_values[1,17],
                                                   tornado_values[1,19]))
    
    tornado_HCA_attributable <- data.frame("variable" = c("Discount Rate",
                                                          "Population Growth Rate",
                                                          "Willingness-to-Pay Threshold",
                                                          "Productivity Growth Rate",
                                                          "Labour Force Participation Rate",
                                                          "Incidence of Bacterial Infection",
                                                          "Rate of AMR",
                                                          "AMR Growth Rate",
                                                          "Livestock Industrialisation"),
                                           
                                           "min" = c(tornado_values[2,2],
                                                     tornado_values[2,4],
                                                     tornado_values[2,6],
                                                     tornado_values[2,8],
                                                     tornado_values[2,10],
                                                     tornado_values[2,12],
                                                     tornado_values[2,14],
                                                     tornado_values[2,16],
                                                     tornado_values[2,18]),
                                           
                                           "max" = c(tornado_values[2,3],
                                                     tornado_values[2,5],
                                                     tornado_values[2,7],
                                                     tornado_values[2,9],
                                                     tornado_values[2,11],
                                                     tornado_values[2,13],
                                                     tornado_values[2,15],
                                                     tornado_values[2,17],
                                                     tornado_values[2,19]))
    
    tornado_FCA_associated <- data.frame("variable" = c("Discount Rate",
                                                        "Population Growth Rate",
                                                        "Willingness-to-Pay Threshold",
                                                        "Productivity Growth Rate",
                                                        "Labour Force Participation Rate",
                                                        "Incidence of Bacterial Infection",
                                                        "Rate of AMR",
                                                        "AMR Growth Rate",
                                                        "Livestock Industrialisation"),
                                         
                                         "min" = c(tornado_values[3,2],
                                                   tornado_values[3,4],
                                                   tornado_values[3,6],
                                                   tornado_values[3,8],
                                                   tornado_values[3,10],
                                                   tornado_values[3,12],
                                                   tornado_values[3,14],
                                                   tornado_values[3,16],
                                                   tornado_values[3,18]),
                                         
                                         "max" = c(tornado_values[3,3],
                                                   tornado_values[3,5],
                                                   tornado_values[3,7],
                                                   tornado_values[3,9],
                                                   tornado_values[3,11],
                                                   tornado_values[3,13],
                                                   tornado_values[3,15],
                                                   tornado_values[3,17],
                                                   tornado_values[3,19]))
    
    tornado_FCA_attributable <- data.frame("variable" = c("Discount Rate",
                                                          "Population Growth Rate",
                                                          "Willingness-to-Pay Threshold",
                                                          "Productivity Growth Rate",
                                                          "Labour Force Participation Rate",
                                                          "Incidence of Bacterial Infection",
                                                          "Rate of AMR",
                                                          "AMR Growth Rate",
                                                          "Livestock Industrialisation"),
                                           
                                           "min" = c(tornado_values[4,2],
                                                     tornado_values[4,4],
                                                     tornado_values[4,6],
                                                     tornado_values[4,8],
                                                     tornado_values[4,10],
                                                     tornado_values[4,12],
                                                     tornado_values[4,14],
                                                     tornado_values[4,16],
                                                     tornado_values[4,18]),
                                           
                                           "max" = c(tornado_values[4,3],
                                                     tornado_values[4,5],
                                                     tornado_values[4,7],
                                                     tornado_values[4,9],
                                                     tornado_values[4,11],
                                                     tornado_values[4,13],
                                                     tornado_values[4,15],
                                                     tornado_values[4,17],
                                                     tornado_values[4,19]))
    
    #create the plots, display them together, and save
    
    tornado_HCA_associated$burden <- "associated"
    tornado_HCA_associated$productivity <- "HCA"
    
    tornado_HCA_attributable$burden <- "attributable"
    tornado_HCA_attributable$productivity <- "HCA"
    
    tornado_FCA_associated$burden <- "associated"
    tornado_FCA_associated$productivity <- "FCA"
    
    tornado_FCA_attributable$burden <- "attributable"
    tornado_FCA_attributable$productivity <- "FCA"
    
    width = 1
    base.value = 0
    
    #split into four dataframes to (hopefully) avoid duplication error
    
    order.parameters.HCA.associated <- tornado_df %>% 
      filter(productivity == "HCA", burden == "associated") %>%
      arrange(difference) %>%
      mutate(Parameter = factor(x = variable, levels = variable)) %>%
      select(Parameter) %>% 
      unlist() %>% 
      levels()
    
    order.parameters.HCA.attributable <- tornado_df %>% 
      filter(productivity == "HCA", burden == "attributable") %>%
      arrange(difference) %>%
      mutate(Parameter = factor(x = variable, levels = variable)) %>%
      select(Parameter) %>% 
      unlist() %>% 
      levels()
    
    order.parameters.FCA.associated <- tornado_df %>% 
      filter(productivity == "FCA", burden == "associated") %>%
      arrange(difference) %>%
      mutate(Parameter = factor(x = variable, levels = variable)) %>%
      select(Parameter) %>% 
      unlist() %>% 
      levels()
    
    order.parameters.FCA.attributable <- tornado_df %>% 
      filter(productivity == "FCA", burden == "attributable") %>%
      arrange(difference) %>%
      mutate(Parameter = factor(x = variable, levels = variable)) %>%
      select(Parameter) %>% 
      unlist() %>% 
      levels()
    
    tornado_df_HCA_associated <- tornado_HCA_associated %>%
      mutate(difference = abs(max - min))
    
    tornado_df_HCA_associated_2 <- tornado_df_HCA_associated %>%
      mutate(ordering = factor(variable, levels = order.parameters.HCA.associated)) %>%
      arrange(ordering) %>%
      mutate(para = seq(1,n())) %>%
      pivot_longer(cols = min:max) %>%
      mutate(ymin = pmin(value, base.value),
             ymax = pmax(value, base.value),
             xmin = as.numeric(para) - width/2,
             xmax = as.numeric(para) + width/2)
    
    tornado_df_HCA_associated_2$ymax_mils <- tornado_df_HCA_associated_2$ymax / 1e+06
    tornado_df_HCA_associated_2$ymin_mils <- tornado_df_HCA_associated_2$ymin / 1e+06
    
    g1 <- ggplot() +
      geom_rect(data <- tornado_df_HCA_associated_2,
                mapping = aes(ymax = ymax_mils, ymin = ymin_mils, xmax = xmax, xmin = xmin, fill = variable)) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title.y = element_blank(), legend.position = 'bottom',
            legend.title = element_blank()) +
      geom_hline(yintercept = base.value) +
      facet_grid(scales = "free") +
      scale_x_continuous(breaks = seq(1:9), labels = order.parameters.HCA.associated) +
      coord_flip() +
      scale_fill_discrete(limits = rev(order.parameters.HCA.associated)) +
      scale_y_continuous(limits = c(-1e+06, 1e+06), "Millions 2019 $USD over 20 Years") +
      ggtitle("Effect on Associated Burden (Human Captial Approach)")
    
    tornado_df_HCA_attributable <- tornado_HCA_attributable %>%
      mutate(difference = abs(max - min))
    
    tornado_df_HCA_attributable_2 <- tornado_df_HCA_attributable %>%
      mutate(ordering = factor(variable, levels = order.parameters.HCA.attributable)) %>%
      arrange(ordering) %>%
      mutate(para = seq(1,n())) %>%
      pivot_longer(cols = min:max) %>%
      mutate(ymin = pmin(value, base.value),
             ymax = pmax(value, base.value),
             xmin = as.numeric(para) - width/2,
             xmax = as.numeric(para) + width/2)
    
    tornado_df_HCA_attributable_2$ymax_mils <- tornado_df_HCA_attributable_2$ymax / 1e+06
    tornado_df_HCA_attributable_2$ymin_mils <- tornado_df_HCA_attributable_2$ymin / 1e+06
    
    g2 <- ggplot() +
      geom_rect(data <- tornado_df_HCA_attributable_2,
                mapping = aes(ymax = ymax_mils, ymin = ymin_mils, xmax = xmax, xmin = xmin, fill = variable)) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title.y = element_blank(), legend.position = 'bottom',
            legend.title = element_blank()) +
      geom_hline(yintercept = base.value) +
      facet_grid(scales = "free") +
      scale_x_continuous(breaks = seq(1:9), labels = order.parameters.HCA.attributable) +
      coord_flip() +
      scale_fill_discrete(limits = rev(order.parameters.HCA.attributable)) +
      scale_y_continuous(limits = c(-1e+06, 1e+06), "Millions 2019 $USD over 20 Years") +
      ggtitle("Effect on Attributable Burden (Human Captial Approach)")
    
    
    tornado_df_FCA_associated <- tornado_FCA_associated %>%
      mutate(difference = abs(max - min))
    
    tornado_df_FCA_associated_2 <- tornado_df_FCA_associated %>%
      mutate(ordering = factor(variable, levels = order.parameters.FCA.associated)) %>%
      arrange(ordering) %>%
      mutate(para = seq(1,n())) %>%
      pivot_longer(cols = min:max) %>%
      mutate(ymin = pmin(value, base.value),
             ymax = pmax(value, base.value),
             xmin = as.numeric(para) - width/2,
             xmax = as.numeric(para) + width/2)
    
    tornado_df_FCA_associated_2$ymax_mils <- tornado_df_FCA_associated_2$ymax / 1e+06
    tornado_df_FCA_associated_2$ymin_mils <- tornado_df_FCA_associated_2$ymin / 1e+06
    
    
    g3 <- ggplot() +
      geom_rect(data <- tornado_df_FCA_associated_2,
                mapping = aes(ymax = ymax_mils, ymin = ymin_mils, xmax = xmax, xmin = xmin, fill = variable)) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title.y = element_blank(), legend.position = 'bottom',
            legend.title = element_blank()) +
      geom_hline(yintercept = base.value) +
      facet_grid(scales = "free") +
      scale_x_continuous(breaks = seq(1:9), labels = order.parameters.FCA.associated) +
      coord_flip() +
      scale_fill_discrete(limits = rev(order.parameters.FCA.associated)) +
      scale_y_continuous(limits = c(-1e+06, 1e+06), "Millions 2019 $USD over 20 Years") +
      ggtitle("Effect on Associated Burden (Friction Cost Approach)")
    
    
    tornado_df_FCA_attributable <- tornado_FCA_attributable %>%
      mutate(difference = abs(max - min))
    
    tornado_df_FCA_attributable_2 <- tornado_df_FCA_attributable %>%
      mutate(ordering = factor(variable, levels = order.parameters.FCA.attributable)) %>%
      arrange(ordering) %>%
      mutate(para = seq(1,n())) %>%
      pivot_longer(cols = min:max) %>%
      mutate(ymin = pmin(value, base.value),
             ymax = pmax(value, base.value),
             xmin = as.numeric(para) - width/2,
             xmax = as.numeric(para) + width/2)
    
    tornado_df_FCA_attributable_2$ymax_mils <- tornado_df_FCA_attributable_2$ymax / 1e+06
    tornado_df_FCA_attributable_2$ymin_mils <- tornado_df_FCA_attributable_2$ymin / 1e+06
    
    g4 <- ggplot() +
      geom_rect(data <- tornado_df_FCA_attributable_2,
                mapping = aes(ymax = ymax_mils, ymin = ymin_mils, xmax = xmax, xmin = xmin, fill = variable)) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title.y = element_blank(), legend.position = 'bottom',
            legend.title = element_blank()) +
      geom_hline(yintercept = base.value) +
      facet_grid(scales = "free") +
      scale_x_continuous(breaks = seq(1:9), labels = order.parameters.FCA.attributable) +
      coord_flip() +
      scale_fill_discrete(limits = rev(order.parameters.FCA.attributable)) +
      scale_y_continuous(limits = c(-1e+06, 1e+06), "Millions 2019 $USD over 20 Years") +
      ggtitle("Effect on Attributable Burden (Friction Cost Approach)")
    
    
    #create model output
    tornadoplot <- ggarrange(g1, g2, g3, g4,
                             ncol = 2,
                             nrow = 2)
    
    return(tornadoplot)
  }
  

# render outputs for the UI -----------------------------------------------

  output$introtext <- renderText(introtext)
  output$inputslink <- renderText(inputslink)
  outputs$chart <- renderPlot(plotfunction(inputs))
  outputs$tornado <- renderPlot(tornadofunction(inputs))
  
  
}      
      
shinyApp(ui, server)      
      
      

    
      
      

      
      
      
