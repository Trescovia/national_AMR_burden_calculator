
#'This file contains the code needed to generate results using the AHHME-B model

# Packages ----------------------------------------------------------------

library("data.table")
library("readxl")
#library("stargazer")
library("tidyverse")
library("tseries") 
library("forecast") 
library("dynlm") 
library("seastests")
library("forecast")
library("TSA")
#library("epiR")
library("extraDistr")
library("MonoInc")
#library("pksensi")
#library("sensitivity")
library("xlsx")
library("gridExtra")
library("ggplot2")
library("reshape2")
library("here")
library("devtools")
#library("multisensi")
#library("rsq")
library("forcats")
library("patchwork")
library("scales")
library("ggpubr")

# Loading inputs and setting scenarios ------------------------------------

main_dir <- here()
sub_dir <- "Outputs"
ifelse(!dir.exists(file.path(main_dir, sub_dir)), dir.create(file.path(main_dir, sub_dir)), F)
rm("main_dir", "sub_dir")

inputs <- read.csv(here("input spreadsheet.csv"))
inputs <- as.data.table(inputs)
colnames(inputs) <- c("parameter", "description", "Value", "country min", "country max",
                      "Distribution", "country param 1", "country param 2")
scenario_prod <- "HCA"


# Loading main model functions --------------------------------------------

source("Functions.R")

# Test run ----------------------------------------------------------------

Model(inputs, "HCA")

# Display and store outputs under HCA and FCA -----------------------------

#Table
main_outputs <- table_function(inputs)
write.xlsx(main_outputs, "Outputs/Main Outputs.xlsx")

#Graphical representation
plot_function(inputs)
ggsave("Outputs/Main Outputs.jpg", width = 20, height = 15)

# Tornado Plots -----------------------------------------------------------

tornado_function(inputs)
ggsave("Outputs/Tornado Plots.jpeg",width = 30, height = 20)

# Monte Carlo Simulation --------------------------------------------------

#create spreadsheet

number_runs <- 100

MCmatrix <- MC_function(inputs, 100)
write.xlsx(MCmatrix, "Outputs/Monte Carlo results.xlsx")


#create summary stats

MCmatrix <- read.xlsx("Outputs/Monte Carlo results.xlsx", 1)[,2:5]
colnames(MCmatrix) <- c("Associated Burden (HCA)", "Attributable Burden(HCA)",
                        "Associated Burdedn (FCA)", "Attributable Burden FCA)")

MonteCarlo_summary_stats <- MC_sumstats_function(MCmatrix)

write.xlsx(MonteCarlo_summary_stats, "Outputs/Monte Carlo summary stats.xlsx")


# create box and whisker --------------------------------------------------

MC_boxplot_function(MCmatrix)
ggsave("Outputs/Boxplot.jpg", width = 12, height = 8)

