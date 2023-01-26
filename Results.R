
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


# Loading main model function ---------------------------------------------

source("Functions.R")


# Test run ----------------------------------------------------------------

Model(inputs, "HCA")



# Display and store outputs under HCA and FCA -----------------------------

#Table
main_outputs <- matrix(rep(0), ncol = 5, nrow = 11)
colnames(main_outputs) = c("Output", "Associated Burden (HCA)", "Associated Burden (FCA)",
                           "Attributable Burden (HCA)", "Attributable Burden (FCA)")

main_outputs[,1] <- Model(inputs, "HCA")[,1]
main_outputs[,2] <- as.numeric(Model(inputs, "HCA")[,2])
main_outputs[,3] <- as.numeric(Model(inputs, "FCA")[,2])
main_outputs[,4] <- as.numeric(Model(inputs, "HCA")[,3])
main_outputs[,5] <- as.numeric(Model(inputs, "FCA")[,3])

write.xlsx(main_outputs, "Outputs/Main Outputs.xlsx")

#Graphical representation

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

main_outputs_plot_HCA

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

main_outputs_plot_FCA

##Together
main_outputs_plot_HCA / main_outputs_plot_FCA
ggsave("Outputs/Main Outputs.jpg", width = 20, height = 15)


# Tornado Plots -----------------------------------------------------------

# Parameters to include in the plot: discount rate, population growth,
# wtp threshold, productivity growth rate, LFPR, incidence of infection, AMR incidence,
#AMR growth rate, portion of animals in industrial farms
# c(3, 5:7, 10, 27:28, 48)

#get the base case scenarios
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

#display and save
(g1 + g2) /
  (g3 + g4)
ggsave("Outputs/Tornado Plots.jpeg",width = 30, height = 20)



# Monte Carlo Simulation --------------------------------------------------

number_runs <- 100

MCmatrix <- matrix(rep(0), nrow = number_runs, ncol = 4)
colnames(MCmatrix) <- c("Associated Burden (HCA)", "Attributable Burden(HCA)",
                        "Associated Burdedn (FCA)", "Attributable Burden FCA)")

inputsPSA <- inputs

set.seed(80085)

for(i in 1:number_runs){
  
  print(i / number_runs)
  
  inputsPSA <- inputs
  
  inputsPSA[parameter == "well_sick", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "well_sick", "country param 1"]), 
                                                        as.numeric(inputsPSA[parameter == "well_sick", "country param 2"]))
  
  inputsPSA[parameter == "portion_res", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "portion_res", "country param 1"]), 
                                                          as.numeric(inputsPSA[parameter == "portion_res", "country param 2"]))
  
  inputsPSA[parameter == "mort_res", "Value"] <- 1.62 * rbeta(1,as.numeric(inputsPSA[parameter == "mort_res", "country param 1"]), 
                                                              as.numeric(inputsPSA[parameter == "mort_res", "country param 2"]))
  
  inputsPSA[parameter == "mort_sus", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "mort_sus", "country param 1"]), 
                                                       as.numeric(inputsPSA[parameter == "mort_sus", "country param 2"]))
  
  inputsPSA[parameter == "los_sus", "Value"] <- (1/365.25) * rlnorm(1, as.numeric(inputsPSA[parameter == "los_sus", "country param 1"]), 
                                                                    as.numeric(inputsPSA[parameter == "los_sus", "country param 2"]))
  
  inputsPSA[parameter == "los_res", "Value"] <- (1.27/365.25) * rlnorm(1, as.numeric(inputsPSA[parameter == "los_res", "country param 1"]), 
                                                                       as.numeric(inputsPSA[parameter == "los_res", "country param 2"]))
  
  inputsPSA[parameter == "bed_day_cost", "Value"] <- rlnorm(1, as.numeric(inputsPSA[parameter == "bed_day_cost", "country param 1"]), 
                                                            as.numeric(inputsPSA[parameter == "bed_day_cost", "country param 2"]))
  
  inputsPSA[parameter == "qol_sick", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "qol_sick", "country param 1"]), 
                                                       as.numeric(inputsPSA[parameter == "qol_sick", "country param 2"]))
  
  inputsPSA[parameter == "qol_seq", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "qol_seq", "country param 1"]), 
                                                      as.numeric(inputsPSA[parameter == "qol_seq", "country param 2"]))
  
  inputsPSA[parameter == "amr_grow", "Value"] <- 0.01 * rgamma(1, as.numeric(inputsPSA[parameter == "amr_grow", "country param 1"]), 
                                                               scale = as.numeric(inputsPSA[parameter == "amr_grow", "country param 2"]))
  
  inputsPSA[parameter == "n_pigs", "Value"] <- rlnorm(1, as.numeric(inputsPSA[parameter == "n_pigs", "country param 1"]), 
                                                      as.numeric(inputsPSA[parameter == "n_pigs", "country param 2"]))
  
  inputsPSA[parameter == "n_chickens", "Value"] <- rlnorm(1, as.numeric(inputsPSA[parameter == "n_chickens", "country param 1"]), 
                                                          as.numeric(inputsPSA[parameter == "n_chickens", "country param 2"]))
  
  inputsPSA[parameter == "n_chickens_farm_ind", "Value"] <- rlnorm(1, as.numeric(inputsPSA[parameter == "n_chickens_farm_ind", "country param 1"]), 
                                                                   as.numeric(inputsPSA[parameter == "n_chickens_farm_ind", "country param 2"]))
  
  inputsPSA[parameter == "n_chickens_farm_small", "Value"] <- rlnorm(1, as.numeric(inputsPSA[parameter == "n_chickens_farm_small", "country param 1"]), 
                                                                     as.numeric(inputsPSA[parameter == "n_chickens_farm_small", "country param 2"]))
  
  inputsPSA[parameter == "n_pigs_farm_ind", "Value"] <- rlnorm(1, as.numeric(inputsPSA[parameter == "n_pigs_farm_ind", "country param 1"]), 
                                                               as.numeric(inputsPSA[parameter == "n_pigs_farm_ind", "country param 2"]))
  
  inputsPSA[parameter == "n_pigs_farm_small", "Value"] <- rlnorm(1, as.numeric(inputsPSA[parameter == "n_pigs_farm_small", "country param 1"]), 
                                                                 as.numeric(inputsPSA[parameter == "n_pigs_farm_small", "country param 2"]))
  
  inputsPSA[parameter == "portion_animals_ind", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "portion_animals_ind", "country param 1"]), 
                                                                  as.numeric(inputsPSA[parameter == "portion_animals_ind", "country param 2"]))
  
  inputsPSA[parameter == "pig_price", "Value"] <- rlnorm(1, as.numeric(inputsPSA[parameter == "pig_price", "country param 1"]), 
                                                         as.numeric(inputsPSA[parameter == "pig_price", "country param 2"]))
  
  inputsPSA[parameter == "chicken_price", "Value"] <- rlnorm(1, as.numeric(inputsPSA[parameter == "chicken_price", "country param 1"]), 
                                                             as.numeric(inputsPSA[parameter == "chicken_price", "country param 2"]))
  
  inputsPSA[parameter == "c_mort_ind", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "c_mort_ind", "country param 1"]), 
                                                         as.numeric(inputsPSA[parameter == "c_mort_ind", "country param 2"]))
  
  inputsPSA[parameter == "c_mort_small", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "c_mort_small", "country param 1"]), 
                                                           as.numeric(inputsPSA[parameter == "c_mort_small", "country param 2"]))
  
  inputsPSA[parameter == "p_mort_ind", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "p_mort_ind", "country param 1"]), 
                                                         as.numeric(inputsPSA[parameter == "p_mort_ind", "country param 2"]))
  
  inputsPSA[parameter == "p_mort_small", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "p_mort_small", "country param 1"]), 
                                                           as.numeric(inputsPSA[parameter == "p_mort_small", "country param 2"]))
  
  inputsPSA[parameter == "seq_sus", "Value"] <- 1.54158796 * rbeta(1, as.numeric(inputsPSA[parameter == "seq_sus", "country param 1"]), 
                                                                   as.numeric(inputsPSA[parameter == "seq_sus", "country param 2"]))
  
  inputsPSA[parameter == "seq_res", "Value"] <- 2.4973725 * rbeta(1, as.numeric(inputsPSA[parameter == "seq_res", "country param 1"]), 
                                                                  as.numeric(inputsPSA[parameter == "seq_res", "country param 2"]))
  
  MCmatrix[i, 1] <- as.numeric(Model(inputsPSA, "HCA")[1,2])
  MCmatrix[i, 2] <- as.numeric(Model(inputsPSA, "HCA")[1,3])
  MCmatrix[i, 3] <- as.numeric(Model(inputsPSA, "FCA")[1,2])
  MCmatrix[i, 4] <- as.numeric(Model(inputsPSA, "FCA")[1,3])
  
}

#save outputs in a spreadsheet
write.xlsx(MCmatrix, "Outputs/Monte Carlo results.xlsx")

#to load existing MC outputs and save redoing the simulation:
# MCmatrix <- read.xlsx("Outputs/Monte Carlo results.xlsx", 1)[,2:5]
# colnames(MCmatrix) <- c("Associated Burden (HCA)", "Attributable Burden(HCA)",
#                         "Associated Burdedn (FCA)", "Attributable Burden FCA)")


#summary stats
MonteCarlo_summary_stats <- matrix(rep(0), ncol = 4, nrow = 5)
colnames(MonteCarlo_summary_stats) = c("Associated Burden (HCA)", "Attributable Burden (HCA)",
                                     "Associated Burden (FCA)", "Attributable Burden FCA)")
rownames(MonteCarlo_summary_stats) <- c("Minimum", "25th Percentile", "Median",
                                        "75th Percentile", "Maximum")
for(i in 1:4){
  MonteCarlo_summary_stats[1,i] <- min(MCmatrix[,i])
  MonteCarlo_summary_stats[2,i] <- as.numeric(quantile(MCmatrix[,i], .25))
  MonteCarlo_summary_stats[3,i] <- median(MCmatrix[,i])
  MonteCarlo_summary_stats[4,i] <- as.numeric(quantile(MCmatrix[,i], .75))
  MonteCarlo_summary_stats[5,i] <- max(MCmatrix[,i])
}

write.xlsx(MonteCarlo_summary_stats, "Outputs/Monte Carlo summary stats.xlsx")

#box and whisker chart
theme_set(theme_bw())

#combine the four columns into one long column with the type of burden as a variable
df_boxplot <- matrix(rep(0), ncol = 2, nrow = 4 * nrow(MCmatrix))

colnames(df_boxplot) <- c("Burden", "Estimation_Method")
df_boxplot[1:number_runs,2] <- "Associated Burden (HCA)"
df_boxplot[(number_runs+1):(2*number_runs),2] <- "Attributable Burden (HCA)"
df_boxplot[(2*number_runs+1):(3*number_runs),2] <- "Associated Burden (FCA)"
df_boxplot[(3*number_runs+1):(4*number_runs),2] <- "Attributable Burden (FCA)"

df_boxplot[1:number_runs,1] <- MCmatrix[,1]
df_boxplot[(number_runs+1):(2*number_runs),1] <- MCmatrix[,2]
df_boxplot[(2*number_runs+1):(3*number_runs),1] <- MCmatrix[,3]
df_boxplot[(3*number_runs+1):(4*number_runs),1] <- MCmatrix[,4]

df_boxplot <- as.data.frame(df_boxplot)

#convert the vars to correct format
df_boxplot$Estimation_Method <- as.factor(df_boxplot$Estimation_Method)
df_boxplot$Burden <- as.numeric(df_boxplot$Burden)

#divide by 1 million
df_boxplot$Burden <- df_boxplot$Burden / 1e+09

boxplot <- ggplot(df_boxplot, aes(x = Estimation_Method, y = Burden)) +
  geom_boxplot() +
  ggtitle("Distribution of AMR Burden over 10,000 Simulations, 2019 $USD Billions") +
  xlab("Estimation Method") 
  
boxplot

ggsave("Outputs/Boxplot.jpg", width = 12, height = 8)
