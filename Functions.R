

# What this script does ---------------------------------------------------

#'This file contains the functions needed to run the AHHME-B model

# Loading inputs (only for playing around with model) ---------------------

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

main_dir <- here()
sub_dir <- "Outputs"
ifelse(!dir.exists(file.path(main_dir, sub_dir)), dir.create(file.path(main_dir, sub_dir)), F)
rm("main_dir", "sub_dir")

inputs <- read.csv(here("input spreadsheet.csv"))
inputs <- as.data.table(inputs)
colnames(inputs) <- c("parameter", "description", "Value", "country min", "country max",
                              "Distribution", "country param 1", "country param 2")
scenario_prod <- "HCA"

# General Model -----------------------------------------------------------

Model <- function(inputs, scenario_prod){
  
  inputs[ , Value := as.numeric(as.character(Value))]
  
  n.t <- inputs[parameter=="n.t",Value] + 1
  
  # Other Functions ---------------------------------------------------------
  
  f_expvalue <- function(x,y,z){
    ## x is the epi matrix
    ## y is the cost matrix
    ## z is the reward matrix
    results <- matrix(c(sum(x*y),sum(x*z)),1,2)
    return(results)
    
  }
  
  f_di <- function(x,y){
    # function to apply a discount rate
    # x is cost
    # y is discount rate 
    x2 <- x - (x*y)
    return(x2)
  }
  
  # Population growth -------------------------------------------------------
  
  population <- c(rep(0,n.t+1))
  population[1] <- inputs[parameter == "pop", Value]
  for(i in 2:length(population)){
    population[i] <- population[i-1]*(1+inputs[parameter=="pop_growth", Value])
  }
  
  popchange <- c(rep(0,n.t))
  for(i in 1:n.t){
    popchange[i] <- population[i+1]-population[i]
  }
  
  #create vector of population relative to initial (for food production growth)
  pop_relative <- rep(0, length(population))
  for(i in 1:length(pop_relative)){
    pop_relative[i] <- population[i] / population[1]
  }
  
  # Human Epi Model ---------------------------------------------------------
  
  #building parameter matrix
  
  state_names <- c("well", "res","sus","dead", "was_well", "seq") ## the compartments
  transition_names  <- c("birth","r","s","mort_r", "mort_s","mort_w", "rec_r","rec_s", "dead_aft", "r_seq", "s_seq")  ## the transition probabilities
  parameter_names <- c(state_names, transition_names)
  
  state_i <- c(inputs[parameter=="pop",Value], rep(0,length=length(state_names)-1))
  #initial state vector
  
  m_param <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names))
  colnames(m_param) <- parameter_names
  rownames(m_param) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  m_param[ , "r_seq"] <- rep(inputs[parameter=="seq_res", Value], n.t) 
  #chance of developing sequelae following resistant infection
  
  m_param[ , "s_seq"] <- rep(inputs[parameter=="seq_sus", Value], n.t) 
  #chance of developing sequelae following susceptible infection
  
  m_param[ , "r"] <- rep(inputs[parameter=="well_sick", Value]*inputs[parameter=="portion_res",Value], n.t)
  #chance of developing a resistant infection in a year
  
  m_param[ , "s"] <- rep(inputs[parameter=="well_sick", Value]*(1-inputs[parameter=="portion_res",Value]), n.t)
  #chance of developing a susceptible infection in a year
  
  m_param[ , "mort_r"] <- rep(inputs[parameter=="mort_res", Value], n.t)
  #fatality from resistant infection
  
  m_param[ , "mort_s"] <- rep(inputs[parameter=="mort_sus", Value], n.t)
  #fatality from susceptible infection
  
  m_param[ , "rec_r"] <- rep(max(0,(1-(m_param[1,"mort_r"]+m_param[1,"r_seq"]))), n.t)
  #chance of recovering from a resistant infection
  
  m_param[ , "rec_s"] <- rep(max(0,(1-(m_param[1,"mort_s"]+m_param[1,"s_seq"]))), n.t)
  #chance of recovering from a susceptible infection
  
  m_param[ , "mort_w"] <- rep(0, n.t) 
  #chance of dying without infection from 'well', set to zero because background mortality is included in net births
  
  m_param[ , "dead_aft"] <- rep(1, n.t) 
  #'the 'aft' compartment is simply used for counting deaths, and individuals go there
  #'after they can no longer transition between states and all impacts have been 
  #'accounted for (i.e. they die or are in the sequelae state for the remainder
  #'of their life)
  
  m_param[ , "birth"] <- popchange[1:n.t]
  #predicted net births in a given year
  
  
  m_param[1, 1:length(state_names)] <- state_i 
  #adding initial cycle 0 values
  
  #'have growing AMR prevalence
  #'note that we will later ensure that the total disease incidence does not change,
  #'only the portion of infections from resistant bacteria
  #'this may overestimate the impact of our intervention if improvements in living
  #'standards cause the overall incidence of disease to fall over time
  amr_growth <- inputs[parameter=="amr_grow", Value]
  
  #set the maximum portion of resistant infections
  max_res <- inputs[parameter == "max_r", Value]
  
  for (i in 2:(n.t)){
    m_param[i, "r"] <- m_param[i-1, "r"]*amr_growth
  }
  
  for(i in 1:n.t){
    if(m_param[i, "r"] > max_res *(m_param[1,"r"]+m_param[1,"s"])){ 
      m_param[i, "r"] <- max_res *(m_param[1,"r"]+m_param[1,"s"]) 
    }
    m_param[i, "s"] <- m_param[1,"r"]+m_param[1,"s"] - m_param[i, "r"] 
  }
  #made it so that, while the incidence of resistant infections increases,
  #the total number of infections doesn't increase. However, if we wanted to 
  #have a changing number of infections over time, we could do that, but would
  #need to replace "m_param[1,"r"]+m_param[1,"s"]" with a parameter called 
  #something like "chance_sick" which can change each period, and have that 
  #grow before allowing the AMR incidence to grow
  
  #make sure that the transition probabilities don't exceed 1
  for(i in 1:n.t){
    m_param[i, "mort_r"] <- m_param[i, "mort_r"] / (m_param[i, "mort_r"] + m_param[i, "rec_r"] + m_param[i, "r_seq"])
    m_param[i, "rec_r"] <- m_param[i, "rec_r"] / (m_param[i, "mort_r"] + m_param[i, "rec_r"] + m_param[i, "r_seq"])
    m_param[i, "r_seq"] <- m_param[i, "r_seq"] / (m_param[i, "mort_r"] + m_param[i, "rec_r"] + m_param[i, "r_seq"])
    
    m_param[i, "mort_s"] <- m_param[i, "mort_s"] / (m_param[i, "mort_s"] + m_param[i, "rec_s"] + m_param[i, "s_seq"])
    m_param[i, "rec_s"] <- m_param[i, "rec_s"] / (m_param[i, "mort_s"] + m_param[i, "rec_s"] + m_param[i, "s_seq"])
    m_param[i, "s_seq"] <- m_param[i, "s_seq"] / (m_param[i, "mort_s"] + m_param[i, "rec_s"] + m_param[i, "s_seq"])
  }
  
  m_param_blank <- m_param
  
  ##set the initial state
  
  #born
  m_param[1, 1:length(state_names)] <- state_i
  m_param[1, "well"] <- m_param[1, "well"] + m_param[1, "birth"]
  
  #transition out of well
  m_param[1, "was_well"] <- m_param[1, "well"]
  m_param[1, "well"] <- m_param[1, "was_well"] * (1 - m_param[1, "r"] - m_param[1, "s"] - m_param[1, "mort_w"])
  
  m_param[1, "res"] <- m_param[1, "was_well"] * m_param[1, "r"]
  m_param[1, "sus"] <- m_param[1, "was_well"] * m_param[1, "s"]
  m_param[1, "dead"] <- m_param[1, "was_well"] * m_param[1, "mort_w"]
  
  #transition out of sick
  m_param[1, "well"] <- m_param[1, "well"] + 
    (m_param[1, "rec_r"] * m_param[1, "res"]) + 
    (m_param[1, "rec_s"] * m_param[1, "sus"])
  
  m_param[1, "dead"] <- m_param[1, "dead"] +
    (m_param[1, "mort_s"] * m_param[1, "sus"]) + 
    (m_param[1, "mort_r"] * m_param[1, "res"])
  
  m_param[1, "seq"] <- m_param[1, "seq"] +
    (m_param[1, "s_seq"] * m_param[1, "sus"]) + 
    (m_param[1, "r_seq"] * m_param[1, "res"])
  
  ##difference equation
  
  f_human_epi <- function(m_param, n.t){
    
    n.t.val <- n.t
    
    for(i in 2:n.t.val){
      
      #carry over from last period and be born
      m_param[i, "well"] <- m_param[i-1, "well"] + m_param[i, "birth"]
      
      m_param[i, "was_well"] <- m_param[i, "well"]
      
      #transition out of well
      m_param[i, "well"] <- m_param[i, "was_well"] * (1 - m_param[i, "r"] - m_param[i, "s"] - m_param[i, "mort_w"])
      
      m_param[i, "res"] <- m_param[i, "was_well"] * m_param[i, "r"]
      m_param[i, "sus"] <- m_param[i, "was_well"] * m_param[i, "s"]
      m_param[i, "dead"] <- m_param[i, "was_well"] * m_param[i, "mort_w"]
      
      #transition out of sick
      m_param[i, "well"] <- m_param[i, "well"] + 
        (m_param[i, "rec_r"] * m_param[i, "res"]) + 
        (m_param[i, "rec_s"] * m_param[i, "sus"])
      
      m_param[i, "dead"] <- m_param[i, "dead"] +
        (m_param[i, "mort_s"] * m_param[i, "sus"]) + 
        (m_param[i, "mort_r"] * m_param[i, "res"])
      
      m_param[i, "seq"] <- m_param[i, "seq"] +
        (m_param[i, "s_seq"] * m_param[i, "sus"]) + 
        (m_param[i, "r_seq"] * m_param[i, "res"])
      
    }
    
    return(m_param)
  }
  
  m_param <- f_human_epi(m_param,n.t) 
  
  # Healthcare Costs --------------------------------------------------------
  
  dr <- inputs[parameter == "dr", Value]
  
  m_cost <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names))
  colnames(m_cost) <- parameter_names
  rownames(m_cost) <- paste("cycle", 0:(n.t-1), sep  =  "")  
  
  #cost of hospital stay for res and sus infection
  
  c_r <- inputs[parameter == "los_sus", Value] * inputs[parameter == "bed_day_cost", Value] * 365.25
  c_s <- inputs[parameter == "los_res", Value] * inputs[parameter == "bed_day_cost", Value] * 365.25
  
  cost_i <- c(0,c_r,c_s,0,0,0) #only infections incur a healthcare cost here
  
  ## start at cycle 1 so you do not multiply initial state vector 
  m_cost[2, 1:length(state_names)] <- cost_i
  
  for (j in 1:length(state_names)) {
    for (i in 3:(n.t)){
      m_cost[i,j] <- f_di(m_cost[i-1,j],dr)
    }  
  }
  
  
  # Healthcare Rewards ------------------------------------------------------
  
  m_rwd <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names))
  colnames(m_rwd) <- parameter_names
  rownames(m_rwd) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #calculate the present value of expected remaining life years for a) a perfectly
  #healthy person and b) someone with sequelae. The negative difference is the 'reward' 
  #of being in the sequelae state 
  
  remaining_ly <- inputs[parameter=="remaining_ly", Value]
  
  #pv_fut_life <- c(rep(0,n.t-1))
  pv_fut_life <- c(rep(0, remaining_ly))
  
  #for (i in 1:n.t-1){
  for(i in 1:length(pv_fut_life)){
    pv_fut_life[i] <- inputs[parameter=="background_qol", Value] * (1-dr)^(i-1)
  }
  pv_life <- sum(pv_fut_life)
  
  #pv_fut_life_seq <- c(rep(0,n.t-1))
  pv_fut_life_seq <- c(rep(0,remaining_ly))
  
  #for (i in 1:n.t-1){
  for(i in 1:length(pv_fut_life_seq)){
    pv_fut_life_seq[i] <- inputs[parameter=="qol_seq", Value] * (1-dr)^(i-1)
  }
  pv_life_seq <- sum(pv_fut_life_seq)
  
  #the 'reward' for death is not sweet release, but in fact the negative of the present value of remaining life 
  #years, and the 'reward' for being infected is the loss of welfare while hospitalised
  #therefore in a scenario with more deaths, infections and sequelae will have a negative
  #'reward' of larger absolute value
  
  r_s <- inputs[parameter == "los_sus", Value] *
    (inputs[parameter == "qol_sick", Value] - inputs[parameter == "background_qol", Value])
  
  r_r <- inputs[parameter == "los_res", Value] * 
    (inputs[parameter == "qol_sick", Value] - inputs[parameter == "background_qol", Value])
  
  r_d <- -1 * pv_life #discounted QoL loss from death
  r_seq <- pv_life_seq - pv_life # fixed this because we were previously assigning a benefit to sequelae (the subtraction was the wrong way around)
  
  rwd_i <- c(0,r_r,r_s,r_d,0,r_seq) 
  
  
  ## start at cycle 1 so you do not multiply initial state vector 
  m_rwd[2, 1:length(state_names)] <- rwd_i
  
  ### accounting for discounting
  for (j in 1:length(state_names)) {
    for (i in 3:(n.t)){
      m_rwd[i,j] <- f_di(m_rwd[i-1,j],dr)
    }  
  }
  
  
  # Productivity Costs -------------------------------------------------------
  
  #for HCA and FCA, we only care about the losses in productivity, so we set the 
  #reward for 'well' to be zero. Going into the 'dead' state incurs a productivity
  #loss equal to the discounted value of forgone future earnings, going into the
  #'res' or 'sus states incurs a loss equal to the earnings that would have
  #been made during the time in hospital. After 1 period, all people in 'dead' go
  #to the 'aft' compartment, which has a reward of zero. 
  #Using the HCA, the forgone future earnings are those of expected remaining economically
  #active years. For FCA, it is the forgone earnings during the time that it takes
  #to find a replacement worker
  
  #all productivity rewards are set to zero, and we only see a difference between
  #the reward matrices of different scenarios
  m_cost_prod <- matrix(rep(0),nrow = n.t, ncol = length(parameter_names))
  colnames(m_cost_prod) <- parameter_names
  rownames(m_cost_prod) <- paste("cycle", 0:(n.t-1), sep = "")
  
  cost_i_prod <- rep(0,(length(state_names)))
  
  ## start at cycle 1 so you do not multiply initial state vector
  m_cost_prod[2,1:length(state_names)] <- cost_i_prod
  
  for (j in 1:length(state_names)) {
    for (i in 3:(n.t)){
      m_cost_prod[i,j] <- f_di(m_cost_prod[i-1,j],dr)
    }
  }
  
  
  # Productivity Rewards ----------------------------------------------------
  
  #the 'reward' for being infected is equal to the negative forgone productivity
  #while hospitalised
  
  m_rwd_prod <- matrix(rep(0), nrow = n.t, ncol = length(parameter_names))
  colnames(m_rwd_prod) <- parameter_names
  rownames(m_rwd_prod) <- paste("cycle", 0:(n.t-1), sep = "")
  
  r_r_prod <- -1 * inputs[parameter == "los_res", Value] * 
    ((inputs[parameter == "prod_pc", Value]*inputs[parameter == "lfpr", Value])* 
       inputs[parameter == "unpaid_prod_pc", Value])
  
  r_s_prod <- -1 * inputs[parameter == "los_sus", Value] *
    ((inputs[parameter == "prod_pc", Value]*inputs[parameter == "lfpr", Value])*
       inputs[parameter == "unpaid_prod_pc", Value])
  
  r_w_prod <- 0
  
  r_aft_prod <- 0
  
  r_seq_prod <- 0 #importantly assumes that people with sequelae are equally productive
  
  #the reward for dead is the present discounted value of future work (either for the
  #remainder of economically active life, or for the 6 months needed to find a replacement)
  
  remaining_work_years <- inputs[parameter == "remaining_work_years", Value]
  
  yearly_prod <- inputs[parameter == "prod_pc", Value]*inputs[parameter == "lfpr", Value]*
    inputs[parameter == "unpaid_prod_pc", Value]
  
  pv_fut_prod <- c(rep(0,remaining_work_years))
  
  dr_pgrowth <- dr - inputs[parameter=="prod_growth", Value] 
  ##discount rate net of productivity growth (in theory it can be negative)
  
  for (i in 1:remaining_work_years){ 
    pv_fut_prod[i] <- yearly_prod * (1-dr_pgrowth)^(i-1)
  }
  
  pv_life_prod <- sum(pv_fut_prod)
  
  #'productivity effect depends on whether we use the human capital approach (HCA) 
  #'or friction cost approach (FCA) to estimating productivity losses (below)
  #' 
  #' if we use the HCA, then a death incurs a loss equal to all of the labour 
  #' productivity that the dying person would have contributed during the remainder
  #' of their life had they not died. If we use the FCA, then it is assumed that
  #' there exists a reserve army of labour from which new workers can be drawn in 
  #' the event of somebody's death, meaning that we lose the productivity that
  #' the dying person would have contributed during the time that it takes to find
  #' a replacement worker (this is typically assumed to be six months)
  
  
  if(scenario_prod == "HCA"){
    r_d_prod <- -1 * pv_life_prod
  } else if(scenario_prod == "FCA"){
    r_d_prod <- -0.5 * yearly_prod
  } else{
    paste("ERROR: PLEASE CHOOSE AN APPROACH TO ESTIMATING PRODUCTIVITY OUTCOMES")
  }
  
  rwd_i_prod <- c(r_w_prod, r_r_prod, r_s_prod, r_d_prod, r_aft_prod, r_seq_prod)
  
  ## start at cycle 1 so you do not multiply initial state vector
  m_rwd_prod[2,1:length(state_names)] <- rwd_i_prod 
  
  ### discount, but also account for labour productivity growth
  
  for (j in 1:length(state_names)) {
    for (i in 3:(n.t)){
      m_rwd_prod[i,j] <- f_di(m_rwd_prod[i-1,j],dr_pgrowth)
    }
  }
  
  # Human epi (associated burden counterfactual) ----------------------------
  
  #'here, we propose a counterfactual where all resistant infections are counted
  #'as healthy
  
  #building parameter matrix
  
  # state_names <- c("well", "res","sus","dead", "was_well", "seq") ## the compartments
  # transition_names  <- c("birth","r","s","mort_r", "mort_s","mort_w", "rec_r","rec_s", "dead_aft", "r_seq", "s_seq")  ## the transition probabilities
  # parameter_names <- c(state_names, transition_names)
  
  # state_i <- c(inputs[parameter=="pop",Value], rep(0,length=length(state_names)-1))
  # #initial state vector
  
  #' m_param_associated <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names))
  #' colnames(m_param_associated) <- parameter_names
  #' rownames(m_param_associated) <- paste("cycle", 0:(n.t-1), sep  =  "")
  #' 
  #' m_param_associated[ , "r_seq"] <- rep(inputs[parameter=="seq_res", Value], n.t) 
  #' #chance of developing sequelae following resistant infection
  #' 
  #' m_param_associated[ , "s_seq"] <- rep(inputs[parameter=="seq_sus", Value], n.t) 
  #' #chance of developing sequelae following susceptible infection
  #' 
  #' m_param_associated[ , "r"] <- rep(inputs[parameter=="well_sick", Value]*inputs[parameter=="portion_res",Value], n.t)
  #' #chance of developing a resistant infection in a year (we keep this here, but during the state transition we actually 
  #' #put them in 'well')
  #' 
  #' m_param_associated[ , "s"] <- rep(inputs[parameter=="well_sick", Value]*(1-inputs[parameter=="portion_res",Value]), n.t)
  #' #chance of developing a susceptible infection in a year
  #' 
  #' m_param_associated[ , "mort_r"] <- rep(inputs[parameter=="mort_res", Value], n.t)
  #' #fatality from resistant infection
  #' 
  #' m_param_associated[ , "mort_s"] <- rep(inputs[parameter=="mort_sus", Value], n.t)
  #' #fatality from susceptible infection
  #' 
  #' m_param_associated[ , "rec_r"] <- rep(max(0,(1-(m_param[1,"mort_r"]+m_param[1,"r_seq"]))), n.t)
  #' #chance of recovering from a resistant infection
  #' 
  #' m_param_associated[ , "rec_s"] <- rep(max(0,(1-(m_param[1,"mort_s"]+m_param[1,"s_seq"]))), n.t)
  #' #chance of recovering from a susceptible infection
  #' 
  #' m_param_associated[ , "mort_w"] <- rep(0, n.t) 
  #' #chance of dying without infection from 'well', set to zero because background mortality is included in net births
  #' 
  #' m_param_associated[ , "dead_aft"] <- rep(1, n.t) 
  #' #'the 'aft' compartment is simply used for counting deaths, and individuals go there
  #' #'after they can no longer transition between states and all impacts have been 
  #' #'accounted for (i.e. they die or are in the sequelae state for the remainder
  #' #'of their life)
  #' 
  #' m_param_associated[ , "birth"] <- popchange[1:n.t]
  #' #predicted net births in a given year
  #' 
  #' state_i
  #' 
  #' m_param_associated[1, 1:length(state_names)] <- state_i 
  #' #adding initial cycle 0 values
  #' 
  #' #' #'have growing AMR prevalence
  #' #' #'note that we will later ensure that the total disease incidence does not change,
  #' #' #'only the portion of infections from resistant bacteria
  #' #' #'this may overestimate the impact of our intervention if improvements in living
  #' #' #'standards cause the overall incidence of disease to fall over time
  #' #' amr_growth <- inputs[parameter=="amr_grow", Value]
  #' #' 
  #' #' #set the maximum portion of resistant infections
  #' #' max_res <- inputs[parameter == "max_r", Value]
  #' 
  #' for (i in 2:(n.t)){
  #'   m_param_associated[i, "r"] <- m_param_associated[i-1, "r"]*amr_growth
  #' }
  #' 
  #' for(i in 1:n.t){
  #'   if(m_param_associated[i, "r"] > max_res *(m_param_associated[1,"r"]+m_param_associated[1,"s"])){ 
  #'     m_param_associated[i, "r"] <- max_res *(m_param_associated[1,"r"]+m_param_associated[1,"s"]) 
  #'   }
  #'   m_param_associated[i, "s"] <- m_param_associated[1,"r"]+m_param_associated[1,"s"] - m_param_associated[i, "r"] 
  #' }
  #' #made it so that, while the incidence of resistant infections increases,
  #' #the total number of infections doesn't increase. However, if we wanted to 
  #' #have a changing number of infections over time, we could do that, but would
  #' #need to replace "m_param[1,"r"]+m_param[1,"s"]" with a parameter called 
  #' #something like "chance_sick" which can change each period, and have that 
  #' #grow before allowing the AMR incidence to grow
  #' 
  #' #make sure that the transition probabilities don't exceed 1
  #' for(i in 1:n.t){
  #'   m_param_associated[i, "mort_r"] <- m_param_associated[i, "mort_r"] / (m_param_associated[i, "mort_r"] + m_param_associated[i, "rec_r"] + m_param_associated[i, "r_seq"])
  #'   m_param_associated[i, "rec_r"] <- m_param_associated[i, "rec_r"] / (m_param_associated[i, "mort_r"] + m_param_associated[i, "rec_r"] + m_param_associated[i, "r_seq"])
  #'   m_param_associated[i, "r_seq"] <- m_param_associated[i, "r_seq"] / (m_param_associated[i, "mort_r"] + m_param_associated[i, "rec_r"] + m_param_associated[i, "r_seq"])
  #'   
  #'   m_param_associated[i, "mort_s"] <- m_param_associated[i, "mort_s"] / (m_param_associated[i, "mort_s"] + m_param_associated[i, "rec_s"] + m_param_associated[i, "s_seq"])
  #'   m_param_associated[i, "rec_s"] <- m_param_associated[i, "rec_s"] / (m_param_associated[i, "mort_s"] + m_param_associated[i, "rec_s"] + m_param_associated[i, "s_seq"])
  #'   m_param_associated[i, "s_seq"] <- m_param_associated[i, "s_seq"] / (m_param_associated[i, "mort_s"] + m_param_associated[i, "rec_s"] + m_param_associated[i, "s_seq"])
  #' }
  
  # m_param_blank <- m_param_associated
  #so we can use this as the basis for the attributable scenario and save ourselves 100 lines of code
  
  m_param_associated <- m_param_blank
  
  #make it so that the chance of res infection is zero
  for(i in 1:nrow(m_param_associated)){
    m_param_associated[i, "r"] <- 0
  }
  
  #set the initial state
  
  #born
  m_param_associated[1, 1:length(state_names)] <- state_i
  m_param_associated[1, "well"] <- m_param_associated[1, "well"] + m_param_associated[1, "birth"]
  
  #transition out of well
  m_param_associated[1, "was_well"] <- m_param_associated[1, "well"]
  m_param_associated[1, "well"] <- m_param_associated[1, "was_well"] * (1 - m_param_associated[1, "r"] - m_param_associated[1, "s"] - m_param_associated[1, "mort_w"])
  
  m_param_associated[1, "res"] <- m_param_associated[1, "was_well"] * m_param_associated[1, "r"]
  m_param_associated[1, "sus"] <- m_param_associated[1, "was_well"] * m_param_associated[1, "s"]
  m_param_associated[1, "dead"] <- m_param_associated[1, "was_well"] * m_param_associated[1, "mort_w"]
  
  #transition out of sick
  m_param_associated[1, "well"] <- m_param_associated[1, "well"] + 
    (m_param_associated[1, "rec_r"] * m_param_associated[1, "res"]) + 
    (m_param_associated[1, "rec_s"] * m_param_associated[1, "sus"])
  
  m_param_associated[1, "dead"] <- m_param_associated[1, "dead"] +
    (m_param_associated[1, "mort_s"] * m_param_associated[1, "sus"]) + 
    (m_param_associated[1, "mort_r"] * m_param_associated[1, "res"])
  
  m_param_associated[1, "seq"] <- m_param_associated[1, "seq"] +
    (m_param_associated[1, "s_seq"] * m_param_associated[1, "sus"]) + 
    (m_param_associated[1, "r_seq"] * m_param_associated[1, "res"])
  
  #apply the difference equation
  m_param_associated <- f_human_epi(m_param_associated,n.t)
  
  # ##set the initial state
  # 
  # #born
  # m_param_associated[1, 1:length(state_names)] <- state_i
  # m_param_associated[1, "well"] <- m_param_associated[1, "well"] + m_param_associated[1, "birth"]
  # 
  # #transition out of well
  # m_param_associated[1, "was_well"] <- m_param_associated[1, "well"]
  # #m_param_associated[1, "well"] <- m_param_associated[1, "was_well"] * (1 - m_param_associated[1, "r"] - m_param_associated[1, "s"] - m_param_associated[1, "mort_w"])
  # m_param_associated[1, "well"] <- m_param_associated[1, "was_well"] * (1 - m_param_associated[1, "s"] - m_param_associated[1, "mort_w"])
  # 
  # #m_param_associated[1, "res"] <- m_param_associated[1, "was_well"] * m_param_associated[1, "r"]
  # m_param_associated[1, "sus"] <- m_param_associated[1, "was_well"] * m_param_associated[1, "s"]
  # m_param_associated[1, "dead"] <- m_param_associated[1, "was_well"] * m_param_associated[1, "mort_w"]
  # 
  # #transition out of sick
  # m_param_associated[1, "well"] <- m_param_associated[1, "well"] + 
  #   (m_param_associated[1, "rec_r"] * m_param_associated[1, "res"]) + 
  #   (m_param_associated[1, "rec_s"] * m_param_associated[1, "sus"])
  # 
  # m_param_associated[1, "dead"] <- m_param_associated[1, "dead"] +
  #   (m_param_associated[1, "mort_s"] * m_param_associated[1, "sus"]) + 
  #   (m_param_associated[1, "mort_r"] * m_param_associated[1, "res"])
  # 
  # m_param_associated[1, "seq"] <- m_param_associated[1, "seq"] +
  #   (m_param_associated[1, "s_seq"] * m_param_associated[1, "sus"]) + 
  #   (m_param_associated[1, "r_seq"] * m_param_associated[1, "res"])
  # 
  # ##difference equation
  # 
  # f_human_epi_associated <- function(m_param_associated, n.t){
  #   
  #   n.t.val <- n.t
  #   
  #   for(i in 2:n.t.val){
  #     
  #     #carry over from last period and be born
  #     m_param_associated[i, "well"] <- m_param_associated[i-1, "well"] + m_param_associated[i, "birth"]
  #     
  #     m_param_associated[i, "was_well"] <- m_param_associated[i, "well"]
  #     
  #     #transition out of well
  #     #m_param[i, "well"] <- m_param[i, "was_well"] * (1 - m_param[i, "r"] - m_param[i, "s"] - m_param[i, "mort_w"])
  #     m_param_associated[i, "well"] <- m_param_associated[i, "was_well"] * (1 - m_param_associated[i, "s"] - m_param_associated[i, "mort_w"])
  #     
  #     #m_param[i, "res"] <- m_param[i, "was_well"] * m_param[i, "r"]
  #     m_param_associated[i, "sus"] <- m_param_associated[i, "was_well"] * m_param_associated[i, "s"]
  #     m_param_associated[i, "dead"] <- m_param_associated[i, "was_well"] * m_param_associated[i, "mort_w"]
  #     
  #     #transition out of sick
  #     m_param_associated[i, "well"] <- m_param_associated[i, "well"] + 
  #       (m_param_associated[i, "rec_r"] * m_param_associated[i, "res"]) + 
  #       (m_param_associated[i, "rec_s"] * m_param_associated[i, "sus"])
  #     
  #     m_param_associated[i, "dead"] <- m_param_associated[i, "dead"] +
  #       (m_param_associated[i, "mort_s"] * m_param_associated[i, "sus"]) + 
  #       (m_param_associated[i, "mort_r"] * m_param_associated[i, "res"])
  #     
  #     m_param_associated[i, "seq"] <- m_param_associated[i, "seq"] +
  #       (m_param_associated[i, "s_seq"] * m_param_associated[i, "sus"]) + 
  #       (m_param_associated[i, "r_seq"] * m_param_associated[i, "res"])
  #     
  #   }
  #   
  #   return(m_param_associated)
  # }
  
  # m_param_associated <- f_human_epi_associated(m_param_associated,n.t) 
  
  #' okay so it partially works as there are no res infections and the deaths
  #' change appropriately. The number of sus infections is higher in the associated 
  #' scenario, but this actually makes sense because there are slightly more healthy people
  #' who can transition into 'sus' (so basically because we assume that you can't have both).
  
  # Human epi (attributable burden counterfactual) --------------------------
  
  #'here, we propose a counterfactual where all resistant infections are counted
  #'as susceptible infections
  m_param_attributable <- m_param_blank
  
  #make it so that res infections go to sus
  for(i in 1:nrow(m_param_attributable)){
    m_param_attributable[i, "s"] <- m_param_attributable[i, "s"] + m_param_attributable[i, "r"]
    m_param_attributable[i, "r"] <- 0
  }
  
  #set the initial state
  
  #born
  m_param_attributable[1, 1:length(state_names)] <- state_i
  m_param_attributable[1, "well"] <- m_param_attributable[1, "well"] + m_param_attributable[1, "birth"]
  
  #transition out of well
  m_param_attributable[1, "was_well"] <- m_param_attributable[1, "well"]
  m_param_attributable[1, "well"] <- m_param_attributable[1, "was_well"] * (1 - m_param_attributable[1, "r"] - m_param_attributable[1, "s"] - m_param_attributable[1, "mort_w"])
  
  m_param_attributable[1, "res"] <- m_param_attributable[1, "was_well"] * m_param_attributable[1, "r"]
  m_param_attributable[1, "sus"] <- m_param_attributable[1, "was_well"] * m_param_attributable[1, "s"]
  m_param_attributable[1, "dead"] <- m_param_attributable[1, "was_well"] * m_param_attributable[1, "mort_w"]
  
  #transition out of sick
  m_param_attributable[1, "well"] <- m_param_attributable[1, "well"] + 
    (m_param_attributable[1, "rec_r"] * m_param_attributable[1, "res"]) + 
    (m_param_attributable[1, "rec_s"] * m_param_attributable[1, "sus"])
  
  m_param_attributable[1, "dead"] <- m_param_attributable[1, "dead"] +
    (m_param_attributable[1, "mort_s"] * m_param_attributable[1, "sus"]) + 
    (m_param_attributable[1, "mort_r"] * m_param_attributable[1, "res"])
  
  m_param_attributable[1, "seq"] <- m_param_attributable[1, "seq"] +
    (m_param_attributable[1, "s_seq"] * m_param_attributable[1, "sus"]) + 
    (m_param_attributable[1, "r_seq"] * m_param_attributable[1, "res"])
  
  #apply the difference equation
  m_param_attributable <- f_human_epi(m_param_attributable,n.t)
  
  # #born
  # m_param_attributable[1, 1:length(state_names)] <- state_i
  # m_param_attributable[1, "well"] <- m_param_attributable[1, "well"] + m_param_attributable[1, "birth"]
  # 
  # #transition out of well
  # m_param_attributable[1, "was_well"] <- m_param_attributable[1, "well"]
  # m_param_attributable[1, "well"] <- m_param_attributable[1, "was_well"] * (1 - m_param_attributable[1, "r"] - m_param_attributable[1, "s"] - m_param_attributable[1, "mort_w"])
  # 
  # #m_param_attributable[1, "res"] <- m_param_attributable[1, "was_well"] * m_param_attributable[1, "r"]
  # m_param_attributable[1, "sus"] <- m_param_attributable[1, "was_well"] * (m_param_attributable[1, "s"] + m_param_attributable[1, "r"])
  # m_param_attributable[1, "dead"] <- m_param_attributable[1, "was_well"] * m_param_attributable[1, "mort_w"]
  # 
  # #transition out of sick
  # m_param_attributable[1, "well"] <- m_param_attributable[1, "well"] + 
  #   (m_param_attributable[1, "rec_r"] * m_param_attributable[1, "res"]) + 
  #   (m_param_attributable[1, "rec_s"] * m_param_attributable[1, "sus"])
  # 
  # m_param_attributable[1, "dead"] <- m_param_attributable[1, "dead"] +
  #   (m_param_attributable[1, "mort_s"] * m_param_attributable[1, "sus"]) + 
  #   (m_param_attributable[1, "mort_r"] * m_param_attributable[1, "res"])
  # 
  # m_param_attributable[1, "seq"] <- m_param_attributable[1, "seq"] +
  #   (m_param_attributable[1, "s_seq"] * m_param_attributable[1, "sus"]) + 
  #   (m_param_attributable[1, "r_seq"] * m_param_attributable[1, "res"])
  # 
  # ##difference equation
  # 
  # f_human_epi_attributable <- function(m_param_attributable, n.t){
  #   
  #   n.t.val <- n.t
  #   
  #   for(i in 2:n.t.val){
  #     
  #     #carry over from last period and be born
  #     m_param_attributable[i, "well"] <- m_param_attributable[i-1, "well"] + m_param_attributable[i, "birth"]
  #     
  #     m_param_attributable[i, "was_well"] <- m_param_attributable[i, "well"]
  #     
  #     #transition out of well
  #     m_param_attributable[i, "well"] <- m_param_attributable[i, "was_well"] * (1 - m_param_attributable[i, "r"] - m_param_attributable[i, "s"] - m_param_attributable[i, "mort_w"])
  #     
  #     #m_param_attributable[i, "res"] <- m_param_attributable[i, "was_well"] * m_param_attributable[i, "r"]
  #     m_param_attributable[i, "sus"] <- m_param_attributable[i, "was_well"] * (m_param_attributable[i, "s"] + m_param_attributable[i, "r"])
  #     m_param_attributable[i, "dead"] <- m_param_attributable[i, "was_well"] * m_param_attributable[i, "mort_w"]
  #     
  #     #transition out of sick
  #     m_param_attributable[i, "well"] <- m_param_attributable[i, "well"] + 
  #       (m_param_attributable[i, "rec_r"] * m_param_attributable[i, "res"]) + 
  #       (m_param_attributable[i, "rec_s"] * m_param_attributable[i, "sus"])
  #     
  #     m_param_attributable[i, "dead"] <- m_param_attributable[i, "dead"] +
  #       (m_param_attributable[i, "mort_s"] * m_param_attributable[i, "sus"]) + 
  #       (m_param_attributable[i, "mort_r"] * m_param_attributable[i, "res"])
  #     
  #     m_param_attributable[i, "seq"] <- m_param_attributable[i, "seq"] +
  #       (m_param_attributable[i, "s_seq"] * m_param_attributable[i, "sus"]) + 
  #       (m_param_attributable[i, "r_seq"] * m_param_attributable[i, "res"])
  #     
  #   }
  #   
  #   return(m_param_attributable)
  # }
  # 
  # m_param_attributable <- f_human_epi_attributable(m_param_attributable,n.t) 
  
  
  # Animal Epi Model --------------------------------------------------------
  
  state_names_a <- c("well", "res","sus","fallen","sold") ## the compartments
  transition_names_a  <- c("birth","r","s","mort_r", "mort_s","mort_w", "rec_r","rec_s","w_sold")  ## the rates
  parameter_names_a <- c(state_names_a, transition_names_a)
  
  f_animal_epi <- function(m_param_a_base, n.t, scenario_animal){
    
    n.t.val <- n.t
    
    if(scenario_animal == "chicken_small"){
      n_animals_farm <- inputs[parameter == "n_chickens_farm_small", Value]
      annual_cycles <- inputs[parameter == "pcycles_chicken_small", Value]
    } else if (scenario_animal == "chicken_ind"){
      n_animals_farm <- inputs[parameter == "n_chickens_farm_ind", Value]
      annual_cycles <- inputs[parameter == "pcycles_chicken_ind", Value]
    } else if (scenario_animal == "pig_small"){
      n_animals_farm <- inputs[parameter == "n_pigs_farm_small", Value]
      annual_cycles <- inputs[parameter == "pcycles_pig_small", Value]
    } else if (scenario_animal == "pig_ind"){
      n_animals_farm <- inputs[parameter == "n_pigs_farm_ind", Value]
      annual_cycles <- inputs[parameter == "pcycles_pig_ind", Value]  
    }    
    
    m_param_a_temp <- m_param_a_base[1:4,]
    rownames(m_param_a_temp) <- NULL ##removing rownames
    
    i <- 2 ## getting sick
    m_param_a_temp[i,"well"] <- m_param_a_temp[i-1,"well"] -(m_param_a_temp[i-1,"r"]*m_param_a_temp[i-1,"well"]) -
      (m_param_a_temp[i-1,"s"]*m_param_a_temp[i-1,"well"])
    m_param_a_temp[i,"res"] <- m_param_a_temp[i-1,"res"] + (m_param_a_temp[i-1,"r"]*m_param_a_temp[i-1,"well"]) 
    m_param_a_temp[i,"sus"] <- m_param_a_temp[i-1,"sus"] + (m_param_a_temp[i-1,"s"]*m_param_a_temp[i-1,"well"])
    
    i <- 3 ## dying and recovering
    m_param_a_temp[i,"fallen"] <- (m_param_a_temp[i-1,"mort_w"]*m_param_a_temp[i-1,"well"]) +
      (m_param_a_temp[i-1,"mort_r"]*m_param_a_temp[i-1,"res"]) + 
      (m_param_a_temp[i-1,"mort_s"]*m_param_a_temp[i-1,"sus"])
    m_param_a_temp[i, "well"] <- m_param_a_temp[i-1,"well"] - (m_param_a_temp[i-1,"well"]*m_param_a_temp[i-1,"mort_w"]) +
      (m_param_a_temp[i-1,"rec_r"]*m_param_a_temp[i-1,"res"])+ 
      (m_param_a_temp[i-1,"rec_s"]*m_param_a_temp[i-1,"sus"])
    
    i <- 4 ##sold
    m_param_a_temp[i,"sold"] <- (m_param_a_temp[i-1,"w_sold"]*m_param_a_temp[i-1,"well"])
    
    ## final states
    m_a_sum <- m_param_a_temp[4,]
    m_a_sum["res"] <- m_param_a_temp[2,"res"]
    m_a_sum["sus"] <- m_param_a_temp[2,"sus"]
    m_a_sum["well"] <- n_animals_farm - m_a_sum["res"] - m_a_sum["sus"] ## reset the number in 'well'
    
    m_a_sum[1:5] <- annual_cycles * m_a_sum[1:5] #multiply by the number of annual cycles
    
    m_param_a <- matrix(rep(m_a_sum), nrow=n.t.val, ncol =length(parameter_names_a))
    m_param_a <- t(replicate(n.t.val,m_a_sum))
    colnames(m_param_a) <- parameter_names_a
    rownames(m_param_a) <- paste("cycle", 0:(n.t-1), sep  =  "")
    
    #adjust farm outputs for the growth in agricultural output
    for(i in 1:nrow(m_param_a)){
      m_param_a[i,] <- m_param_a[i,] * pop_relative[i]
    }
    
    m_param_a
    
    return(m_param_a)
    
  }
  
  # Chicken Epi Module - Smallholder -------------------------------------------------------
  
  scenario_animal <- "chicken_small"
  
  state_i_c_s <- c(inputs[parameter == "n_chickens_farm_small", Value], rep(0,length=length(state_names_a)-1))
  
  #make the regular smallholder chicken matrix
  m_param_c_s_base <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_c_s_base) <- parameter_names_a
  rownames(m_param_c_s_base) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  m_param_c_s_base[ , "r"] <- rep(inputs[parameter=="c_res_small", Value], n.t)
  m_param_c_s_base[ , "s"] <- rep(inputs[parameter=="c_sus_small", Value], n.t)
  m_param_c_s_base[ , "mort_s"] <- rep(inputs[parameter=="c_mort_small_sus", Value], n.t) 
  m_param_c_s_base[ , "mort_r"] <- rep(inputs[parameter=="c_mort_small_res", Value], n.t)
  m_param_c_s_base[ , "rec_r"] <- rep(1-(m_param_c_s_base[1,"mort_r"]), n.t)
  m_param_c_s_base[ , "rec_s"] <- rep(1-(m_param_c_s_base[1,"mort_s"]), n.t)
  m_param_c_s_base[ , "birth"] <- rep(1, n.t)
  m_param_c_s_base[ , "mort_w"] <- rep(inputs[parameter=="c_mort_small_well", Value], n.t)
  m_param_c_s_base[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_c_s_base[i, "r"] > max_res * (m_param_c_s_base[1, "r"] + m_param_c_s_base[1, "s"])){
      m_param_c_s_base[i, "r"] <- max_res * (m_param_c_s_base[1, "r"] + m_param_c_s_base[1, "s"])
    }
    m_param_c_s_base[i, "s"] <- m_param_c_s_base[1, "r"] + m_param_c_s_base[1, "s"] - m_param_c_s_base[i, "r"]
  }
  
  m_param_c_s_base[1, 1:length(state_names_a)] <- state_i_c_s
  
  m_param_c_s <- f_animal_epi(m_param_c_s_base,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals
  
  
  
  
  
  
  #make the attributable smallholder chicken matrix
  m_param_c_s_base_attributable <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_c_s_base_attributable) <- parameter_names_a
  rownames(m_param_c_s_base_attributable) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #m_param_c_s_base_attributable[ , "r"] <- rep(inputs[parameter=="c_res_small", Value], n.t)
  m_param_c_s_base_attributable[ , "s"] <- rep(inputs[parameter=="c_sus_small", Value] + inputs[parameter=="c_res_small", Value], n.t)
  m_param_c_s_base_attributable[ , "mort_s"] <- rep(inputs[parameter=="c_mort_small_sus", Value], n.t) 
  m_param_c_s_base_attributable[ , "mort_r"] <- rep(inputs[parameter=="c_mort_small_res", Value], n.t)
  m_param_c_s_base_attributable[ , "rec_r"] <- rep(1-(m_param_c_s_base_attributable[1,"mort_r"]), n.t)
  m_param_c_s_base_attributable[ , "rec_s"] <- rep(1-(m_param_c_s_base_attributable[1,"mort_s"]), n.t)
  m_param_c_s_base_attributable[ , "birth"] <- rep(1, n.t)
  m_param_c_s_base_attributable[ , "mort_w"] <- rep(inputs[parameter=="c_mort_small_well", Value], n.t)
  m_param_c_s_base_attributable[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_c_s_base_attributable[i, "r"] > max_res * (m_param_c_s_base_attributable[1, "r"] + m_param_c_s_base_attributable[1, "s"])){
      m_param_c_s_base_attributable[i, "r"] <- max_res * (m_param_c_s_base_attributable[1, "r"] + m_param_c_s_base_attributable[1, "s"])
    }
    m_param_c_s_base_attributable[i, "s"] <- m_param_c_s_base_attributable[1, "r"] + m_param_c_s_base_attributable[1, "s"] - m_param_c_s_base_attributable[i, "r"]
  }
  
  m_param_c_s_base_attributable[1, 1:length(state_names_a)] <- state_i_c_s
  
  m_param_c_s_attributable <- f_animal_epi(m_param_c_s_base_attributable,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals
  
  
  
  
  
  
  
  #make the associated smallholder chicken matrix
  m_param_c_s_base_associated <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_c_s_base_associated) <- parameter_names_a
  rownames(m_param_c_s_base_associated) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #m_param_c_s_base_associated[ , "r"] <- rep(inputs[parameter=="c_res_small", Value], n.t)
  m_param_c_s_base_associated[ , "s"] <- rep(inputs[parameter=="c_sus_small", Value], n.t)
  m_param_c_s_base_associated[ , "mort_s"] <- rep(inputs[parameter=="c_mort_small_sus", Value], n.t) 
  m_param_c_s_base_associated[ , "mort_r"] <- rep(inputs[parameter=="c_mort_small_res", Value], n.t)
  m_param_c_s_base_associated[ , "rec_r"] <- rep(1-(m_param_c_s_base_associated[1,"mort_r"]), n.t)
  m_param_c_s_base_associated[ , "rec_s"] <- rep(1-(m_param_c_s_base_associated[1,"mort_s"]), n.t)
  m_param_c_s_base_associated[ , "birth"] <- rep(1, n.t)
  m_param_c_s_base_associated[ , "mort_w"] <- rep(inputs[parameter=="c_mort_small_well", Value], n.t)
  m_param_c_s_base_associated[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_c_s_base_associated[i, "r"] > max_res * (m_param_c_s_base_associated[1, "r"] + m_param_c_s_base_associated[1, "s"])){
      m_param_c_s_base_associated[i, "r"] <- max_res * (m_param_c_s_base_associated[1, "r"] + m_param_c_s_base_associated[1, "s"])
    }
    m_param_c_s_base_associated[i, "s"] <- m_param_c_s_base_associated[1, "r"] + m_param_c_s_base_associated[1, "s"] - m_param_c_s_base_associated[i, "r"]
  }
  
  m_param_c_s_base_associated[1, 1:length(state_names_a)] <- state_i_c_s
  
  m_param_c_s_associated <- f_animal_epi(m_param_c_s_base_associated,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals
  
  # Chicken Epi Module - Industrial -----------------------------------------
  
  scenario_animal <- "chicken_ind"
  
  state_i_c_i <- c(inputs[parameter == "n_chickens_farm_ind", Value], rep(0,length=length(state_names_a)-1))
  
  
  #make the regular industrial chicken matrix
  m_param_c_i_base <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_c_i_base) <- parameter_names_a
  rownames(m_param_c_i_base) <- paste("cycle", 0:(n.t-1), sep  =  "")

  m_param_c_i_base[ , "r"] <- rep(inputs[parameter=="c_res_ind", Value], n.t)
  m_param_c_i_base[ , "s"] <- rep(inputs[parameter=="c_sus_ind", Value], n.t)
  m_param_c_i_base[ , "mort_s"] <- rep(inputs[parameter=="c_mort_ind_sus", Value], n.t) 
  m_param_c_i_base[ , "mort_r"] <- rep(inputs[parameter=="c_mort_ind_res", Value], n.t)
  m_param_c_i_base[ , "rec_r"] <- rep(1-(m_param_c_i_base[1,"mort_r"]), n.t)
  m_param_c_i_base[ , "rec_s"] <- rep(1-(m_param_c_i_base[1,"mort_s"]), n.t)
  m_param_c_i_base[ , "birth"] <- rep(1, n.t)
  m_param_c_i_base[ , "mort_w"] <- rep(inputs[parameter=="c_mort_ind_well", Value], n.t)
  m_param_c_i_base[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_c_i_base[i, "r"] > max_res * (m_param_c_i_base[1, "r"] + m_param_c_i_base[1, "s"])){
      m_param_c_i_base[i, "r"] <- max_res * (m_param_c_i_base[1, "r"] + m_param_c_i_base[1, "s"])
    }
    m_param_c_i_base[i, "s"] <- m_param_c_i_base[1, "r"] + m_param_c_i_base[1, "s"] - m_param_c_i_base[i, "r"]
  }
  
  m_param_c_i_base[1, 1:length(state_names_a)] <- state_i_c_i
  
  m_param_c_i <- f_animal_epi(m_param_c_i_base,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals
  
  
  
  
  
  
  
  
  #make the attributable industrial chicken matrix
  
  m_param_c_i_base_attributable <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_c_i_base_attributable) <- parameter_names_a
  rownames(m_param_c_i_base_attributable) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #m_param_c_i_base_attributable[ , "r"] <- rep(inputs[parameter=="c_res_ind", Value], n.t)
  m_param_c_i_base_attributable[ , "s"] <- rep(inputs[parameter=="c_sus_ind", Value] + inputs[parameter=="c_res_ind", Value], n.t)
  m_param_c_i_base_attributable[ , "mort_s"] <- rep(inputs[parameter=="c_mort_ind_sus", Value], n.t) 
  m_param_c_i_base_attributable[ , "mort_r"] <- rep(inputs[parameter=="c_mort_ind_res", Value], n.t)
  m_param_c_i_base_attributable[ , "rec_r"] <- rep(1-(m_param_c_i_base_attributable[1,"mort_r"]), n.t)
  m_param_c_i_base_attributable[ , "rec_s"] <- rep(1-(m_param_c_i_base_attributable[1,"mort_s"]), n.t)
  m_param_c_i_base_attributable[ , "birth"] <- rep(1, n.t)
  m_param_c_i_base_attributable[ , "mort_w"] <- rep(inputs[parameter=="c_mort_ind_well", Value], n.t)
  m_param_c_i_base_attributable[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_c_i_base_attributable[i, "r"] > max_res * (m_param_c_i_base_attributable[1, "r"] + m_param_c_i_base_attributable[1, "s"])){
      m_param_c_i_base_attributable[i, "r"] <- max_res * (m_param_c_i_base_attributable[1, "r"] + m_param_c_i_base_attributable[1, "s"])
    }
    m_param_c_i_base_attributable[i, "s"] <- m_param_c_i_base_attributable[1, "r"] + m_param_c_i_base_attributable[1, "s"] - m_param_c_i_base_attributable[i, "r"]
  }
  
  m_param_c_i_base_attributable[1, 1:length(state_names_a)] <- state_i_c_i
  
  m_param_c_i_attributable <- f_animal_epi(m_param_c_i_base_attributable,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals
  
  
  
  
  
  
  
  #make the associated industrial chicken matrix
  
  m_param_c_i_base_associated <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_c_i_base_associated) <- parameter_names_a
  rownames(m_param_c_i_base_associated) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #m_param_c_i_base_associated[ , "r"] <- rep(inputs[parameter=="c_res_ind", Value], n.t)
  m_param_c_i_base_associated[ , "s"] <- rep(inputs[parameter=="c_sus_ind", Value], n.t)
  m_param_c_i_base_associated[ , "mort_s"] <- rep(inputs[parameter=="c_mort_ind_sus", Value], n.t) 
  m_param_c_i_base_associated[ , "mort_r"] <- rep(inputs[parameter=="c_mort_ind_res", Value], n.t)
  m_param_c_i_base_associated[ , "rec_r"] <- rep(1-(m_param_c_i_base_associated[1,"mort_r"]), n.t)
  m_param_c_i_base_associated[ , "rec_s"] <- rep(1-(m_param_c_i_base_associated[1,"mort_s"]), n.t)
  m_param_c_i_base_associated[ , "birth"] <- rep(1, n.t)
  m_param_c_i_base_associated[ , "mort_w"] <- rep(inputs[parameter=="c_mort_ind_well", Value], n.t)
  m_param_c_i_base_associated[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_c_i_base_associated[i, "r"] > max_res * (m_param_c_i_base_associated[1, "r"] + m_param_c_i_base_associated[1, "s"])){
      m_param_c_i_base_associated[i, "r"] <- max_res * (m_param_c_i_base_associated[1, "r"] + m_param_c_i_base_associated[1, "s"])
    }
    m_param_c_i_base_associated[i, "s"] <- m_param_c_i_base_associated[1, "r"] + m_param_c_i_base_associated[1, "s"] - m_param_c_i_base_associated[i, "r"]
  }
  
  m_param_c_i_base_associated[1, 1:length(state_names_a)] <- state_i_c_i
  
  m_param_c_i_associated <- f_animal_epi(m_param_c_i_base_associated,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals
  
  
  # Pig Epi Module - Smallholder --------------------------------------------
  
  scenario_animal <- "pig_small"
  
  state_i_p_s <- c(inputs[parameter == "n_pigs_farm_small", Value], rep(0,length=length(state_names_a)-1))
  
  
  #make the regular smallholder pig matrix
  m_param_p_s_base <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_p_s_base) <- parameter_names_a
  rownames(m_param_p_s_base) <- paste("cycle", 0:(n.t-1), sep  =  "")

  m_param_p_s_base[ , "r"] <- rep(inputs[parameter=="p_res_small", Value], n.t)
  m_param_p_s_base[ , "s"] <- rep(inputs[parameter=="p_sus_small", Value], n.t)
  m_param_p_s_base[ , "mort_s"] <- rep(inputs[parameter=="p_mort_small_sus", Value], n.t) 
  m_param_p_s_base[ , "mort_r"] <- rep(inputs[parameter=="p_mort_small_res", Value], n.t)
  m_param_p_s_base[ , "rec_r"] <- rep(1-(m_param_p_s_base[1,"mort_r"]), n.t)
  m_param_p_s_base[ , "rec_s"] <- rep(1-(m_param_p_s_base[1,"mort_s"]), n.t)
  m_param_p_s_base[ , "birth"] <- rep(1, n.t)
  m_param_p_s_base[ , "mort_w"] <- rep(inputs[parameter=="p_mort_small_well", Value], n.t)
  m_param_p_s_base[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_p_s_base[i, "r"] > max_res * (m_param_p_s_base[1, "r"] + m_param_p_s_base[1, "s"])){
      m_param_p_s_base[i, "r"] <- max_res * (m_param_p_s_base[1, "r"] + m_param_p_s_base[1, "s"])
    }
    m_param_p_s_base[i, "s"] <- m_param_p_s_base[1, "r"] + m_param_p_s_base[1, "s"] - m_param_p_s_base[i, "r"]
  }
  
  m_param_p_s_base[1, 1:length(state_names_a)] <- state_i_p_s
  
  m_param_p_s <- f_animal_epi(m_param_p_s_base,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals
  
  
  
  
  
  #make the attributable smallholder pig matrix
  m_param_p_s_base_attributable <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_p_s_base_attributable) <- parameter_names_a
  rownames(m_param_p_s_base_attributable) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #m_param_p_s_base_attributable[ , "r"] <- rep(inputs[parameter=="p_res_small", Value], n.t)
  m_param_p_s_base_attributable[ , "s"] <- rep(inputs[parameter=="p_sus_small", Value] + inputs[parameter=="p_res_small", Value], n.t)
  m_param_p_s_base_attributable[ , "mort_s"] <- rep(inputs[parameter=="p_mort_small_sus", Value], n.t) 
  m_param_p_s_base_attributable[ , "mort_r"] <- rep(inputs[parameter=="p_mort_small_res", Value], n.t)
  m_param_p_s_base_attributable[ , "rec_r"] <- rep(1-(m_param_p_s_base_attributable[1,"mort_r"]), n.t)
  m_param_p_s_base_attributable[ , "rec_s"] <- rep(1-(m_param_p_s_base_attributable[1,"mort_s"]), n.t)
  m_param_p_s_base_attributable[ , "birth"] <- rep(1, n.t)
  m_param_p_s_base_attributable[ , "mort_w"] <- rep(inputs[parameter=="p_mort_small_well", Value], n.t)
  m_param_p_s_base_attributable[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_p_s_base_attributable[i, "r"] > max_res * (m_param_p_s_base_attributable[1, "r"] + m_param_p_s_base_attributable[1, "s"])){
      m_param_p_s_base_attributable[i, "r"] <- max_res * (m_param_p_s_base_attributable[1, "r"] + m_param_p_s_base_attributable[1, "s"])
    }
    m_param_p_s_base_attributable[i, "s"] <- m_param_p_s_base_attributable[1, "r"] + m_param_p_s_base_attributable[1, "s"] - m_param_p_s_base_attributable[i, "r"]
  }
  
  m_param_p_s_base_attributable[1, 1:length(state_names_a)] <- state_i_p_s
  
  m_param_p_s_attributable <- f_animal_epi(m_param_p_s_base_attributable,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals
  
  
  
  #make the associated smallholder pig matrix
  m_param_p_s_base_associated <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_p_s_base_associated) <- parameter_names_a
  rownames(m_param_p_s_base_associated) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #m_param_p_s_base_associated[ , "r"] <- rep(inputs[parameter=="p_res_small", Value], n.t)
  m_param_p_s_base_associated[ , "s"] <- rep(inputs[parameter=="p_sus_small", Value], n.t)
  m_param_p_s_base_associated[ , "mort_s"] <- rep(inputs[parameter=="p_mort_small_sus", Value], n.t) 
  m_param_p_s_base_associated[ , "mort_r"] <- rep(inputs[parameter=="p_mort_small_res", Value], n.t)
  m_param_p_s_base_associated[ , "rec_r"] <- rep(1-(m_param_p_s_base_associated[1,"mort_r"]), n.t)
  m_param_p_s_base_associated[ , "rec_s"] <- rep(1-(m_param_p_s_base_associated[1,"mort_s"]), n.t)
  m_param_p_s_base_associated[ , "birth"] <- rep(1, n.t)
  m_param_p_s_base_associated[ , "mort_w"] <- rep(inputs[parameter=="p_mort_small_well", Value], n.t)
  m_param_p_s_base_associated[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_p_s_base_associated[i, "r"] > max_res * (m_param_p_s_base_associated[1, "r"] + m_param_p_s_base_associated[1, "s"])){
      m_param_p_s_base_associated[i, "r"] <- max_res * (m_param_p_s_base_associated[1, "r"] + m_param_p_s_base_associated[1, "s"])
    }
    m_param_p_s_base_associated[i, "s"] <- m_param_p_s_base_associated[1, "r"] + m_param_p_s_base_associated[1, "s"] - m_param_p_s_base_associated[i, "r"]
  }
  
  m_param_p_s_base_associated[1, 1:length(state_names_a)] <- state_i_p_s
  
  m_param_p_s_associated <- f_animal_epi(m_param_p_s_base_associated,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals
  
  
  
  
  
  # Pig Epi Module - Industrial ---------------------------------------------
  
  scenario_animal <- "pig_ind"
  
  state_i_p_i <- c(inputs[parameter == "n_pigs_farm_ind", Value], rep(0,length=length(state_names_a)-1))
  
  
  #make the regular industrial pig matrix
  m_param_p_i_base <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_p_i_base) <- parameter_names_a
  rownames(m_param_p_i_base) <- paste("cycle", 0:(n.t-1), sep  =  "")

  m_param_p_i_base[ , "r"] <- rep(inputs[parameter=="p_res_ind", Value], n.t)
  m_param_p_i_base[ , "s"] <- rep(inputs[parameter=="p_sus_ind", Value], n.t)
  m_param_p_i_base[ , "mort_s"] <- rep(inputs[parameter=="p_mort_ind_sus", Value], n.t) 
  m_param_p_i_base[ , "mort_r"] <- rep(inputs[parameter=="p_mort_ind_res", Value], n.t)
  m_param_p_i_base[ , "rec_r"] <- rep(1-(m_param_p_i_base[1,"mort_r"]), n.t)
  m_param_p_i_base[ , "rec_s"] <- rep(1-(m_param_p_i_base[1,"mort_s"]), n.t)
  m_param_p_i_base[ , "birth"] <- rep(1, n.t)
  m_param_p_i_base[ , "mort_w"] <- rep(inputs[parameter=="p_mort_ind_well", Value], n.t)
  m_param_p_i_base[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_p_i_base[i, "r"] > max_res * (m_param_p_i_base[1, "r"] + m_param_p_i_base[1, "s"])){
      m_param_p_i_base[i, "r"] <- max_res * (m_param_p_i_base[1, "r"] + m_param_p_i_base[1, "s"])
    }
    m_param_p_i_base[i, "s"] <- m_param_p_i_base[1, "r"] + m_param_p_i_base[1, "s"] - m_param_p_i_base[i, "r"]
  }
  
  m_param_p_i_base[1, 1:length(state_names_a)] <- state_i_p_i
  
  m_param_p_i <- f_animal_epi(m_param_p_i_base,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals  
  
  
  
  
  #make the attributable industrial pig matrix
  m_param_p_i_base_attributable <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_p_i_base_attributable) <- parameter_names_a
  rownames(m_param_p_i_base_attributable) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #m_param_p_i_base_attributable[ , "r"] <- rep(inputs[parameter=="p_res_ind", Value], n.t)
  m_param_p_i_base_attributable[ , "s"] <- rep(inputs[parameter=="p_sus_ind", Value] + inputs[parameter=="p_res_ind", Value], n.t)
  m_param_p_i_base_attributable[ , "mort_s"] <- rep(inputs[parameter=="p_mort_ind_sus", Value], n.t) 
  m_param_p_i_base_attributable[ , "mort_r"] <- rep(inputs[parameter=="p_mort_ind_res", Value], n.t)
  m_param_p_i_base_attributable[ , "rec_r"] <- rep(1-(m_param_p_i_base_attributable[1,"mort_r"]), n.t)
  m_param_p_i_base_attributable[ , "rec_s"] <- rep(1-(m_param_p_i_base_attributable[1,"mort_s"]), n.t)
  m_param_p_i_base_attributable[ , "birth"] <- rep(1, n.t)
  m_param_p_i_base_attributable[ , "mort_w"] <- rep(inputs[parameter=="p_mort_ind_well", Value], n.t)
  m_param_p_i_base_attributable[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_p_i_base_attributable[i, "r"] > max_res * (m_param_p_i_base_attributable[1, "r"] + m_param_p_i_base_attributable[1, "s"])){
      m_param_p_i_base_attributable[i, "r"] <- max_res * (m_param_p_i_base_attributable[1, "r"] + m_param_p_i_base_attributable[1, "s"])
    }
    m_param_p_i_base_attributable[i, "s"] <- m_param_p_i_base_attributable[1, "r"] + m_param_p_i_base_attributable[1, "s"] - m_param_p_i_base_attributable[i, "r"]
  }
  
  m_param_p_i_base_attributable[1, 1:length(state_names_a)] <- state_i_p_i
  
  m_param_p_i_attributable <- f_animal_epi(m_param_p_i_base_attributable,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals  
  
  
  
  
  
  #make the associated industrial pig matrix
  m_param_p_i_base_associated <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_param_p_i_base_associated) <- parameter_names_a
  rownames(m_param_p_i_base_associated) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #m_param_p_i_base_associated[ , "r"] <- rep(inputs[parameter=="p_res_ind", Value], n.t)
  m_param_p_i_base_associated[ , "s"] <- rep(inputs[parameter=="p_sus_ind", Value], n.t)
  m_param_p_i_base_associated[ , "mort_s"] <- rep(inputs[parameter=="p_mort_ind_sus", Value], n.t) 
  m_param_p_i_base_associated[ , "mort_r"] <- rep(inputs[parameter=="p_mort_ind_res", Value], n.t)
  m_param_p_i_base_associated[ , "rec_r"] <- rep(1-(m_param_p_i_base_associated[1,"mort_r"]), n.t)
  m_param_p_i_base_associated[ , "rec_s"] <- rep(1-(m_param_p_i_base_associated[1,"mort_s"]), n.t)
  m_param_p_i_base_associated[ , "birth"] <- rep(1, n.t)
  m_param_p_i_base_associated[ , "mort_w"] <- rep(inputs[parameter=="p_mort_ind_well", Value], n.t)
  m_param_p_i_base_associated[ , "w_sold"] <- rep(1, n.t) 
  
  #make it so that the total incidence of infections stays the same, and only the
  #portion of them that are resistant changes
  for(i in 1:n.t){
    if(m_param_p_i_base_associated[i, "r"] > max_res * (m_param_p_i_base_associated[1, "r"] + m_param_p_i_base_associated[1, "s"])){
      m_param_p_i_base_associated[i, "r"] <- max_res * (m_param_p_i_base_associated[1, "r"] + m_param_p_i_base_associated[1, "s"])
    }
    m_param_p_i_base_associated[i, "s"] <- m_param_p_i_base_associated[1, "r"] + m_param_p_i_base_associated[1, "s"] - m_param_p_i_base_associated[i, "r"]
  }
  
  m_param_p_i_base_associated[1, 1:length(state_names_a)] <- state_i_p_i
  
  m_param_p_i_associated <- f_animal_epi(m_param_p_i_base_associated,n.t, scenario_animal)
  ### ignore totals of transition probs etc. as they are over counted etc.
  ## just want to focus on health state totals  
  
  
  
  
  
  #'All of the farm costs are zero here, as any changes to profits per animal
  #'from the intervention will be reflected in changes to the 'farm rewards' 
  
  # Chicken Farm Costs - Smallholder ----------------------------------------
  
  m_cost_c_s <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_cost_c_s) <- parameter_names_a
  rownames(m_cost_c_s) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  # Chicken Farm Costs - Industrial -----------------------------------------
  
  m_cost_c_i <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_cost_c_i) <- parameter_names_a
  rownames(m_cost_c_i) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  # Pig Farm Costs - Smallholder --------------------------------------------
  
  m_cost_p_s <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_cost_p_s) <- parameter_names_a
  rownames(m_cost_p_s) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  # Pig Farm Costs - Industrial ---------------------------------------------
  
  m_cost_p_i <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_cost_p_i) <- parameter_names_a
  rownames(m_cost_p_i) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  
  # Chicken Farm Rewards - Smallholder --------------------------------------
  
  m_rwd_c_s <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_rwd_c_s) <- parameter_names_a
  rownames(m_rwd_c_s) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #only get a reward for selling animals
  r_sold_c_s <- inputs[parameter == "chicken_weight", Value] * inputs[parameter == "chicken_price", Value]
  rwd_i_c_s <- c(0,0,0,0,r_sold_c_s)
  
  ## start at cycle 1 so you do not multiply initial state vector 
  m_rwd_c_s[2, 1:length(state_names_a)] <- rwd_i_c_s
  
  #discount
  for (j in 1:length(state_names_a)) {
    for (i in 3:(n.t)){
      m_rwd_c_s[i,j] <- f_di(m_rwd_c_s[i-1,j],dr)
    }  
  } 
  
  # Chicken Farm Rewards - Industrial ---------------------------------------
  
  m_rwd_c_i <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_rwd_c_i) <- parameter_names_a
  rownames(m_rwd_c_i) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #only get a reward for selling animals
  r_sold_c_i <- inputs[parameter == "chicken_weight", Value] * inputs[parameter == "chicken_price", Value]
  rwd_i_c_i <- c(0,0,0,0,r_sold_c_i)
  
  ## start at cycle 1 so you do not multiply initial state vector 
  m_rwd_c_i[2, 1:length(state_names_a)] <- rwd_i_c_i
  
  #discount
  for (j in 1:length(state_names_a)) {
    for (i in 3:(n.t)){
      m_rwd_c_i[i,j] <- f_di(m_rwd_c_i[i-1,j],dr)
    }  
  } 
  
  # Pig Farm Rewards - Smallholder --------------------------------------
  
  m_rwd_p_s <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_rwd_p_s) <- parameter_names_a
  rownames(m_rwd_p_s) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #only get a reward for selling animals
  r_sold_p_s <- inputs[parameter == "pig_weight", Value] * inputs[parameter == "pig_price", Value]
  rwd_i_p_s <- c(0,0,0,0,r_sold_p_s)
  
  ## start at cycle 1 so you do not multiply initial state vector 
  m_rwd_p_s[2, 1:length(state_names_a)] <- rwd_i_p_s
  
  #discount
  for (j in 1:length(state_names_a)) {
    for (i in 3:(n.t)){
      m_rwd_p_s[i,j] <- f_di(m_rwd_p_s[i-1,j],dr)
    }  
  } 
  
  # Pig Farm Rewards - Industrial ---------------------------------------
  
  m_rwd_p_i <- matrix(rep(0), nrow=n.t, ncol =length(parameter_names_a))
  colnames(m_rwd_p_i) <- parameter_names_a
  rownames(m_rwd_p_i) <- paste("cycle", 0:(n.t-1), sep  =  "")
  
  #only get a reward for selling animals
  r_sold_p_i <- inputs[parameter == "pig_weight", Value] * inputs[parameter == "pig_price", Value]
  rwd_i_p_i <- c(0,0,0,0,r_sold_p_i)
  
  ## start at cycle 1 so you do not multiply initial state vector 
  m_rwd_p_i[2, 1:length(state_names_a)] <- rwd_i_p_i
  
  #discount
  for (j in 1:length(state_names_a)) {
    for (i in 3:(n.t)){
      m_rwd_p_i[i,j] <- f_di(m_rwd_p_i[i-1,j],dr)
    }  
  }  
  



# Results matrices --------------------------------------------------------

  #Get number of each type of farm nationally
  n_farms_chicken_small <- inputs[parameter == "n_chickens", Value] * 
    (1 - inputs[parameter == "portion_animals_ind", Value]) / 
    inputs[parameter == "n_chickens_farm_small", Value]
  
  n_farms_chicken_ind <- inputs[parameter == "n_chickens", Value] * 
    (inputs[parameter == "portion_animals_ind", Value]) / 
    inputs[parameter == "n_chickens_farm_ind", Value]
  
  n_farms_pig_small <- inputs[parameter == "n_pigs", Value] * 
    (1 - inputs[parameter == "portion_animals_ind", Value]) / 
    inputs[parameter == "n_pigs_farm_small", Value]
  
  n_farms_pig_ind <- inputs[parameter == "n_pigs", Value] * 
    (inputs[parameter == "portion_animals_ind", Value]) / 
    inputs[parameter == "n_pigs_farm_ind", Value]
  
  #results matrices for healthcare
  results_base_h <- f_expvalue(m_param,m_cost,m_rwd)
  results_attributable_h <- f_expvalue(m_param_attributable,m_cost,m_rwd)
  results_associated_h <- f_expvalue(m_param_associated,m_cost,m_rwd)
  
  #results matrices for productivity
  results_base_prod <- f_expvalue(m_param,m_cost_prod,m_rwd_prod)
  results_attributable_prod <- f_expvalue(m_param_attributable,m_cost_prod,m_rwd_prod)
  results_associated_prod <- f_expvalue(m_param_associated,m_cost_prod,m_rwd_prod)
  
  #results matrix for industrial chicken farms
  results_base_c_i <- f_expvalue(m_param_c_i,m_cost_c_i,m_rwd_c_i)
  results_attributable_c_i <- f_expvalue(m_param_c_i_attributable,m_cost_c_i,m_rwd_c_i)
  results_associated_c_i <- f_expvalue(m_param_c_i_associated,m_cost_c_i,m_rwd_c_i)
  
  #results matrix for smallholder chicken farms
  results_base_c_s <- f_expvalue(m_param_c_s,m_cost_c_s,m_rwd_c_s)
  results_attributable_c_s <- f_expvalue(m_param_c_s_attributable,m_cost_c_s,m_rwd_c_s)
  results_associated_c_s <- f_expvalue(m_param_c_s_associated,m_cost_c_s,m_rwd_c_s)
  
  #results matrix for industrial pig farms
  results_base_p_i <- f_expvalue(m_param_p_i,m_cost_p_i,m_rwd_p_i)
  results_attributable_p_i <- f_expvalue(m_param_p_i_attributable,m_cost_p_i,m_rwd_p_i)
  results_associated_p_i <- f_expvalue(m_param_p_i_associated,m_cost_p_i,m_rwd_p_i)
  
  #results matrix for smallholder pig farms
  results_base_p_s <- f_expvalue(m_param_p_s,m_cost_p_s,m_rwd_p_s)
  results_attributable_p_s <- f_expvalue(m_param_p_s_attributable,m_cost_p_s,m_rwd_p_s)
  results_associated_p_s <- f_expvalue(m_param_p_s_associated,m_cost_p_s,m_rwd_p_s)

# Associated burden -------------------------------------------------------
  
  associated_burden_HC <- matrix(rep(0), nrow=2, ncol=2)
  colnames(associated_burden_HC) <- c("Costs ($)", "QALYs")
  rownames(associated_burden_HC) <- c("Base Case", "Associated Burden Counterfactual")
  
  associated_burden_HC[1,] <- results_base_h[1,]
  associated_burden_HC[2,] <- results_associated_h[1,]
  
  wtp <- inputs[parameter == "wtp", Value]
  
  incr_cost_health_associated <- (results_associated_h[1,1] - results_base_h[1,1])
  QALYs_saved_associated <-  (results_associated_h[1,2]-results_base_h[1,2])
  NMB_health_associated <- (QALYs_saved_associated*wtp)-(incr_cost_health_associated)
  
  
  
  associated_burden_prod <- matrix(rep(0), nrow=2, ncol=2)
  colnames(associated_burden_prod) <- c("Productivity", "QALYs")
  rownames(associated_burden_prod) <- c("Base Case", "Associated Burden Counterfactual")
  
  associated_burden_prod[1,2] <- results_base_h[,2]
  associated_burden_prod[2,2] <- results_associated_h[,2]
  associated_burden_prod[1,1] <- results_base_prod[,2]   #will be negative
  associated_burden_prod[2,1] <- results_associated_prod[,2] #will be negative but hopefully closer to zero
  
  incr_cost_prod_associated <- associated_burden_prod[1,1] - associated_burden_prod[2,1] #hopefully negative
  incr_benefit_prod_associated <- associated_burden_prod[2,2] - associated_burden_prod[1,2] #hopefully positive
  NMB_prod_associated <- associated_burden_prod[2,1] - associated_burden_prod[1,1] #hopefully positive
  
  
  
  associated_burden_c_i<- matrix(rep(0), nrow=2, ncol=2)
  colnames(associated_burden_c_i) <- c("Costs ($)", "Benefits ($)")
  rownames(associated_burden_c_i) <- c("Base Case", "Associated Burden Counterfactual")
  
  incr_cost_c_i_associated <- (results_associated_c_i[1,1] - results_base_c_i[1,1])
  incr_benefit_c_i_associated <-  (results_associated_c_i[1,2] - results_base_c_i[1,2])
   
  associated_burden_c_i[1,] <- results_base_c_i[1,] 
  associated_burden_c_i[2,] <- results_associated_c_i[1,] 
  
  NMB_c_i_associated <- (incr_benefit_c_i_associated - incr_cost_c_i_associated) * n_farms_chicken_ind
  
  
  
  
  associated_burden_c_s<- matrix(rep(0), nrow=2, ncol=2)
  colnames(associated_burden_c_s) <- c("Costs ($)", "Benefits ($)")
  rownames(associated_burden_c_s) <- c("Base Case", "Associated Burden Counterfactual")
  
  incr_cost_c_s_associated <- (results_associated_c_s[1,1] - results_base_c_s[1,1])
  incr_benefit_c_s_associated <-  (results_associated_c_s[1,2] - results_base_c_s[1,2])
  
  associated_burden_c_s[1,] <- results_base_c_s[1,] 
  associated_burden_c_s[2,] <- results_associated_c_s[1,] 
  
  NMB_c_s_associated <- (incr_benefit_c_s_associated - incr_cost_c_s_associated) * n_farms_chicken_small
  
  
  
  
  associated_burden_p_i<- matrix(rep(0), nrow=2, ncol=2)
  colnames(associated_burden_p_i) <- c("Costs ($)", "Benefits ($)")
  rownames(associated_burden_p_i) <- c("Base Case", "Associated Burden Counterfactual")
  
  incr_cost_p_i_associated <- (results_associated_p_i[1,1] - results_base_p_i[1,1])
  incr_benefit_p_i_associated <-  (results_associated_p_i[1,2] - results_base_p_i[1,2])
  
  associated_burden_p_i[1,] <- results_base_p_i[1,] 
  associated_burden_p_i[2,] <- results_associated_p_i[1,] 
  
  NMB_p_i_associated <- (incr_benefit_p_i_associated - incr_cost_p_i_associated) * n_farms_pig_ind
  
  
  
  
  associated_burden_p_s<- matrix(rep(0), nrow=2, ncol=2)
  colnames(associated_burden_p_s) <- c("Costs ($)", "Benefits ($)")
  rownames(associated_burden_p_s) <- c("Base Case", "Associated Burden Counterfactual")
  
  incr_cost_p_s_associated <- (results_associated_p_s[1,1] - results_base_p_s[1,1])
  incr_benefit_p_s_associated <-  (results_associated_p_s[1,2] - results_base_p_s[1,2])
  
  associated_burden_p_s[1,] <- results_base_p_s[1,] 
  associated_burden_p_s[2,] <- results_associated_p_s[1,] 
  
  NMB_p_s_associated <- (incr_benefit_p_s_associated - incr_cost_p_s_associated) * n_farms_pig_small
  
  
  
  
  total_associated_burden_discounted <- NMB_c_i_associated + NMB_c_s_associated +
    NMB_p_i_associated + NMB_p_s_associated + NMB_health_associated + NMB_prod_associated
  
  #convert maximum total cost to maximum annual cost
  discount_vector <- rep(0,n.t)
  discount_vector[1] <- 1
  for(i in 2:length(discount_vector)){
    discount_vector[i] = discount_vector[i-1]*(1-dr)
  }
  discount_sum <- sum(discount_vector)
  
  annual_equivalent_associated_burden <- total_associated_burden_discounted / discount_sum
  
  money_saved_health_associated <- -1 * incr_cost_health_associated
  valuation_QALYs_associated <- QALYs_saved_associated * wtp
  
  #Final outputs
  outputs_associated <- data.table("Associated Burden over Time Period" = total_associated_burden_discounted,
                                   "Equivalent Annual Burden"= annual_equivalent_associated_burden,
                                   "Burden to Productivtiy"=NMB_prod_associated,
                                   "Burden to Healthcare Sector"=money_saved_health_associated,
                                   "Burden to Human Life Years"=valuation_QALYs_associated,
                                   "Burden to Smallholder Pig Farms"=NMB_p_s_associated,
                                   "Burden to Industrial Pig Farms"=NMB_p_i_associated,
                                   "Burden to Smallholder Chicken Farms"=NMB_c_s_associated,
                                   "Burden to Industrial Chicken Farms"=NMB_c_i_associated,
                                   "Burden to Agriculture (Total)" = NMB_p_s_associated + NMB_p_i_associated + NMB_c_s_associated + NMB_c_i_associated,
                                   "QALYs Lost"=QALYs_saved_associated)
  outputs_associated
  
# Attributable burden -----------------------------------------------------

  attributable_burden_HC <- matrix(rep(0), nrow=2, ncol=2)
  colnames(attributable_burden_HC) <- c("Costs ($)", "QALYs")
  rownames(attributable_burden_HC) <- c("Base Case", "Attributable Burden Counterfactual")
  
  attributable_burden_HC[1,] <- results_base_h[1,]
  attributable_burden_HC[2,] <- results_attributable_h[1,]
  
  wtp <- inputs[parameter == "wtp", Value]
  
  incr_cost_health_attributable <- (results_attributable_h[1,1] - results_base_h[1,1])
  QALYs_saved_attributable <-  (results_attributable_h[1,2]-results_base_h[1,2])
  NMB_health_attributable <- (QALYs_saved_attributable*wtp)-(incr_cost_health_attributable)
  
  
  
  attributable_burden_prod <- matrix(rep(0), nrow=2, ncol=2)
  colnames(attributable_burden_prod) <- c("Productivity", "QALYs")
  rownames(attributable_burden_prod) <- c("Base Case", "Attributable Burden Counterfactual")
  
  attributable_burden_prod[1,2] <- results_base_h[,2]
  attributable_burden_prod[2,2] <- results_attributable_h[,2]
  attributable_burden_prod[1,1] <- results_base_prod[,2]   #will be negative
  attributable_burden_prod[2,1] <- results_attributable_prod[,2] #will be negative but hopefully closer to zero
  
  incr_cost_prod_attributable <- attributable_burden_prod[1,1] - attributable_burden_prod[2,1] #hopefully negative
  incr_benefit_prod_attributable <- attributable_burden_prod[2,2] - attributable_burden_prod[1,2] #hopefully positive
  NMB_prod_attributable <- attributable_burden_prod[2,1] - attributable_burden_prod[1,1] #hopefully positive
  
  
  
  attributable_burden_c_i<- matrix(rep(0), nrow=2, ncol=2)
  colnames(attributable_burden_c_i) <- c("Costs ($)", "Benefits ($)")
  rownames(attributable_burden_c_i) <- c("Base Case", "Attributable Burden Counterfactual")
  
  incr_cost_c_i_attributable <- (results_attributable_c_i[1,1] - results_base_c_i[1,1])
  incr_benefit_c_i_attributable <-  (results_attributable_c_i[1,2] - results_base_c_i[1,2])
  
  attributable_burden_c_i[1,] <- results_base_c_i[1,] 
  attributable_burden_c_i[2,] <- results_attributable_c_i[1,] 
  
  NMB_c_i_attributable <- (incr_benefit_c_i_attributable - incr_cost_c_i_attributable) * n_farms_chicken_ind
  
  
  
  
  attributable_burden_c_s<- matrix(rep(0), nrow=2, ncol=2)
  colnames(attributable_burden_c_s) <- c("Costs ($)", "Benefits ($)")
  rownames(attributable_burden_c_s) <- c("Base Case", "Attributable Burden Counterfactual")
  
  incr_cost_c_s_attributable <- (results_attributable_c_s[1,1] - results_base_c_s[1,1])
  incr_benefit_c_s_attributable <-  (results_attributable_c_s[1,2] - results_base_c_s[1,2])
  
  attributable_burden_c_s[1,] <- results_base_c_s[1,] 
  attributable_burden_c_s[2,] <- results_attributable_c_s[1,] 
  
  NMB_c_s_attributable <- (incr_benefit_c_s_attributable - incr_cost_c_s_attributable) * n_farms_chicken_small
  
  
  
  
  attributable_burden_p_i<- matrix(rep(0), nrow=2, ncol=2)
  colnames(attributable_burden_p_i) <- c("Costs ($)", "Benefits ($)")
  rownames(attributable_burden_p_i) <- c("Base Case", "Attributable Burden Counterfactual")
  
  incr_cost_p_i_attributable <- (results_attributable_p_i[1,1] - results_base_p_i[1,1])
  incr_benefit_p_i_attributable <-  (results_attributable_p_i[1,2] - results_base_p_i[1,2])
  
  attributable_burden_p_i[1,] <- results_base_p_i[1,] 
  attributable_burden_p_i[2,] <- results_attributable_p_i[1,] 
  
  NMB_p_i_attributable <- (incr_benefit_p_i_attributable - incr_cost_p_i_attributable) * n_farms_pig_ind
  
  
  
  
  attributable_burden_p_s<- matrix(rep(0), nrow=2, ncol=2)
  colnames(attributable_burden_p_s) <- c("Costs ($)", "Benefits ($)")
  rownames(attributable_burden_p_s) <- c("Base Case", "Attributable Burden Counterfactual")
  
  incr_cost_p_s_attributable <- (results_attributable_p_s[1,1] - results_base_p_s[1,1])
  incr_benefit_p_s_attributable <-  (results_attributable_p_s[1,2] - results_base_p_s[1,2])
  
  attributable_burden_p_s[1,] <- results_base_p_s[1,] 
  attributable_burden_p_s[2,] <- results_attributable_p_s[1,] 
  
  NMB_p_s_attributable <- (incr_benefit_p_s_attributable - incr_cost_p_s_attributable) * n_farms_pig_small
  
  
  
  
  total_attributable_burden_discounted <- NMB_c_i_attributable + NMB_c_s_attributable +
    NMB_p_i_attributable + NMB_p_s_attributable + NMB_health_attributable + NMB_prod_attributable
  
  #convert maximum total cost to maximum annual cost
  discount_vector <- rep(0,n.t)
  discount_vector[1] <- 1
  for(i in 2:length(discount_vector)){
    discount_vector[i] = discount_vector[i-1]*(1-dr)
  }
  discount_sum <- sum(discount_vector)
  
  annual_equivalent_attributable_burden <- total_attributable_burden_discounted / discount_sum
  
  money_saved_health_attributable <- -1 * incr_cost_health_attributable
  valuation_QALYs_attributable <- QALYs_saved_attributable * wtp
  
  #Final outputs
  outputs_attributable <- data.table("Attributable Burden over Time Period" = total_attributable_burden_discounted,
                                   "Equivalent Annual Burden"= annual_equivalent_attributable_burden,
                                   "Burden to Productivtiy"=NMB_prod_attributable,
                                   "Burden to Healthcare Sector"=money_saved_health_attributable,
                                   "Burden to Human Life Years"=valuation_QALYs_attributable,
                                   "Burden to Smallholder Pig Farms"=NMB_p_s_attributable,
                                   "Burden to Industrial Pig Farms"=NMB_p_i_attributable,
                                   "Burden to Smallholder Chicken Farms"=NMB_c_s_attributable,
                                   "Burden to Industrial Chicken Farms"=NMB_c_i_attributable,
                                   "Burden to Agriculture (Total)" = NMB_p_s_attributable + NMB_p_i_attributable + NMB_c_s_attributable + NMB_c_i_attributable,
                                   "QALYs Lost"=QALYs_saved_attributable)
  outputs_attributable

  # Final Output Display ----------------------------------------------------

  Outputs <- matrix(rep(0), nrow = ncol(outputs_associated), ncol = 3)
  colnames(Outputs) = c("Output", "Associated Burden", "Attributable Burden")
  Outputs[,1] <- names(outputs_associated) 
  Outputs[,2] <- as.numeric(outputs_associated)
  Outputs[,3] <-as.numeric(outputs_attributable)

  
  return(Outputs)

}



# Tornado Function --------------------------------------------------------



tornado_function <- function(inputs){
  
  
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
  
  for(i in c(3, 5:7, 10, 27:28, 38, 48)){
    
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
  
  # #do livestock industrialisation manually bc it didn't work
  # inputs_tornado <- inputs
  # 
  # inputs_tornado[48,3] <- inputs_tornado[48,4]
  # tornado_values[1,18] <- as.numeric(Model(inputs_tornado, "HCA")[1,2]) - tornado_values[1,1]
  # tornado_values[2,18] <- as.numeric(Model(inputs_tornado, "HCA")[1,3]) - tornado_values[2,1]
  # tornado_values[3,18] <- as.numeric(Model(inputs_tornado, "FCA")[1,2]) - tornado_values[3,1]
  # tornado_values[4,18] <- as.numeric(Model(inputs_tornado, "FCA")[1,3]) - tornado_values[4,1]
  # 
  # inputs_tornado[48,3] <- inputs_tornado[48,5]
  # tornado_values[1,19] <- as.numeric(Model(inputs_tornado, "HCA")[1,2]) - tornado_values[1,1]
  # tornado_values[2,19] <- as.numeric(Model(inputs_tornado, "HCA")[1,3]) - tornado_values[2,1]
  # tornado_values[3,19] <- as.numeric(Model(inputs_tornado, "FCA")[1,2]) - tornado_values[3,1]
  # tornado_values[4,19] <- as.numeric(Model(inputs_tornado, "FCA")[1,3]) - tornado_values[4,1]
  
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
                                       
                                       # "min" = c(tornado_values[1,2],
                                       #         tornado_values[1,4],
                                       #         tornado_values[1,6],
                                       #         tornado_values[1,8],
                                       #         tornado_values[1,10],
                                       #         tornado_values[1,12],
                                       #         tornado_values[1,14],
                                       #         tornado_values[1,16],
                                       #         tornado_values[1,18]),
                                       # 
                                       # "max" = c(tornado_values[1,3],
                                       #         tornado_values[1,5],
                                       #         tornado_values[1,7],
                                       #         tornado_values[1,9],
                                       #         tornado_values[1,11],
                                       #         tornado_values[1,13],
                                       #         tornado_values[1,15],
                                       #         tornado_values[1,17],
                                       #         tornado_values[1,19])
                                       
                                       "min" = c(
                                         min(tornado_values[1,2], tornado_values[1,3]),
                                         min(tornado_values[1,4], tornado_values[1,5]),
                                         min(tornado_values[1,6], tornado_values[1,7]),
                                         min(tornado_values[1,8], tornado_values[1,9]),
                                         min(tornado_values[1,10], tornado_values[1,11]),
                                         min(tornado_values[1,12], tornado_values[1,13]),
                                         min(tornado_values[1,14], tornado_values[1,15]),
                                         min(tornado_values[1,16], tornado_values[1,17]),
                                         min(tornado_values[1,18], tornado_values[1,19])
                                       ),
                                       
                                       "max" = c(
                                         max(tornado_values[1,2], tornado_values[1,3]),
                                         max(tornado_values[1,4], tornado_values[1,5]),
                                         max(tornado_values[1,6], tornado_values[1,7]),
                                         max(tornado_values[1,8], tornado_values[1,9]),
                                         max(tornado_values[1,10], tornado_values[1,11]),
                                         max(tornado_values[1,12], tornado_values[1,13]),
                                         max(tornado_values[1,14], tornado_values[1,15]),
                                         max(tornado_values[1,16], tornado_values[1,17]),
                                         max(tornado_values[1,18], tornado_values[1,19])
                                       )
                                       
  )
  
  tornado_HCA_attributable <- data.frame("variable" = c("Discount Rate",
                                                        "Population Growth Rate",
                                                        "Willingness-to-Pay Threshold",
                                                        "Productivity Growth Rate",
                                                        "Labour Force Participation Rate",
                                                        "Incidence of Bacterial Infection",
                                                        "Rate of AMR",
                                                        "AMR Growth Rate",
                                                        "Livestock Industrialisation"),
                                         
                                         # "min" = c(tornado_values[2,2],
                                         #           tornado_values[2,4],
                                         #           tornado_values[2,6],
                                         #           tornado_values[2,8],
                                         #           tornado_values[2,10],
                                         #           tornado_values[2,12],
                                         #           tornado_values[2,14],
                                         #           tornado_values[2,16],
                                         #           tornado_values[2,18]),
                                         # 
                                         # "max" = c(tornado_values[2,3],
                                         #           tornado_values[2,5],
                                         #           tornado_values[2,7],
                                         #           tornado_values[2,9],
                                         #           tornado_values[2,11],
                                         #           tornado_values[2,13],
                                         #           tornado_values[2,15],
                                         #           tornado_values[2,17],
                                         #           tornado_values[2,19])
                                         
                                         "min" = c(
                                           min(tornado_values[2,2], tornado_values[2,3]),
                                           min(tornado_values[2,4], tornado_values[2,5]),
                                           min(tornado_values[2,6], tornado_values[2,7]),
                                           min(tornado_values[2,8], tornado_values[2,9]),
                                           min(tornado_values[2,10], tornado_values[2,11]),
                                           min(tornado_values[2,12], tornado_values[2,13]),
                                           min(tornado_values[2,14], tornado_values[2,15]),
                                           min(tornado_values[2,16], tornado_values[2,17]),
                                           min(tornado_values[2,18], tornado_values[2,19])
                                         ),
                                         
                                         "max" = c(
                                           max(tornado_values[2,2], tornado_values[2,3]),
                                           max(tornado_values[2,4], tornado_values[2,5]),
                                           max(tornado_values[2,6], tornado_values[2,7]),
                                           max(tornado_values[2,8], tornado_values[2,9]),
                                           max(tornado_values[2,10], tornado_values[2,11]),
                                           max(tornado_values[2,12], tornado_values[2,13]),
                                           max(tornado_values[2,14], tornado_values[2,15]),
                                           max(tornado_values[2,16], tornado_values[2,17]),
                                           max(tornado_values[2,18], tornado_values[2,19])
                                         )
                                         
                                         
                                         
  )
  
  tornado_FCA_associated <- data.frame("variable" = c("Discount Rate",
                                                      "Population Growth Rate",
                                                      "Willingness-to-Pay Threshold",
                                                      "Productivity Growth Rate",
                                                      "Labour Force Participation Rate",
                                                      "Incidence of Bacterial Infection",
                                                      "Rate of AMR",
                                                      "AMR Growth Rate",
                                                      "Livestock Industrialisation"),
                                       
                                       # "min" = c(tornado_values[3,2],
                                       #           tornado_values[3,4],
                                       #           tornado_values[3,6],
                                       #           tornado_values[3,8],
                                       #           tornado_values[3,10],
                                       #           tornado_values[3,12],
                                       #           tornado_values[3,14],
                                       #           tornado_values[3,16],
                                       #           tornado_values[3,18]),
                                       # 
                                       # "max" = c(tornado_values[3,3],
                                       #           tornado_values[3,5],
                                       #           tornado_values[3,7],
                                       #           tornado_values[3,9],
                                       #           tornado_values[3,11],
                                       #           tornado_values[3,13],
                                       #           tornado_values[3,15],
                                       #           tornado_values[3,17],
                                       #           tornado_values[3,19])
                                       
                                       "min" = c(
                                         min(tornado_values[3,2], tornado_values[3,3]),
                                         min(tornado_values[3,4], tornado_values[3,5]),
                                         min(tornado_values[3,6], tornado_values[3,7]),
                                         min(tornado_values[3,8], tornado_values[3,9]),
                                         min(tornado_values[3,10], tornado_values[3,11]),
                                         min(tornado_values[3,12], tornado_values[3,13]),
                                         min(tornado_values[3,14], tornado_values[3,15]),
                                         min(tornado_values[3,16], tornado_values[3,17]),
                                         min(tornado_values[3,18], tornado_values[3,19])
                                       ),
                                       
                                       "max" = c(
                                         max(tornado_values[3,2], tornado_values[3,3]),
                                         max(tornado_values[3,4], tornado_values[3,5]),
                                         max(tornado_values[3,6], tornado_values[3,7]),
                                         max(tornado_values[3,8], tornado_values[3,9]),
                                         max(tornado_values[3,10], tornado_values[3,11]),
                                         max(tornado_values[3,12], tornado_values[3,13]),
                                         max(tornado_values[3,14], tornado_values[3,15]),
                                         max(tornado_values[3,16], tornado_values[3,17]),
                                         max(tornado_values[3,18], tornado_values[3,19])
                                       )
                                       
  )
  
  tornado_FCA_attributable <- data.frame("variable" = c("Discount Rate",
                                                        "Population Growth Rate",
                                                        "Willingness-to-Pay Threshold",
                                                        "Productivity Growth Rate",
                                                        "Labour Force Participation Rate",
                                                        "Incidence of Bacterial Infection",
                                                        "Rate of AMR",
                                                        "AMR Growth Rate",
                                                        "Livestock Industrialisation"),
                                         
                                         # "min" = c(tornado_values[4,2],
                                         #           tornado_values[4,4],
                                         #           tornado_values[4,6],
                                         #           tornado_values[4,8],
                                         #           tornado_values[4,10],
                                         #           tornado_values[4,12],
                                         #           tornado_values[4,14],
                                         #           tornado_values[4,16],
                                         #           tornado_values[4,18]),
                                         # 
                                         # "max" = c(tornado_values[4,3],
                                         #           tornado_values[4,5],
                                         #           tornado_values[4,7],
                                         #           tornado_values[4,9],
                                         #           tornado_values[4,11],
                                         #           tornado_values[4,13],
                                         #           tornado_values[4,15],
                                         #           tornado_values[4,17],
                                         #           tornado_values[4,19])
                                         
                                         "min" = c(
                                           min(tornado_values[4,2], tornado_values[4,3]),
                                           min(tornado_values[4,4], tornado_values[4,5]),
                                           min(tornado_values[4,6], tornado_values[4,7]),
                                           min(tornado_values[4,8], tornado_values[4,9]),
                                           min(tornado_values[4,10], tornado_values[4,11]),
                                           min(tornado_values[4,12], tornado_values[4,13]),
                                           min(tornado_values[4,14], tornado_values[4,15]),
                                           min(tornado_values[4,16], tornado_values[4,17]),
                                           min(tornado_values[4,18], tornado_values[4,19])
                                         ),
                                         
                                         "max" = c(
                                           max(tornado_values[4,2], tornado_values[4,3]),
                                           max(tornado_values[4,4], tornado_values[4,5]),
                                           max(tornado_values[4,6], tornado_values[4,7]),
                                           max(tornado_values[4,8], tornado_values[4,9]),
                                           max(tornado_values[4,10], tornado_values[4,11]),
                                           max(tornado_values[4,12], tornado_values[4,13]),
                                           max(tornado_values[4,14], tornado_values[4,15]),
                                           max(tornado_values[4,16], tornado_values[4,17]),
                                           max(tornado_values[4,18], tornado_values[4,19])
                                         )
                                         
  )
  
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
  
  tornado_df_HCA_associated <- tornado_HCA_associated %>%
    mutate(difference = abs(max - min))
  
  tornado_df_HCA_attributable <- tornado_HCA_attributable %>%
    mutate(difference = abs(max - min))
  
  tornado_df_FCA_associated <- tornado_FCA_associated %>%
    mutate(difference = abs(max - min))
  
  tornado_df_FCA_attributable <- tornado_FCA_attributable %>%
    mutate(difference = abs(max - min))
  
  tornado_df <- rbind(tornado_df_HCA_associated, 
                      tornado_df_HCA_attributable,
                      tornado_df_FCA_associated,
                      tornado_df_FCA_attributable)
  
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
  
  #tornado_df$lim <- max(abs(tornado_df$max), max(tornado_df$min))
  lim <- max(abs(tornado_df$max), max(tornado_df$min))
  lim <- lim / 1e+06
  
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
    scale_y_continuous(limits = c(-1.1 * lim, 1.1 * lim), "Millions 2019 $USD over 20 Years") +
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
    scale_y_continuous(limits = c(-1.1 * lim, 1.1 * lim), "Millions 2019 $USD over 20 Years") +
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
    scale_y_continuous(limits = c(-1.1 * lim, 1.1 * lim), "Millions 2019 $USD over 20 Years") +
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
    scale_y_continuous(limits = c(-1.1 * lim, 1.1 * lim), "Millions 2019 $USD over 20 Years") +
    ggtitle("Effect on Attributable Burden (Friction Cost Approach)")
  
  #create model output
  tornadoplot <- ggarrange(g1, g2, g3, g4,
                           ncol = 2,
                           nrow = 2)
  
  return(tornadoplot)
  
}



# table function ----------------------------------------------------------

table_function <- function(inputs){
  
  main_outputs <- matrix(rep(0), ncol = 5, nrow = 11)
  colnames(main_outputs) = c("Output", "Associated Burden (HCA)", "Associated Burden (FCA)",
                             "Attributable Burden (HCA)", "Attributable Burden (FCA)")
  
  main_outputs[,1] <- Model(inputs, "HCA")[,1]
  main_outputs[,2] <- as.numeric(Model(inputs, "HCA")[,2])
  main_outputs[,3] <- as.numeric(Model(inputs, "FCA")[,2])
  main_outputs[,4] <- as.numeric(Model(inputs, "HCA")[,3])
  main_outputs[,5] <- as.numeric(Model(inputs, "FCA")[,3])
  
  return(main_outputs)
  
}

# plot function -----------------------------------------------------------

plot_function <- function(inputs){
  
  
  #Table
  main_outputs <- matrix(rep(0), ncol = 5, nrow = 11)
  colnames(main_outputs) = c("Output", "Associated Burden (HCA)", "Associated Burden (FCA)",
                             "Attributable Burden (HCA)", "Attributable Burden (FCA)")
  
  main_outputs[,1] <- Model(inputs, "HCA")[,1]
  main_outputs[,2] <- as.numeric(Model(inputs, "HCA")[,2])
  main_outputs[,3] <- as.numeric(Model(inputs, "FCA")[,2])
  main_outputs[,4] <- as.numeric(Model(inputs, "HCA")[,3])
  main_outputs[,5] <- as.numeric(Model(inputs, "FCA")[,3])
  
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
  
  plot <- ggarrange(main_outputs_plot_HCA, main_outputs_plot_FCA,
                    ncol = 1, nrow = 2)
  
  return(plot)
  
}


# Monte Carlo function ----------------------------------------------------

MC_function <- function(inputs, number_runs){
  
  
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
    
    
    
    
    
    
    inputsPSA[parameter == "c_mort_ind_well", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "c_mort_ind_well", "country param 1"]), 
                                                           as.numeric(inputsPSA[parameter == "c_mort_ind_well", "country param 2"]))
    
    inputsPSA[parameter == "c_mort_small_well", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "c_mort_small_well", "country param 1"]), 
                                                             as.numeric(inputsPSA[parameter == "c_mort_small_well", "country param 2"]))
    
    inputsPSA[parameter == "p_mort_ind_well", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "p_mort_ind_well", "country param 1"]), 
                                                           as.numeric(inputsPSA[parameter == "p_mort_ind_well", "country param 2"]))
    
    inputsPSA[parameter == "p_mort_small_well", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "p_mort_small_well", "country param 1"]), 
                                                             as.numeric(inputsPSA[parameter == "p_mort_small_well", "country param 2"]))

    inputsPSA[parameter == "c_mort_ind_res", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "c_mort_ind_res", "country param 1"]), 
                                                           as.numeric(inputsPSA[parameter == "c_mort_ind_res", "country param 2"]))
    
    inputsPSA[parameter == "c_mort_small_res", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "c_mort_small_res", "country param 1"]), 
                                                             as.numeric(inputsPSA[parameter == "c_mort_small_res", "country param 2"]))
    
    inputsPSA[parameter == "p_mort_ind_res", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "p_mort_ind_res", "country param 1"]), 
                                                           as.numeric(inputsPSA[parameter == "p_mort_ind_res", "country param 2"]))
    
    inputsPSA[parameter == "p_mort_small_res", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "p_mort_small_res", "country param 1"]), 
                                                             as.numeric(inputsPSA[parameter == "p_mort_small_res", "country param 2"]))
    
    inputsPSA[parameter == "c_mort_ind_sus", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "c_mort_ind_sus", "country param 1"]), 
                                                           as.numeric(inputsPSA[parameter == "c_mort_ind_sus", "country param 2"]))
    
    inputsPSA[parameter == "c_mort_small_sus", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "c_mort_small_sus", "country param 1"]), 
                                                             as.numeric(inputsPSA[parameter == "c_mort_small_sus", "country param 2"]))
    
    inputsPSA[parameter == "p_mort_ind_sus", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "p_mort_ind_sus", "country param 1"]), 
                                                           as.numeric(inputsPSA[parameter == "p_mort_ind_sus", "country param 2"]))
    
    inputsPSA[parameter == "p_mort_small_sus", "Value"] <- rbeta(1, as.numeric(inputsPSA[parameter == "p_mort_small_sus", "country param 1"]), 
                                                             as.numeric(inputsPSA[parameter == "p_mort_small_sus", "country param 2"]))
    
    
    
    
    
    inputsPSA[parameter == "seq_sus", "Value"] <- 1.54158796 * rbeta(1, as.numeric(inputsPSA[parameter == "seq_sus", "country param 1"]), 
                                                                     as.numeric(inputsPSA[parameter == "seq_sus", "country param 2"]))
    
    inputsPSA[parameter == "seq_res", "Value"] <- 2.4973725 * rbeta(1, as.numeric(inputsPSA[parameter == "seq_res", "country param 1"]), 
                                                                    as.numeric(inputsPSA[parameter == "seq_res", "country param 2"]))
    
    MCmatrix[i, 1] <- as.numeric(Model(inputsPSA, "HCA")[1,2])
    MCmatrix[i, 2] <- as.numeric(Model(inputsPSA, "HCA")[1,3])
    MCmatrix[i, 3] <- as.numeric(Model(inputsPSA, "FCA")[1,2])
    MCmatrix[i, 4] <- as.numeric(Model(inputsPSA, "FCA")[1,3])
    
  }
  
  
  return(MCmatrix)
  
}


# Monte Carlo summary stats -----------------------------------------------

MC_sumstats_function <- function(MCmatrix){
  
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
  
  return(MonteCarlo_summary_stats)
  
}


# Monte Carlo box and whisker plot ----------------------------------------

MC_boxplot_function <- function(MCmatrix){
  
  number_runs <- nrow(MCmatrix)
  
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
  
  return(boxplot)
  
}
