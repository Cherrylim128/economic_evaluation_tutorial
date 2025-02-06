#'Anti-microbial resistance in LMIC- practical applications in R
#'Ben S Cooper, Cherry Lim
#'
#'1. Runs decision tree model and populates parameter values.
#'2. Estimates and interpret ICER.
#'3. Performs deterministic sensitivity analysis.
#'4. (optional) Performs probability sensitivity analysis.

# set your working directory by replacing "xxx" below with your working file path 
setwd("xxx")
setwd("/Users/Cherry/2025/HDS CDT")
## Example: setwd("/Users/Liz/practicals")
rm(list=ls())

# packages ----
library(dplyr)
library(ggplot2)
library(brms)

# Define the basic population to infer ----
N = 1000 # number of the individual patients in the analysis population
num_bottle = 2 # number of blood bottles collected for microbiology testing
length_empiric = 3 # number of days of empiric antibiotic treatment
length_targeted = 7 # number of days of targetted antibiotic treatment
year_lost =  11  
daly_weight = 0.133 #adults (GBD) 0.133

# Decision tree model ----
source("practical_cea_decisiontree.R")
## define the parameter values in the decision tree model ----
## Active microbiology laboratory service branch
### Probability of each transition path for the chance node
p.bact = 0.40            # Probability of bacterial infection
p.nobact = 1-p.bact        # Probability of non-bacterial infection 

p.bact_culture_pos = 0.8           # Probability of culture positive given bacterial infection 
p.nobact_culture_pos = 0.05        # Probability of culture positive given non-bacterial infection

## pathways of culture positives
p.CefGen_culpos_stepdown = 0.61/10              # Probability of stepping down from ceftriaxone/gentamicin to oral amoxicillin for patient with culture positive result
p.CefGen_culpos_stepup = 0.39                   # Probability of antibiotic stepping up from ceftriaxone/gentamicin for patient with culture positive result
p.CefGen_culpos_nochange = (1-(p.CefGen_culpos_stepdown)-p.CefGen_culpos_stepup)   # Probability of no change in ceftriaxone/gentamicin for patient with culture positive result

## pathways of culture negatives 
p.CefGen_bact_culneg_deter = 0.24/5             # Probability of patient with bacterial infection  whose condition deteriorates under ceftriaxone/gentamicin treatment
p.CefGen_bact_culneg_impro = (1-(0.24/5)-0.16)  # Probability of patient with bacterial infection  whose condition improved under ceftriaxone/gentamicin treatment
p.CefGen_bact_culneg_nocha = 0.16               # Probability of patient with bacterial infection  whose condition remains unchanged under ceftriaxone/gentamicin treatment
p.CefGen_nobact_culneg_deter = 0.07             # Probability of patient with non-bacterial infection  whose condition deteriorates under ceftriaxone/gentamicin treatment
p.CefGen_nobact_culneg_impro = 0.37             # Probability of patient with non-bacterial infection  whose condition improved under ceftriaxone/gentamicin treatment
p.CefGen_nobact_culneg_nocha = 0.56             # Probability of patient with non-bacterial infection  whose condition remains unchanged under ceftriaxone/gentamicin treatment

p.CefGen_culneg_deter_stepup = 0.70             # Probability of stepping up from the current antibiotic treatment when clinical condition deteriorates
p.CefGen_culneg_deter_nochange = 0.30           # Probability of not changing current antibiotic treatment when clinical condition deteriorates
p.CefGen_culneg_impro_stepdown = 0.22           # Probability of stepping down from the current antibiotic treatment when clinical condition improves
p.CefGen_culneg_impro_nochange = 0.78           # Probability of not changing from the current antibiotic treatment when clinical condition improves
p.CefGen_culneg_nocha_stepup = 0.56             # Probability of stepping up from the current antibiotic treatment when clinical condition remains unchanged
p.CefGen_culneg_nocha_nochange = 0.44           # Probability of not changing from the current antibiotic treatment when clinical condition remains unchanged

p.culneg_deter_stepup = 0.012                # Probability of stepping up from the current antibiotic treatment when clinical condition deteriorates
p.culneg_deter_nochange = 0.988                # Probability of not changing current antibiotic treatment when clinical condition deteriorates
p.culneg_impro_stepdown = 0.80                 # Probability of stepping down from the current antibiotic treatment when clinical condition improves
p.culneg_impro_nochange = 0.20                 # Probability of not changing from the current antibiotic treatment when clinical condition improves
p.culneg_nocha_stepup = 0.013                   # Probability of stepping up from the current antibiotic treatment when clinical condition remains unchanged
p.culneg_nocha_nochange = 0.987                 # Probability of not changing from the current antibiotic treatment when clinical condition remains unchanged

## Branches of no culture performed (informed by expert elicitation)
p.CefGen_bact_nocul_deter = 0.24                # Probability of patient with bacterial infection  whose condition deteriorates under ceftriaxone/gentamicin treatment
p.CefGen_bact_nocul_impro = 0.60                # Probability of patient with bacterial infection  whose condition improved under ceftriaxone/gentamicin treatment
p.CefGen_bact_nocul_nocha = 0.16                # Probability of patient with bacterial infection  whose condition remains unchanged under ceftriaxone/gentamicin treatment
p.CefGen_nobact_nocul_deter = 0.07              # Probability of patient with non-bacterial infection  whose condition deteriorates under ceftriaxone/gentamicin treatment
p.CefGen_nobact_nocul_impro = 0.37              # Probability of patient with non-bacterial infection  whose condition improved under ceftriaxone/gentamicin treatment
p.CefGen_nobact_nocul_nocha = 0.56              # Probability of patient with non-bacterial infection  whose condition remains unchanged under ceftriaxone/gentamicin treatment

p.CefGen_nocul_deter_stepup = 0.70             # Probability of stepping up from the current antibiotic treatment when clinical condition deteriorates
p.CefGen_nocul_deter_nochange = 0.30          # Probability of not changing current antibiotic treatment when clinical condition deteriorates
p.CefGen_nocul_impro_stepdown = 0.22           # Probability of stepping down from the current antibiotic treatment when clinical condition improves
p.CefGen_nocul_impro_nochange = 0.78           # Probability of not changing from the current antibiotic treatment when clinical condition improves
p.CefGen_nocul_nocha_stepup = 0.56             # Probability of stepping up from the current antibiotic treatment when clinical condition remains unchanged
p.CefGen_nocul_nocha_nochange = 0.44           # Probability of not changing from the current antibiotic treatment when clinical condition remains unchanged

## Cost of each path
cost_henced_microlab = 300*num_bottle  # cost of maintaining microbiology laboratory
cost_CefGen = 1.76                             # cost of 1-day ceftriaxone/gentamicin
cost_Merope = 35.16                            # cost of 1-day meropenem/vancomycin
cost_oralatb = 0.14                          # cost of stepping down to oral amoxicillin 1g

cost_cultpos = 60*num_bottle                # cost of processing microbiology testing per blood culture positive sample
cost_cultneg = 10*num_bottle                 # cost of processing microbiology testing per blood culture negative sample

cost_careunchan = 158.82                      # additional cost of care per patient per day whose condition remains unchanged
cost_caredetrio = 203.19                      # additional cost of care per patient per day whose condition deteriorates
cost_careimprov = 14.44                       # additional cost of care per patient per day whose condition improved

## LOS of each path 
## Culture positive: in-hospital mortality after step-down  step-up  or no change
los_bact_culpos_stepdown = 7 
los_bact_culpos_stepup = los_bact_culpos_stepdown*1.5 
los_bact_culpos_nochange = los_bact_culpos_stepdown 
## Culture negative in-hospital mortality based on clinical condition
los_bact_culneg_stepup_deter = los_bact_culpos_stepup*2 
los_bact_culneg_nochange_deter = los_bact_culneg_stepup_deter*1.5 
los_bact_culneg_stepdown_impro = los_bact_culpos_stepdown 
los_bact_culneg_nochange_impro = los_bact_culpos_stepdown 
los_bact_culneg_stepup_uncha = los_bact_culpos_nochange 
los_bact_culneg_nochange_uncha = los_bact_culpos_nochange 
## No microbiology testing performed: in-hospital mortality based on clinical condition
los_bact_nocul_stepup_deter = los_bact_culneg_stepup_deter*1.25 
los_bact_nocul_nochange_deter = los_bact_culneg_stepup_deter*2 
los_bact_nocul_stepdown_impro = los_bact_culneg_stepup_deter*1.5 
los_bact_nocul_nochange_impro = los_bact_nocul_stepup_deter 
los_bact_nocul_stepup_uncha = los_bact_nocul_stepup_deter 
los_bact_nocul_nochange_uncha = los_bact_culneg_stepup_deter*1.5 

# Patient with no bacterial infection
## Culture positive: in-hospital los after step-down  step-up  or no change
los_nobact_culpos_stepdown = los_bact_culpos_stepdown/2 
los_nobact_culpos_stepup = los_bact_culpos_stepup/2 
los_nobact_culpos_nochange = los_bact_culpos_nochange/2 
## Culture negative: in-hospital los based on clinical condition (for culture negatives)
los_nobact_culneg_stepup_deter = los_bact_culneg_stepup_deter/2 
los_nobact_culneg_nochange_deter = los_bact_culneg_nochange_deter/2 
los_nobact_culneg_stepdown_impro = los_bact_culneg_stepdown_impro/2 
los_nobact_culneg_nochange_impro = los_bact_culneg_nochange_impro/2 
los_nobact_culneg_stepup_uncha = los_bact_culneg_stepup_uncha/2 
los_nobact_culneg_nochange_uncha = los_bact_culneg_nochange_uncha/2 
## No microbiology testing performed: in-hospital los based on clinical condition
los_nobact_nocul_stepup_deter = los_nobact_culneg_stepup_deter 
los_nobact_nocul_nochange_deter = los_nobact_culneg_nochange_deter 
los_nobact_nocul_stepdown_impro = los_nobact_culneg_stepdown_impro 
los_nobact_nocul_nochange_impro = los_nobact_culneg_nochange_impro 
los_nobact_nocul_stepup_uncha = los_nobact_culneg_stepup_uncha 
los_nobact_nocul_nochange_uncha = los_nobact_culneg_nochange_uncha 

# mortality of each path
# Patient with bacterial infection; culture positive 
## Culture positive: in-hospital mortality after step-down  step-up  or no change
mort_bact_culpos_stepdown = 0.2 
mort_bact_culpos_stepup = mort_bact_culpos_stepdown*1.25 
mort_bact_culpos_nochange = mort_bact_culpos_stepdown*1.25 
## Culture negative in-hospital mortality based on clinical condition
mort_bact_culneg_stepup_deter = mort_bact_culpos_stepdown*1.25 
mort_bact_culneg_nochange_deter = mort_bact_culneg_stepup_deter*1.5 
mort_bact_culneg_stepdown_impro = mort_bact_culpos_stepdown 
mort_bact_culneg_nochange_impro = mort_bact_culpos_stepdown 
mort_bact_culneg_stepup_uncha = mort_bact_culpos_nochange 
mort_bact_culneg_nochange_uncha = mort_bact_culpos_nochange 
## No microbiology testing performed: in-hospital mortality based on clinical condition
mort_bact_nocul_stepup_deter = mort_bact_culneg_stepup_deter*1.25 
mort_bact_nocul_nochange_deter = mort_bact_culneg_stepup_deter*2 
mort_bact_nocul_stepdown_impro = mort_bact_culneg_stepup_deter*1.5 
mort_bact_nocul_nochange_impro = mort_bact_nocul_stepup_deter 
mort_bact_nocul_stepup_uncha = mort_bact_nocul_stepup_deter 
mort_bact_nocul_nochange_uncha = mort_bact_culneg_stepup_deter*1.5 

# Patient with no bacterial infection
## Culture positive: in-hospital mortality after step-down  step-up  or no change
mort_nobact_culpos_stepdown = mort_bact_culpos_stepdown/2 
mort_nobact_culpos_stepup = mort_bact_culpos_stepup/2 
mort_nobact_culpos_nochange = mort_bact_culpos_nochange/2 
## Culture negative: in-hospital mortality based on clinical condition (for culture negatives)
mort_nobact_culneg_stepup_deter = mort_bact_culneg_stepup_deter/2 
mort_nobact_culneg_nochange_deter = mort_bact_culneg_nochange_deter/2 
mort_nobact_culneg_stepdown_impro = mort_bact_culneg_stepdown_impro/2 
mort_nobact_culneg_nochange_impro = mort_bact_culneg_nochange_impro/2 
mort_nobact_culneg_stepup_uncha = mort_bact_culneg_stepup_uncha/2 
mort_nobact_culneg_nochange_uncha = mort_bact_culneg_nochange_uncha/2 
## No microbiology testing performed: in-hospital mortality based on clinical condition
mort_nobact_nocul_stepup_deter = mort_nobact_culneg_stepup_deter 
mort_nobact_nocul_nochange_deter = mort_nobact_culneg_nochange_deter 
mort_nobact_nocul_stepdown_impro = mort_nobact_culneg_stepdown_impro 
mort_nobact_nocul_nochange_impro = mort_nobact_culneg_nochange_impro 
mort_nobact_nocul_stepup_uncha = mort_nobact_culneg_stepup_uncha 
mort_nobact_nocul_nochange_uncha = mort_nobact_culneg_nochange_uncha 

## for the choice of standard surveillance
## Culture positive: in-hospital mortality after step-down  step-up  or no change
yll_bact_culpos_stepdown = mort_bact_culpos_stepdown*year_lost 
yll_bact_culpos_stepup = mort_bact_culpos_stepup*year_lost       
yll_bact_culpos_nochange = mort_bact_culpos_nochange*year_lost    
## Culture negative in-hospital mortality based on clinical condition
yll_bact_culneg_stepup_deter = mort_bact_culneg_stepup_deter*year_lost                 
yll_bact_culneg_nochange_deter = mort_bact_culneg_nochange_deter*year_lost  
yll_bact_culneg_stepdown_impro = mort_bact_culneg_stepdown_impro*year_lost      
yll_bact_culneg_nochange_impro = mort_bact_culneg_nochange_impro*year_lost  
yll_bact_culneg_stepup_uncha = mort_bact_culneg_stepup_uncha*year_lost 
yll_bact_culneg_nochange_uncha = mort_bact_culneg_nochange_uncha*year_lost    
## No microbiology testing performed: in-hospital mortality based on clinical condition
yll_bact_nocul_stepup_deter = mort_bact_nocul_stepup_deter*year_lost                
yll_bact_nocul_nochange_deter = mort_bact_nocul_nochange_deter*year_lost  
yll_bact_nocul_stepdown_impro = mort_bact_nocul_stepdown_impro*year_lost      
yll_bact_nocul_nochange_impro = mort_bact_nocul_nochange_impro*year_lost  
yll_bact_nocul_stepup_uncha = mort_bact_nocul_stepup_uncha*year_lost 
yll_bact_nocul_nochange_uncha = mort_bact_nocul_nochange_uncha*year_lost    

# Patient with no bacterial infection
## Culture positive: in-hospital mortality after step-down  step-up  or no change
yll_nobact_culpos_stepdown = mort_nobact_culpos_stepdown*year_lost 
yll_nobact_culpos_stepup = mort_nobact_culpos_stepup*year_lost       
yll_nobact_culpos_nochange = mort_nobact_culpos_nochange*year_lost    
## Culture negative: in-hospital mortality based on clinical condition (for culture negatives)
yll_nobact_culneg_stepup_deter = mort_nobact_culneg_stepup_deter*year_lost              
yll_nobact_culneg_nochange_deter = mort_nobact_culneg_nochange_deter*year_lost  
yll_nobact_culneg_stepdown_impro = mort_nobact_culneg_stepdown_impro*year_lost      
yll_nobact_culneg_nochange_impro = mort_nobact_culneg_nochange_impro*year_lost  
yll_nobact_culneg_stepup_uncha = mort_nobact_culneg_stepup_uncha*year_lost    
yll_nobact_culneg_nochange_uncha = mort_nobact_culneg_nochange_uncha*year_lost    
## No microbiology testing performed: in-hospital mortality based on clinical condition
yll_nobact_nocul_stepup_deter = mort_nobact_nocul_stepup_deter*year_lost                 
yll_nobact_nocul_nochange_deter = mort_nobact_nocul_nochange_deter*year_lost  
yll_nobact_nocul_stepdown_impro = mort_nobact_nocul_stepdown_impro*year_lost      
yll_nobact_nocul_nochange_impro = mort_nobact_nocul_nochange_impro*year_lost  
yll_nobact_nocul_stepup_uncha = mort_nobact_nocul_stepup_uncha*year_lost    
yll_nobact_nocul_nochange_uncha = mort_nobact_nocul_nochange_uncha*year_lost    

# Patient with bacterial infection
## Culture positive: in-hospital length of hospitalisation after step-down  step-up  or no change
yld_bact_culpos_stepdown = los_bact_culpos_stepdown*daly_weight 
yld_bact_culpos_stepup = los_bact_culpos_stepup*daly_weight         
yld_bact_culpos_nochange = los_bact_culpos_nochange*daly_weight    
## Culture negative in-hospital length of hospitalisation based on clinical condition
yld_bact_culneg_stepup_deter = los_bact_culneg_stepup_deter*daly_weight                 
yld_bact_culneg_nochange_deter = los_bact_culneg_nochange_deter*daly_weight  
yld_bact_culneg_stepdown_impro = los_bact_culneg_stepdown_impro*daly_weight      
yld_bact_culneg_nochange_impro = los_bact_culneg_nochange_impro*daly_weight  
yld_bact_culneg_stepup_uncha = los_bact_culneg_stepup_uncha*daly_weight    
yld_bact_culneg_nochange_uncha = los_bact_culneg_nochange_uncha*daly_weight    
## No microbiology testing performed: in-hospital length of hospitalisation based on clinical condition
yld_bact_nocul_stepup_deter = los_bact_nocul_stepup_deter*daly_weight                 
yld_bact_nocul_nochange_deter = los_bact_nocul_nochange_deter*daly_weight  
yld_bact_nocul_stepdown_impro = los_bact_nocul_stepdown_impro*daly_weight      
yld_bact_nocul_nochange_impro = los_bact_nocul_nochange_impro*daly_weight  
yld_bact_nocul_stepup_uncha = los_bact_nocul_stepup_uncha*daly_weight    
yld_bact_nocul_nochange_uncha = los_bact_nocul_nochange_uncha*daly_weight    

# Patient with no bacterial infection
## Culture positive: in-hospital length of hospitalisation after step-down  step-up  or no change
yld_nobact_culpos_stepdown = los_nobact_culpos_stepdown*daly_weight     
yld_nobact_culpos_stepup = los_nobact_culpos_stepup*daly_weight          
yld_nobact_culpos_nochange = los_nobact_culpos_nochange*daly_weight    
## Culture negative: in-hospital length of hospitalisation based on clinical condition (for culture negatives)
yld_nobact_culneg_stepup_deter = los_nobact_culneg_stepup_deter*daly_weight                 
yld_nobact_culneg_nochange_deter = los_nobact_culneg_nochange_deter*daly_weight  
yld_nobact_culneg_stepdown_impro = los_nobact_culneg_stepdown_impro*daly_weight     
yld_nobact_culneg_nochange_impro = los_nobact_culneg_nochange_impro*daly_weight  
yld_nobact_culneg_stepup_uncha = los_nobact_culneg_stepup_uncha*daly_weight 
yld_nobact_culneg_nochange_uncha = los_nobact_culneg_nochange_uncha*daly_weight    
## No microbiology testing performed: in-hospital length of hospitalisation based on clinical condition
yld_nobact_nocul_stepup_deter = los_nobact_nocul_stepup_deter*daly_weight                 
yld_nobact_nocul_nochange_deter = los_nobact_nocul_nochange_deter*daly_weight  
yld_nobact_nocul_stepdown_impro = los_nobact_nocul_stepdown_impro*daly_weight       
yld_nobact_nocul_nochange_impro = los_nobact_nocul_nochange_impro*daly_weight  
yld_nobact_nocul_stepup_uncha = los_nobact_nocul_stepup_uncha*daly_weight 
yld_nobact_nocul_nochange_uncha = los_nobact_nocul_nochange_uncha*daly_weight

## populate the parameter values in the decision tree model ----
results_deter <- as.data.frame(t(dec_tree(data.frame(
  ## Active microbiology laboratory service branch
  ### Probability of each transition path for the chance node
  p.bact, p.nobact,p.bact_culture_pos,p.nobact_culture_pos,
  p.CefGen_culpos_stepdown,p.CefGen_culpos_stepup,p.CefGen_culpos_nochange,
  ## pathways of culture negatives 
  p.CefGen_bact_culneg_deter,p.CefGen_bact_culneg_impro,p.CefGen_bact_culneg_nocha,              # Probability of patient with bacterial infection, whose condition remains unchanged under ceftriaxone/gentamicin treatment
  p.CefGen_nobact_culneg_deter,p.CefGen_nobact_culneg_impro,p.CefGen_nobact_culneg_nocha,            # Probability of patient with non-bacterial infection, whose condition remains unchanged under ceftriaxone/gentamicin treatment
  p.CefGen_culneg_deter_stepup,p.CefGen_culneg_deter_nochange,p.CefGen_culneg_impro_stepdown,          # Probability of stepping down from the current antibiotic treatment when clinical condition improves
  p.CefGen_culneg_impro_nochange,p.CefGen_culneg_nocha_stepup,p.CefGen_culneg_nocha_nochange,          # Probability of not changing from the current antibiotic treatment when clinical condition remains unchanged
  p.culneg_deter_stepup,p.culneg_deter_nochange,p.culneg_impro_stepdown,                # Probability of stepping down from the current antibiotic treatment when clinical condition improves
  p.culneg_impro_nochange,p.culneg_nocha_stepup,p.culneg_nocha_nochange,                # Probability of not changing from the current antibiotic treatment when clinical condition remains unchanged
  
  ## Branches of no culture performed (informed by expert elicitation)
  p.CefGen_bact_nocul_deter,p.CefGen_bact_nocul_impro,p.CefGen_bact_nocul_nocha,               # Probability of patient with bacterial infection, whose condition remains unchanged under ceftriaxone/gentamicin treatment
  p.CefGen_nobact_nocul_deter,p.CefGen_nobact_nocul_impro,p.CefGen_nobact_nocul_nocha,             # Probability of patient with non-bacterial infection, whose condition remains unchanged under ceftriaxone/gentamicin treatment
  
  p.CefGen_nocul_deter_stepup, p.CefGen_nocul_deter_nochange, p.CefGen_nocul_impro_stepdown,          # Probability of stepping down from the current antibiotic treatment when clinical condition improves
  p.CefGen_nocul_impro_nochange,p.CefGen_nocul_nocha_stepup,p.CefGen_nocul_nocha_nochange,          # Probability of not changing from the current antibiotic treatment when clinical condition remains unchanged
  
  ## Cost of each path
  cost_henced_microlab,cost_CefGen,cost_Merope,cost_oralatb,                         # cost of stepping down to oral amoxicillin 1g
  cost_cultpos, cost_cultneg,                # cost of processing microbiology testing per blood culture negative sample
  cost_careunchan,cost_caredetrio, cost_careimprov,                      # additional cost of care per patient per day whose condition improved
  
  ## LOS of each path 
  ## Culture positive: in-hospital mortality after step-down, step-up, or no change
  los_bact_culpos_stepdown,los_bact_culpos_stepup,los_bact_culpos_nochange,
  ## Culture negative in-hospital mortality based on clinical condition
  los_bact_culneg_stepup_deter,los_bact_culneg_nochange_deter,
  los_bact_culneg_stepdown_impro,los_bact_culneg_nochange_impro,
  los_bact_culneg_stepup_uncha,los_bact_culneg_nochange_uncha,
  ## No microbiology testing performed: in-hospital mortality based on clinical condition
  los_bact_nocul_stepup_deter,los_bact_nocul_nochange_deter,
  los_bact_nocul_stepdown_impro,los_bact_nocul_nochange_impro,
  los_bact_nocul_stepup_uncha,los_bact_nocul_nochange_uncha,
  
  # Patient with no bacterial infection
  ## Culture positive: in-hospital los after step-down, step-up, or no change
  los_nobact_culpos_stepdown,los_nobact_culpos_stepup,los_nobact_culpos_nochange,
  ## Culture negative: in-hospital los based on clinical condition (for culture negatives)
  los_nobact_culneg_stepup_deter,los_nobact_culneg_nochange_deter,
  los_nobact_culneg_stepdown_impro,los_nobact_culneg_nochange_impro,
  los_nobact_culneg_stepup_uncha,los_nobact_culneg_nochange_uncha,
  ## No microbiology testing performed: in-hospital los based on clinical condition
  los_nobact_nocul_stepup_deter,los_nobact_nocul_nochange_deter,
  los_nobact_nocul_stepdown_impro,los_nobact_nocul_nochange_impro,
  los_nobact_nocul_stepup_uncha,los_nobact_nocul_nochange_uncha,
  
  # mortality of each path
  # Patient with bacterial infection; culture positive 
  ## Culture positive: in-hospital mortality after step-down, step-up, or no change
  mort_bact_culpos_stepdown,mort_bact_culpos_stepup,mort_bact_culpos_nochange,
  ## Culture negative in-hospital mortality based on clinical condition
  mort_bact_culneg_stepup_deter,mort_bact_culneg_nochange_deter,
  mort_bact_culneg_stepdown_impro,mort_bact_culneg_nochange_impro,
  mort_bact_culneg_stepup_uncha,mort_bact_culneg_nochange_uncha,
  ## No microbiology testing performed: in-hospital mortality based on clinical condition
  mort_bact_nocul_stepup_deter,mort_bact_nocul_nochange_deter,
  mort_bact_nocul_stepdown_impro,mort_bact_nocul_nochange_impro,
  mort_bact_nocul_stepup_uncha,mort_bact_nocul_nochange_uncha,
  
  # Patient with no bacterial infection
  ## Culture positive: in-hospital mortality after step-down, step-up, or no change
  mort_nobact_culpos_stepdown,mort_nobact_culpos_stepup,mort_nobact_culpos_nochange,
  ## Culture negative: in-hospital mortality based on clinical condition (for culture negatives)
  mort_nobact_culneg_stepup_deter,mort_nobact_culneg_nochange_deter,
  mort_nobact_culneg_stepdown_impro, mort_nobact_culneg_nochange_impro,
  mort_nobact_culneg_stepup_uncha,mort_nobact_culneg_nochange_uncha,
  ## No microbiology testing performed: in-hospital mortality based on clinical condition
  mort_nobact_nocul_stepup_deter,mort_nobact_nocul_nochange_deter,
  mort_nobact_nocul_stepdown_impro,mort_nobact_nocul_nochange_impro,
  mort_nobact_nocul_stepup_uncha,mort_nobact_nocul_nochange_uncha,
  
  ## for the choice of standard surveillance
  ## Culture positive: in-hospital mortality after step-down, step-up, or no change
  yll_bact_culpos_stepdown,yll_bact_culpos_stepup, yll_bact_culpos_nochange,   
  ## Culture negative in-hospital mortality based on clinical condition
  yll_bact_culneg_stepup_deter,yll_bact_culneg_nochange_deter, 
  yll_bact_culneg_stepdown_impro,yll_bact_culneg_nochange_impro, 
  yll_bact_culneg_stepup_uncha,yll_bact_culneg_nochange_uncha,   
  ## No microbiology testing performed: in-hospital mortality based on clinical condition
  yll_bact_nocul_stepup_deter, yll_bact_nocul_nochange_deter, 
  yll_bact_nocul_stepdown_impro, yll_bact_nocul_nochange_impro, 
  yll_bact_nocul_stepup_uncha,yll_bact_nocul_nochange_uncha,   
  
  # Patient with no bacterial infection
  ## Culture positive: in-hospital mortality after step-down, step-up, or no change
  yll_nobact_culpos_stepdown,yll_nobact_culpos_stepup,yll_nobact_culpos_nochange,   
  ## Culture negative: in-hospital mortality based on clinical condition (for culture negatives)
  yll_nobact_culneg_stepup_deter,yll_nobact_culneg_nochange_deter, 
  yll_nobact_culneg_stepdown_impro, yll_nobact_culneg_nochange_impro, 
  yll_nobact_culneg_stepup_uncha,yll_nobact_culneg_nochange_uncha,   
  ## No microbiology testing performed: in-hospital mortality based on clinical condition
  yll_nobact_nocul_stepup_deter, yll_nobact_nocul_nochange_deter, 
  yll_nobact_nocul_stepdown_impro,yll_nobact_nocul_nochange_impro, 
  yll_nobact_nocul_stepup_uncha,yll_nobact_nocul_nochange_uncha,   
  
  # Patient with bacterial infection
  ## Culture positive: in-hospital length of hospitalisation after step-down, step-up, or no change
  yld_bact_culpos_stepdown,yld_bact_culpos_stepup,yld_bact_culpos_nochange,   
  ## Culture negative in-hospital length of hospitalisation based on clinical condition
  yld_bact_culneg_stepup_deter,yld_bact_culneg_nochange_deter, 
  yld_bact_culneg_stepdown_impro,yld_bact_culneg_nochange_impro, 
  yld_bact_culneg_stepup_uncha,yld_bact_culneg_nochange_uncha,   
  ## No microbiology testing performed: in-hospital length of hospitalisation based on clinical condition
  yld_bact_nocul_stepup_deter,yld_bact_nocul_nochange_deter, 
  yld_bact_nocul_stepdown_impro,yld_bact_nocul_nochange_impro, 
  yld_bact_nocul_stepup_uncha,yld_bact_nocul_nochange_uncha,   
  
  # Patient with no bacterial infection
  ## Culture positive: in-hospital length of hospitalisation after step-down, step-up, or no change
  yld_nobact_culpos_stepdown,yld_nobact_culpos_stepup,yld_nobact_culpos_nochange,   
  ## Culture negative: in-hospital length of hospitalisation based on clinical condition (for culture negatives)
  yld_nobact_culneg_stepup_deter, yld_nobact_culneg_nochange_deter, 
  yld_nobact_culneg_stepdown_impro, yld_nobact_culneg_nochange_impro, 
  yld_nobact_culneg_stepup_uncha,yld_nobact_culneg_nochange_uncha,   
  ## No microbiology testing performed: in-hospital length of hospitalisation based on clinical condition
  yld_nobact_nocul_stepup_deter, yld_nobact_nocul_nochange_deter, 
  yld_nobact_nocul_stepdown_impro, yld_nobact_nocul_nochange_impro, 
  yld_nobact_nocul_stepup_uncha,yld_nobact_nocul_nochange_uncha))))
  
colnames(results_deter) <- c("diff_im", "diff_cost", "diff_daly", "icer",
                             "mort.noactivelab_CefGen", "mort.activelab_CefGen",
                             "cost.noactivelab_CefGen", "cost.activelab_CefGen",
                             "daly.noactivelab_CefGen", "daly.activelab_CefGen")
# Deterministic sensitivity analysis ----
## write result table for Exercise 3.1
dta.deterministic.sen1.output <- data.frame(
  label = c("high", "high", "base", "low", "low"),
  value = c(65, 50, 40, 30, 15),
  diff_cost = c(310057.60, 355346.50, 385539.10, 415731.70, 461020.60),
  diff_mort = c(64, 48, 38, 28, 13),
  diff_daly = c(2427, 1851, 1468, 1084, 508)
)
dta.deterministic.sen1.output$icer = dta.deterministic.sen1.output$diff_cost/dta.deterministic.sen1.output$diff_daly

dta.deterministic.sen1 <- data.frame(
  type = c("high", "high", "low", "low"),
  value = c(65, 50, 30, 15), 
  xmin = c(0,1,2,3),
  xmax = c(0.95, 1.85, 2.95, 3.95),
  ymin = c(127.7534, 191.9754, 262.6288, 262.6288),
  ymax = c(262.6288, 262.6288, 383.5163, 907.5209)
)

## Plot for results from the fist deterministic sensitivity analysis
### on different values of proportion 
fun_sen_plot <- function(data, text_size){
  ggplot() + 
    geom_rect(data = data, 
              aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, fill=type)) +
    theme_bw() + 
    theme(axis.title.y=element_text(size=text_size), 
          #axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          #panel.grid.major.y = element_blank(),
          #panel.grid.minor.y = element_blank(),
          axis.text=element_text(size=text_size),
          legend.position = 'bottom',
          legend.title = element_blank()) +
    coord_flip()
}

plot.deterministic.sen1 <- fun_sen_plot(dta.deterministic.sen1, 12) + 
  labs(x = "Proportion of bacterial infection among the suspected BSI", y = "Incremental cost-effectiveness ratio", size = 12) +
  scale_x_continuous(breaks = c(0.5, 1.5, 2.5, 3.5), 
                   labels = c("65%", "50%", "30%", "15%"))

## write result table for Exercise 3.2
dta.deterministic.sen2.output <- data.frame(
  label = c("high", "high", "high", "base", "low", "low"),
  value = c(20, 13, 7, 6, 5, 2),
  diff_cost = c(402171.70, 397539.70, 387958.30, 385539.10, 382309.60, 376700.50),
  diff_mort = c(40, 39, 38, 38, 38, 37),
  diff_daly = c(1461, 1457, 1465, 1468, 1472, 1477)
)
dta.deterministic.sen2.output$icer = dta.deterministic.sen2.output$diff_cost/dta.deterministic.sen2.output$diff_daly

dta.deterministic.sen2 <- data.frame(
  type = c("high", "high", "high", "low", "low"),
  value = c(20, 13, 7, 5, 2), 
  xmin = c(0,1,2,3,4),
  xmax = c(0.95, 1.85, 2.95, 3.95, 4.95),
  ymin = c(275.2715, 272.8481, 264.8180, 262.6288, 262.6288),
  ymax = c(262.6288, 262.6288, 262.6288, 259.7212, 255.0443)
)

plot.deterministic.sen2 <- fun_sen_plot(dta.deterministic.sen2, 12) + 
  labs(x = "Proportion of stepping-down given culture positive", y = "Incremental cost-effectiveness ratio", size = 12) +
  scale_x_continuous(breaks = c(0.5, 1.5, 2.5, 3.5, 4.5), 
                     labels = c("20%", "13%", "7%", "5%", "2%"))

# Probabilistic sensitivity analysis ----
results_prob <- data.frame()
niter = 1000
for (i in 1:niter){
  set.seed(2000+i)
  source("practical_cea_probabilisticsen.R")
  print(i)
  temp <- as.data.frame(t(dec_tree(input_prob)))
  colnames(temp) <- c("diff_im", "diff_cost", "diff_daly", "icer",
                      "mort.noactivelab_CefGen", "mort.activelab_CefGen",
                      "cost.noactivelab_CefGen", "cost.activelab_CefGen",
                      "daly.noactivelab_CefGen", "daly.activelab_CefGen")
  results_prob <- rbind(results_prob, temp)
}
mean(results_prob$diff_cost)
mean(results_prob$diff_daly)
mean(results_prob$diff_im)

## Plot CEA
ggplot(results_prob, 
       aes(x=diff_daly, y=diff_cost)) +
  labs(
    y="Incremental costs of maintaining active microbiology service system ($)",
    x="DALY averted per 1,000 patients") +
  xlim(-3000, 3000) +
  ylim(-800000, 800000) + 
  geom_point(size=1) +
  theme_minimal() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=16)) +
  #scale_y_continuous(labels = comma, limits = c(round(min(dta_ceplane_CefGen2$ic), 0)*1.25, 100000)) +
  geom_hline(yintercept=0, colour = "darkgrey") +
  geom_vline(xintercept=0, colour = "darkgrey") +
  geom_abline(slope = 500, linetype=2, colour="darkgrey") +
  annotate("point", x = 1468, y = 385539.10, colour = "red", size=3)
