#'Anti-microbial resistance in LMIC- practical applications in R
#'Ben S Cooper, Cherry Lim
#'
#'4. (optional) Performs probability sensitivity analysis.

## Probabilistic parameters ####
p.infecttype = rbeta(1, 40, 60)
### Microbiology testing service
p.CefGen_culpos_atb = rdirichlet(1, c(6.1,39,54.9))
p.CefGen_bact_culneg_dir = rdirichlet(1, c(4.8,79.2,16))
p.CefGen_nobact_culneg_dir = rdirichlet(1, c(7,37,56))
p.CefGen_culneg_deter_stepup_nochange = rbeta(1, 70, 30)
p.CefGen_culneg_impro_stepdown_nochange = rbeta(1, 22, 78)
p.CefGen_culneg_nocha_stepup_nochange = rbeta(1, 56, 44)

p.culneg_deter_stepup_nochange = rbeta(1, 1.2, 98.8)
p.culneg_impro_stepdown_nochange = rbeta(1, 80, 20)
p.culneg_nocha_stepup_nochange = rbeta(1, 1.3, 98.7)

### No microbiological testing service
p.CefGen_bact_nocul_dir = rdirichlet(1, c(24,60,16))
p.CefGen_nobact_nocul_dir = rdirichlet(1, c(7, 37, 56))
p.CefGen_nocul_deter_stepup_nochange = rbeta(1, 70, 30)
p.CefGen_nocul_impro_stepdown_nochange = rbeta(1, 22, 78)
p.CefGen_nocul_nocha_stepup_nochange = rbeta(1, 56, 44)

input_prob <- data.frame(
  
  ## For the choice of enhanced surveillance
  # Probability of each transition path for the chance node ----
  p.bact = p.infecttype,            # Probability of bacterial infection
  p.nobact = 1-p.infecttype,              # Probability of non-bacterial infection 
  
  p.bact_culture_pos = rbeta(1, 80, 20),
  p.nobact_culture_pos = rbeta(1, 5, 95),
  
  ## Branches of culture positives
  p.CefGen_culpos_stepdown = p.CefGen_culpos_atb[1],    # Probability of stepping down from ceftriaxone/gentamicin to oral amoxicillin for patient with culture positive result
  p.CefGen_culpos_stepup = p.CefGen_culpos_atb[2],      # Probability of antibiotic stepping up from ceftriaxone/gentamicin for patient with culture positive result
  p.CefGen_culpos_nochange = p.CefGen_culpos_atb[3],    # Probability of no change in ceftriaxone/gentamicin for patient with culture positive result

  ## Branches of culture negatives 
  p.CefGen_bact_culneg_deter = p.CefGen_bact_culneg_dir[1],            # Probability of patient with bacterial infection, whose condition deteriorates under ceftriaxone/gentamicin treatment
  p.CefGen_bact_culneg_impro = p.CefGen_bact_culneg_dir[2],            # Probability of patient with bacterial infection, whose condition improved under ceftriaxone/gentamicin treatment
  p.CefGen_bact_culneg_nocha = p.CefGen_bact_culneg_dir[3],            # Probability of patient with bacterial infection, whose condition remains unchanged under ceftriaxone/gentamicin treatment
  p.CefGen_nobact_culneg_deter = p.CefGen_nobact_culneg_dir[1],          # Probability of patient with non-bacterial infection, whose condition deteriorates under ceftriaxone/gentamicin treatment
  p.CefGen_nobact_culneg_impro = p.CefGen_nobact_culneg_dir[2],          # Probability of patient with non-bacterial infection, whose condition improved under ceftriaxone/gentamicin treatment
  p.CefGen_nobact_culneg_nocha = p.CefGen_nobact_culneg_dir[3],          # Probability of patient with non-bacterial infection, whose condition remains unchanged under ceftriaxone/gentamicin treatment
  
  p.CefGen_culneg_deter_stepup = p.CefGen_culneg_deter_stepup_nochange,          # Probability of stepping up from the current antibiotic treatment when clincial condition deteriorates
  p.CefGen_culneg_deter_nochange = 1-p.CefGen_culneg_deter_stepup_nochange,        # Probability of not changing current antibiotic treatment when clincial condition deteriorates
  p.CefGen_culneg_impro_stepdown = p.CefGen_culneg_impro_stepdown_nochange,         # Probability of stepping down from the current antibiotic treatment when clincial condition improves
  p.CefGen_culneg_impro_nochange = 1-p.CefGen_culneg_impro_stepdown_nochange,         # Probability of not changing from the current antibiotic treatment when clincial condition improves
  p.CefGen_culneg_nocha_stepup = p.CefGen_culneg_nocha_stepup_nochange,           # Probability of stepping up from the current antibiotic treatment when clincial condition remains unchanged
  p.CefGen_culneg_nocha_nochange = 1-p.CefGen_culneg_nocha_stepup_nochange,         # Probability of not changing from the current antibiotic treatment when clincial condition remains unchanged
  
  p.culneg_deter_stepup = p.culneg_deter_stepup_nochange,          # Probability of stepping up from the current antibiotic treatment when clincial condition deteriorates
  p.culneg_deter_nochange = 1-p.culneg_deter_stepup_nochange,        # Probability of not changing current antibiotic treatment when clincial condition deteriorates
  p.culneg_impro_stepdown = p.culneg_impro_stepdown_nochange,         # Probability of stepping down from the current antibiotic treatment when clincial condition improves
  p.culneg_impro_nochange = 1-p.culneg_impro_stepdown_nochange,         # Probability of not changing from the current antibiotic treatment when clincial condition improves
  p.culneg_nocha_stepup = p.culneg_nocha_stepup_nochange,           # Probability of stepping up from the current antibiotic treatment when clincial condition remains unchanged
  p.culneg_nocha_nochange = 1-p.culneg_nocha_stepup_nochange,         # Probability of not changing from the current antibiotic treatment when clincial condition remains unchanged
  
  ## Branches of no culture performed
  p.CefGen_bact_nocul_deter = p.CefGen_bact_nocul_dir[1],        # Probability of patient with bacterial infection, whose condition deteriorates under ceftriaxone/gentamicin treatment
  p.CefGen_bact_nocul_impro = p.CefGen_bact_nocul_dir[2],        # Probability of patient with bacterial infection, whose condition improved under ceftriaxone/gentamicin treatment
  p.CefGen_bact_nocul_nocha = p.CefGen_bact_nocul_dir[3],        # Probability of patient with bacterial infection, whose condition remains unchanged under ceftriaxone/gentamicin treatment
  p.CefGen_nobact_nocul_deter = p.CefGen_nobact_nocul_dir[1],      # Probability of patient with non-bacterial infection, whose condition deteriorates under ceftriaxone/gentamicin treatment
  p.CefGen_nobact_nocul_impro = p.CefGen_nobact_nocul_dir[2],      # Probability of patient with non-bacterial infection, whose condition improved under ceftriaxone/gentamicin treatment
  p.CefGen_nobact_nocul_nocha = p.CefGen_nobact_nocul_dir[3],      # Probability of patient with non-bacterial infection, whose condition remains unchanged under ceftriaxone/gentamicin treatment

  p.CefGen_nocul_deter_stepup = p.CefGen_nocul_deter_stepup_nochange,          # Probability of stepping up from the current antibiotic treatment when clinical condition deteriorates
  p.CefGen_nocul_deter_nochange = 1-p.CefGen_nocul_deter_stepup_nochange,        # Probability of not changing current antibiotic treatment when clinical condition deteriorates
  p.CefGen_nocul_impro_stepdown = p.CefGen_nocul_impro_stepdown_nochange,         # Probability of stepping down from the current antibiotic treatment when clinical condition improves
  p.CefGen_nocul_impro_nochange = 1-p.CefGen_nocul_impro_stepdown_nochange,         # Probability of not changing from the current antibiotic treatment when clinical condition improves
  p.CefGen_nocul_nocha_stepup = p.CefGen_nocul_nocha_stepup_nochange,           # Probability of stepping up from the current antibiotic treatment when clinical condition remains unchanged
  p.CefGen_nocul_nocha_nochange = 1-p.CefGen_nocul_nocha_stepup_nochange,         # Probability of not changing from the current antibiotic treatment when clinical condition remains unchanged
  
  ## Cost of each path
  cost_henced_microlab = 300*num_bottle,  # cost of maintaining microbiology laboratory
  cost_CefGen = 1.76,                             # cost of 1-day ceftriaxone/gentamicin
  cost_Merope = 35.16,                            # cost of 1-day meropenem/vancomycin
  cost_oralatb = 0.14,                          # cost of stepping down to oral amoxicillin 1g
  
  cost_cultpos = 60*num_bottle,                # cost of processing microbiology testing per blood culture positive sample
  cost_cultneg = 10*num_bottle,                 # cost of processing microbiology testing per blood culture negative sample
  
  cost_careunchan = rlnorm(1, meanlog = log(158.82), sdlog = 0.25),          # cost of care per patient whose condition remains unchanged
  cost_caredetrio = rlnorm(1, meanlog = log(203.19), sdlog = 0.25),          # cost of care per patient whose condition deteriorates
  cost_careimprov = rlnorm(1, meanlog = log(14.44), sdlog = 0.25),           # cost of care per patient whose condition improved
  
  # LOS of each path ----
  # Patients with bacterial infection
  # using log normal distribution
  ## Culture positive: in-hospital mortality after step-down, step-up, or no change
  los_bact_culpos_stepdown <- rlnorm(1, meanlog = log(7), sdlog = 0.1),
  los_bact_culpos_stepup <- los_bact_culpos_stepdown*1.5,
  los_bact_culpos_nochange <- los_bact_culpos_stepdown,
  ## Culture negative in-hospital mortality based on clinical condition
  los_bact_culneg_stepup_deter <- los_bact_culpos_stepup*2,
  los_bact_culneg_nochange_deter <- los_bact_culneg_stepup_deter*1.5,
  los_bact_culneg_stepdown_impro <- los_bact_culpos_stepdown, # los_bact_culneg_stepup_deter/3,
  los_bact_culneg_nochange_impro <- los_bact_culpos_stepdown, # los_bact_culneg_stepup_deter/3,
  los_bact_culneg_stepup_uncha <- los_bact_culpos_nochange, # los_bact_culneg_stepup_deter/3,
  los_bact_culneg_nochange_uncha <- los_bact_culpos_nochange, # los_bact_culneg_stepup_deter/2,
  ## No microbiology testing performed: in-hospital mortality based on clinical condition
  los_bact_nocul_stepup_deter <- los_bact_culneg_stepup_deter*1.25,
  los_bact_nocul_nochange_deter <- los_bact_culneg_stepup_deter*2,
  los_bact_nocul_stepdown_impro <- los_bact_culneg_stepup_deter*1.5,
  los_bact_nocul_nochange_impro <- los_bact_nocul_stepup_deter,
  los_bact_nocul_stepup_uncha <- los_bact_nocul_stepup_deter,
  los_bact_nocul_nochange_uncha <- los_bact_culneg_stepup_deter*1.5,
  
  # Patient with no bacterial infection
  ## Culture positive: in-hospital losality after step-down, step-up, or no change
  los_nobact_culpos_stepdown <- los_bact_culpos_stepdown/2,
  los_nobact_culpos_stepup <- los_bact_culpos_stepup/2,
  los_nobact_culpos_nochange <- los_bact_culpos_nochange/2,
  ## Culture negative: in-hospital losality based on clinical condition (for culture negatives)
  los_nobact_culneg_stepup_deter <- los_bact_culneg_stepup_deter/2,
  los_nobact_culneg_nochange_deter <- los_bact_culneg_nochange_deter/2,
  los_nobact_culneg_stepdown_impro <- los_bact_culneg_stepdown_impro/2,
  los_nobact_culneg_nochange_impro <- los_bact_culneg_nochange_impro/2,
  los_nobact_culneg_stepup_uncha <- los_bact_culneg_stepup_uncha/2,
  los_nobact_culneg_nochange_uncha <- los_bact_culneg_nochange_uncha/2,
  ## No microbiology testing performed: in-hospital losality based on clinical condition
  los_nobact_nocul_stepup_deter <- los_nobact_culneg_stepup_deter,
  los_nobact_nocul_nochange_deter <- los_nobact_culneg_nochange_deter,
  los_nobact_nocul_stepdown_impro <- los_nobact_culneg_stepdown_impro,
  los_nobact_nocul_nochange_impro <- los_nobact_culneg_nochange_impro,
  los_nobact_nocul_stepup_uncha <- los_nobact_culneg_stepup_uncha,
  los_nobact_nocul_nochange_uncha <- los_nobact_culneg_nochange_uncha,
  
  # Mortality of each path ----
  # Patient with bacterial infection
  ## in-hospital mortality after step-down, step-up, or no change
  ## Culture positive: in-hospital mortality after step-down, step-up, or no change
  mort_bact_culpos_stepdown <- rbeta(1, 20, 80),
  mort_bact_culpos_stepup <- mort_bact_culpos_stepdown*1.25,
  mort_bact_culpos_nochange <- mort_bact_culpos_stepdown*1.25,
  ## Culture negative in-hospital mortality based on clinical condition
  mort_bact_culneg_stepup_deter <- mort_bact_culpos_stepdown*1.25,
  mort_bact_culneg_nochange_deter <- mort_bact_culneg_stepup_deter*1.5,
  mort_bact_culneg_stepdown_impro <- mort_bact_culpos_stepdown, # mort_bact_culneg_stepup_deter/3,
  mort_bact_culneg_nochange_impro <- mort_bact_culpos_stepdown, # mort_bact_culneg_stepup_deter/3,
  mort_bact_culneg_stepup_uncha <- mort_bact_culpos_nochange, # mort_bact_culneg_stepup_deter/3,
  mort_bact_culneg_nochange_uncha <- mort_bact_culpos_nochange, #mort_bact_culneg_stepup_deter/2,
  ## No microbiology testing performed: in-hospital mortality based on clinical condition
  mort_bact_nocul_stepup_deter <- mort_bact_culneg_stepup_deter*1.25,
  mort_bact_nocul_nochange_deter <- mort_bact_culneg_stepup_deter*2,
  mort_bact_nocul_stepdown_impro <- mort_bact_culneg_stepup_deter*1.5,
  mort_bact_nocul_nochange_impro <- mort_bact_nocul_stepup_deter,
  mort_bact_nocul_stepup_uncha <- mort_bact_nocul_stepup_deter,
  mort_bact_nocul_nochange_uncha <- mort_bact_culneg_stepup_deter*1.5,
  
  # Patient with no bacterial infection
  ## Culture positive: in-hospital mortality after step-down, step-up, or no change
  mort_nobact_culpos_stepdown <- mort_bact_culpos_stepdown/2,
  mort_nobact_culpos_stepup <- mort_bact_culpos_stepup/2,
  mort_nobact_culpos_nochange <- mort_bact_culpos_nochange/2,
  ## Culture negative: in-hospital mortality based on clinical condition (for culture negatives)
  mort_nobact_culneg_stepup_deter <- mort_bact_culneg_stepup_deter/2,
  mort_nobact_culneg_nochange_deter <- mort_bact_culneg_nochange_deter/2,
  mort_nobact_culneg_stepdown_impro <- mort_bact_culneg_stepdown_impro/2,
  mort_nobact_culneg_nochange_impro <- mort_bact_culneg_nochange_impro/2,
  mort_nobact_culneg_stepup_uncha <- mort_bact_culneg_stepup_uncha/2,
  mort_nobact_culneg_nochange_uncha <- mort_bact_culneg_nochange_uncha/2,
  ## No microbiology testing performed: in-hospital mortality based on clinical condition
  mort_nobact_nocul_stepup_deter <- mort_nobact_culneg_stepup_deter,
  mort_nobact_nocul_nochange_deter <- mort_nobact_culneg_nochange_deter,
  mort_nobact_nocul_stepdown_impro <- mort_nobact_culneg_stepdown_impro,
  mort_nobact_nocul_nochange_impro <- mort_nobact_culneg_nochange_impro,
  mort_nobact_nocul_stepup_uncha <- mort_nobact_culneg_stepup_uncha,
  mort_nobact_nocul_nochange_uncha <- mort_nobact_culneg_nochange_uncha,
  
  # DALY of each path ----
  # DALY = YLL + YLD
  # YLL is calculated as the number of deaths (n) x the standard life expectancy at age of death (L1). 
  # This measures the reduction in life expectancy.
  # Source: https://www.who.int/data/gho/data/indicators/indicator-details/GHO/gho-ghe-life-tables-by-country
  # Patient with bacterial infection
  ## Culture positive: in-hospital mortality after step-down, step-up, or no change
  yll_bact_culpos_stepdown <- mort_bact_culpos_stepdown*year_lost,       # expectation of life at age 55-59 years
  yll_bact_culpos_stepup <- mort_bact_culpos_stepup*year_lost,         
  yll_bact_culpos_nochange <- mort_bact_culpos_nochange*year_lost,   
  ## Culture negative in-hospital mortality based on clinical condition
  yll_bact_culneg_stepup_deter <- mort_bact_culneg_stepup_deter*year_lost,                
  yll_bact_culneg_nochange_deter <- mort_bact_culneg_nochange_deter*year_lost, 
  yll_bact_culneg_stepdown_impro <- mort_bact_culneg_stepdown_impro*year_lost,     
  yll_bact_culneg_nochange_impro <- mort_bact_culneg_nochange_impro*year_lost, 
  yll_bact_culneg_stepup_uncha <- mort_bact_culneg_stepup_uncha*year_lost,   
  yll_bact_culneg_nochange_uncha <- mort_bact_culneg_nochange_uncha*year_lost,   
  ## No microbiology testing performed: in-hospital mortality based on clinical condition
  yll_bact_nocul_stepup_deter <- mort_bact_nocul_stepup_deter*year_lost,                
  yll_bact_nocul_nochange_deter <- mort_bact_nocul_nochange_deter*year_lost, 
  yll_bact_nocul_stepdown_impro <- mort_bact_nocul_stepdown_impro*year_lost,     
  yll_bact_nocul_nochange_impro <- mort_bact_nocul_nochange_impro*year_lost, 
  yll_bact_nocul_stepup_uncha <- mort_bact_nocul_stepup_uncha*year_lost,   
  yll_bact_nocul_nochange_uncha <- mort_bact_nocul_nochange_uncha*year_lost,   
  
  # Patient with no bacterial infection
  ## Culture positive: in-hospital mortality after step-down, step-up, or no change
  yll_nobact_culpos_stepdown <- mort_nobact_culpos_stepdown*year_lost,       # expectation of life at age 55-59 years
  yll_nobact_culpos_stepup <- mort_nobact_culpos_stepup*year_lost,         
  yll_nobact_culpos_nochange <- mort_nobact_culpos_nochange*year_lost,   
  ## Culture negative: in-hospital mortality based on clinical condition (for culture negatives)
  yll_nobact_culneg_stepup_deter <- mort_nobact_culneg_stepup_deter*year_lost,                
  yll_nobact_culneg_nochange_deter <- mort_nobact_culneg_nochange_deter*year_lost, 
  yll_nobact_culneg_stepdown_impro <- mort_nobact_culneg_stepdown_impro*year_lost,     
  yll_nobact_culneg_nochange_impro <- mort_nobact_culneg_nochange_impro*year_lost, 
  yll_nobact_culneg_stepup_uncha <- mort_nobact_culneg_stepup_uncha*year_lost,   
  yll_nobact_culneg_nochange_uncha <- mort_nobact_culneg_nochange_uncha*year_lost,   
  ## No microbiology testing performed: in-hospital mortality based on clinical condition
  yll_nobact_nocul_stepup_deter <- mort_nobact_nocul_stepup_deter*year_lost,                
  yll_nobact_nocul_nochange_deter <- mort_nobact_nocul_nochange_deter*year_lost, 
  yll_nobact_nocul_stepdown_impro <- mort_nobact_nocul_stepdown_impro*year_lost,     
  yll_nobact_nocul_nochange_impro <- mort_nobact_nocul_nochange_impro*year_lost, 
  yll_nobact_nocul_stepup_uncha <- mort_nobact_nocul_stepup_uncha*year_lost,   
  yll_nobact_nocul_nochange_uncha <- mort_nobact_nocul_nochange_uncha*year_lost,
  
  # Patient with bacterial infection
  ## Culture positive: in-hospital length of hospitalisation after step-down, step-up, or no change
  yld_bact_culpos_stepdown <- los_bact_culpos_stepdown*daly_weight,       # expectation of life at age 55-59 years
  yld_bact_culpos_stepup <- los_bact_culpos_stepup*daly_weight,         
  yld_bact_culpos_nochange <- los_bact_culpos_nochange*daly_weight,   
  ## Culture negative in-hospital length of hospitalisation based on clinical condition
  yld_bact_culneg_stepup_deter <- los_bact_culneg_stepup_deter*daly_weight,                
  yld_bact_culneg_nochange_deter <- los_bact_culneg_nochange_deter*daly_weight, 
  yld_bact_culneg_stepdown_impro <- los_bact_culneg_stepdown_impro*daly_weight,     
  yld_bact_culneg_nochange_impro <- los_bact_culneg_nochange_impro*daly_weight, 
  yld_bact_culneg_stepup_uncha <- los_bact_culneg_stepup_uncha*daly_weight,   
  yld_bact_culneg_nochange_uncha <- los_bact_culneg_nochange_uncha*daly_weight,   
  ## No microbiology testing performed: in-hospital length of hospitalisation based on clinical condition
  yld_bact_nocul_stepup_deter <- los_bact_nocul_stepup_deter*daly_weight,                
  yld_bact_nocul_nochange_deter <- los_bact_nocul_nochange_deter*daly_weight, 
  yld_bact_nocul_stepdown_impro <- los_bact_nocul_stepdown_impro*daly_weight,     
  yld_bact_nocul_nochange_impro <- los_bact_nocul_nochange_impro*daly_weight, 
  yld_bact_nocul_stepup_uncha <- los_bact_nocul_stepup_uncha*daly_weight,   
  yld_bact_nocul_nochange_uncha <- los_bact_nocul_nochange_uncha*daly_weight,   
  
  # Patient with no bacterial infection
  ## Culture positive: in-hospital length of hospitalisation after step-down, step-up, or no change
  yld_nobact_culpos_stepdown <- los_nobact_culpos_stepdown*daly_weight,       # expectation of life at age 55-59 years
  yld_nobact_culpos_stepup <- los_nobact_culpos_stepup*daly_weight,         
  yld_nobact_culpos_nochange <- los_nobact_culpos_nochange*daly_weight,   
  ## Culture negative: in-hospital length of hospitalisation based on clinical condition (for culture negatives)
  yld_nobact_culneg_stepup_deter <- los_nobact_culneg_stepup_deter*daly_weight,                
  yld_nobact_culneg_nochange_deter <- los_nobact_culneg_nochange_deter*daly_weight, 
  yld_nobact_culneg_stepdown_impro <- los_nobact_culneg_stepdown_impro*daly_weight,     
  yld_nobact_culneg_nochange_impro <- los_nobact_culneg_nochange_impro*daly_weight, 
  yld_nobact_culneg_stepup_uncha <- los_nobact_culneg_stepup_uncha*daly_weight,   
  yld_nobact_culneg_nochange_uncha <- los_nobact_culneg_nochange_uncha*daly_weight,   
  ## No microbiology testing performed: in-hospital length of hospitalisation based on clinical condition
  yld_nobact_nocul_stepup_deter <- los_nobact_nocul_stepup_deter*daly_weight,                
  yld_nobact_nocul_nochange_deter <- los_nobact_nocul_nochange_deter*daly_weight, 
  yld_nobact_nocul_stepdown_impro <- los_nobact_nocul_stepdown_impro*daly_weight,     
  yld_nobact_nocul_nochange_impro <- los_nobact_nocul_nochange_impro*daly_weight, 
  yld_nobact_nocul_stepup_uncha <- los_nobact_nocul_stepup_uncha*daly_weight,   
  yld_nobact_nocul_nochange_uncha <- los_nobact_nocul_nochange_uncha*daly_weight
)

