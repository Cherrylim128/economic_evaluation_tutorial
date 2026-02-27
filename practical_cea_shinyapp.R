#'Anti-microbial resistance in LMIC- practical applications in R
#'Ben S Cooper, Cherry Lim
#'
#'ShinyApp for the decision-tree model
#'This script is written by Cherry Lim with helps from ChatGPT

rm(list = ls())

# install and load essential libraries
load.lib<-c("shiny", "bslib", "ggplot2","DT","dplyr","gtools","scales") #gtools: dirichlet
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

# ---------------------------
# dec_tree() — decision tree
# ---------------------------
dec_tree <- function(params){
  with(
    as.list(params),
    {
      # Active microbiology lab service arm
      ep1 <- p.bact * p.bact_culture_pos * p.CefGen_culpos_stepdown
      ep2 <- p.bact * p.bact_culture_pos * p.CefGen_culpos_stepup
      ep3 <- p.bact * p.bact_culture_pos * p.CefGen_culpos_nochange
      ep4 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_deter * p.CefGen_culneg_deter_stepup
      ep5 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_deter * p.CefGen_culneg_deter_nochange
      ep6 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_impro * p.CefGen_culneg_impro_stepdown
      ep7 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_impro * p.CefGen_culneg_impro_nochange
      ep8 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_nocha * p.CefGen_culneg_nocha_stepup
      ep9 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_nocha * p.CefGen_culneg_nocha_nochange
      
      ep10 <- p.nobact * p.nobact_culture_pos * p.CefGen_culpos_stepdown
      ep11 <- p.nobact * p.nobact_culture_pos * p.CefGen_culpos_stepup
      ep12 <- p.nobact * p.nobact_culture_pos * p.CefGen_culpos_nochange
      ep13 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_deter * p.culneg_deter_stepup
      ep14 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_deter * p.culneg_deter_nochange
      ep15 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_impro * p.culneg_impro_stepdown
      ep16 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_impro * p.culneg_impro_nochange
      ep17 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_nocha * p.culneg_nocha_stepup
      ep18 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_nocha * p.culneg_nocha_nochange
      
      activelab_CefGen_ptotal <- ep1+ep2+ep3+ep4+ep5+ep6+
        ep7+ep8+ep9+ep10+ep11+ep12+
        ep13+ep14+ep15+ep16+ep17+ep18
      
      tc1 <- cost_cultpos + cost_oralatb*length_targeted
      tc2 <- cost_cultpos + cost_Merope*length_targeted
      tc3 <- cost_cultpos + cost_CefGen*length_targeted
      tc4 <- cost_cultneg + (cost_caredetrio + cost_Merope)*length_targeted
      tc5 <- cost_cultneg + (cost_caredetrio + cost_CefGen)*length_targeted
      tc6 <- cost_cultneg + (cost_careimprov + cost_oralatb)*length_targeted
      tc7 <- cost_cultneg + (cost_careimprov + cost_CefGen)*length_targeted
      tc8 <- cost_cultneg + (cost_careunchan + cost_Merope)*length_targeted
      tc9 <- cost_cultneg + (cost_careunchan + cost_CefGen)*length_targeted
      
      etc.activelab_CefGen <- ((ep1 * tc1) + (ep2 * tc2) + (ep3 * tc3) +
                                 (ep4 * tc4) + (ep5 * tc5) +
                                 (ep6 * tc6) + (ep7 * tc7) +
                                 (ep8 * tc8) + (ep9 * tc9) +
                                 (ep10 * tc1) + (ep11 * tc2) + (ep12 * tc3) +
                                 (ep13 * tc4) + (ep14 * tc5) +
                                 (ep15 * tc6) + (ep16 * tc7) +
                                 (ep17 * tc8) + (ep18 * tc9) +
                                 (cost_CefGen * length_empiric) + cost_henced_microlab * num_bottle) * N
      
      mort.activelab_CefGen <- ((ep1 * mort_bact_culpos_stepdown) + (ep2 * mort_bact_culpos_stepup) + (ep3 * mort_bact_culpos_nochange) +
                                  (ep4 * mort_bact_culneg_stepup_deter) + (ep5 * mort_bact_culneg_nochange_deter) +
                                  (ep6 * mort_bact_culneg_stepdown_impro) + (ep7 * mort_bact_culneg_nochange_impro) +
                                  (ep8 * mort_bact_culneg_stepup_uncha) + (ep9 * mort_bact_culneg_nochange_uncha) +
                                  (ep10 * mort_nobact_culpos_stepdown) + (ep11 * mort_nobact_culpos_stepup) + (ep12 * mort_nobact_culpos_nochange) +
                                  (ep13 * mort_nobact_culneg_stepup_deter) + (ep14 * mort_nobact_culneg_nochange_deter) +
                                  (ep15 * mort_nobact_culneg_stepdown_impro) + (ep16 * mort_nobact_culneg_nochange_impro) +
                                  (ep17 * mort_nobact_culneg_stepup_uncha) + (ep18 * mort_nobact_culneg_nochange_uncha) ) * N
      
      daly.activelab_CefGen <-
        ((ep1 * (yll_bact_culpos_stepdown+ yld_bact_culpos_stepdown)) + (ep2 * (yll_bact_culpos_stepup+ yld_bact_culpos_stepup)) + (ep3 * (yll_bact_culpos_nochange+ yld_bact_culpos_nochange)) +
           (ep4 * (yll_bact_culneg_stepup_deter+ yld_bact_culneg_stepup_deter)) + (ep5 * (yll_bact_culneg_nochange_deter+ yld_bact_culneg_nochange_deter)) +
           (ep6 * (yll_bact_culneg_stepdown_impro+ yld_bact_culneg_stepdown_impro)) + (ep7 * (yll_bact_culneg_nochange_impro+ yld_bact_culneg_nochange_impro)) +
           (ep8 * (yll_bact_culneg_stepup_uncha+ yld_bact_culneg_stepup_uncha)) + (ep9 * (yll_bact_culneg_nochange_uncha+ yld_bact_culneg_nochange_uncha)) +
           (ep10 * (yll_nobact_culpos_stepdown+ yld_nobact_culpos_stepdown)) + (ep11 * (yll_nobact_culpos_stepup+ yld_nobact_culpos_stepup)) + (ep12 * (yll_nobact_culpos_nochange+ yld_nobact_culpos_nochange)) +
           (ep13 * (yll_nobact_culneg_stepup_deter+ yld_nobact_culneg_stepup_deter)) + (ep14 * (yll_nobact_culneg_nochange_deter+ yld_nobact_culneg_nochange_deter)) +
           (ep15 * (yll_nobact_culneg_stepdown_impro+ yld_nobact_culneg_stepdown_impro)) + (ep16 * (yll_nobact_culneg_nochange_impro+ yld_nobact_culneg_nochange_impro)) +
           (ep17 * (yll_nobact_culneg_stepup_uncha+ yld_nobact_culneg_stepup_uncha)) + (ep18 * (yll_nobact_culneg_nochange_uncha+ yld_nobact_culneg_nochange_uncha)))*N
      
      # No microbiology lab arm
      np1 <- p.bact * p.CefGen_bact_nocul_deter * p.CefGen_nocul_deter_stepup
      np2 <- p.bact * p.CefGen_bact_nocul_deter * p.CefGen_nocul_deter_nochange
      np3 <- p.bact * p.CefGen_bact_nocul_impro * p.CefGen_nocul_impro_stepdown
      np4 <- p.bact * p.CefGen_bact_nocul_impro * p.CefGen_nocul_impro_nochange
      np5 <- p.bact * p.CefGen_bact_nocul_nocha * p.CefGen_nocul_nocha_stepup
      np6 <- p.bact * p.CefGen_bact_nocul_nocha * p.CefGen_nocul_nocha_nochange
      
      np7 <- p.nobact * p.CefGen_nobact_nocul_deter * p.CefGen_nocul_deter_stepup
      np8 <- p.nobact * p.CefGen_nobact_nocul_deter * p.CefGen_nocul_deter_nochange
      np9 <- p.nobact * p.CefGen_nobact_nocul_impro * p.CefGen_nocul_impro_stepdown
      np10 <- p.nobact * p.CefGen_nobact_nocul_impro * p.CefGen_nocul_impro_nochange
      np11 <- p.nobact * p.CefGen_nobact_nocul_nocha * p.CefGen_nocul_nocha_stepup
      np12 <- p.nobact * p.CefGen_nobact_nocul_nocha * p.CefGen_nocul_nocha_nochange
      
      noactivelab_CefGen_ptotal <- np1+np2+np3+np4+np5+np6+
        np7+np8+np9+np10+np11+np12
      
      nstc1 <- (cost_caredetrio + cost_Merope)*length_targeted
      nstc2 <- (cost_caredetrio + cost_CefGen)*length_targeted
      nstc3 <- (cost_careimprov + cost_oralatb)*length_targeted
      nstc4 <- (cost_careimprov + cost_CefGen)*length_targeted
      nstc5 <- (cost_careunchan + cost_Merope)*length_targeted
      nstc6 <- (cost_careunchan + cost_CefGen)*length_targeted
      
      etc.noactivelab_CefGen <-
        ((np1 * nstc1) + (np2 * nstc2) +
           (np3 * nstc3) + (np4 * nstc4) +
           (np5 * nstc5) + (np6 * nstc6) +
           (np7 * nstc1) + (np8 * nstc2) +
           (np9 * nstc3) + (np10 * nstc4) +
           (np11 * nstc5) + (np12 * nstc6) +
           cost_CefGen*length_empiric) * N
      
      mort.noactivelab_CefGen <-
        ((np1 * mort_bact_nocul_stepup_deter) + (np2 * mort_bact_nocul_nochange_deter) +
           (np3 * mort_bact_nocul_stepdown_impro) + (np4 * mort_bact_nocul_nochange_impro) +
           (np5 * mort_bact_nocul_stepup_uncha) + (np6 * mort_bact_nocul_nochange_uncha) +
           (np7 * mort_nobact_nocul_stepup_deter) + (np8 * mort_nobact_nocul_nochange_deter) +
           (np9 * mort_nobact_nocul_stepdown_impro) + (np10 * mort_nobact_nocul_nochange_impro) +
           (np11 * mort_nobact_nocul_stepup_uncha) + (np12 * mort_nobact_nocul_nochange_uncha)) * N
      
      daly.noactivelab_CefGen <-
        ((np1 * (yll_bact_nocul_stepup_deter+ yld_bact_nocul_stepup_deter)) + (np2 * (yll_bact_nocul_nochange_deter+ yld_bact_nocul_nochange_deter)) +
           (np3 * (yll_bact_nocul_stepdown_impro+ yld_bact_nocul_stepdown_impro)) + (np4 * (yll_bact_nocul_nochange_impro+ yld_bact_nocul_nochange_impro)) +
           (np5 * (yll_bact_nocul_stepup_uncha+ yld_bact_nocul_stepup_uncha)) + (np6 * (yll_bact_nocul_nochange_uncha+ yld_bact_nocul_nochange_uncha)) +
           (np7 * (yll_nobact_nocul_stepup_deter+ yld_nobact_nocul_stepup_deter)) + (np8 * (yll_nobact_nocul_nochange_deter+ yld_nobact_nocul_nochange_deter)) +
           (np9 * (yll_nobact_nocul_stepdown_impro+ yld_nobact_nocul_stepdown_impro)) + (np10 * (yll_nobact_nocul_nochange_impro+ yld_nobact_nocul_nochange_impro)) +
           (np11 * (yll_nobact_nocul_stepup_uncha+ yld_nobact_nocul_stepup_uncha)) + (np12 * (yll_nobact_nocul_nochange_uncha+ yld_nobact_nocul_nochange_uncha))) * N
      
      im <- mort.noactivelab_CefGen - mort.activelab_CefGen
      ic <- etc.activelab_CefGen - etc.noactivelab_CefGen
      ie <- daly.noactivelab_CefGen - daly.activelab_CefGen
      icer <- ic/ie
      
      return(c(im, ic, ie, icer,
               mort.noactivelab_CefGen, mort.activelab_CefGen,
               etc.noactivelab_CefGen, etc.activelab_CefGen,
               daly.noactivelab_CefGen, daly.activelab_CefGen))
    }
  )
}

# ---------------------------
# generate_input_prob() — creates sampling for dec_tree
# ---------------------------
generate_input_prob <- function(
    use_sampling = TRUE,
    # infection & culture positivity fixed values (used when use_sampling=FALSE)
    # This function is turned-off!!
    p_infecttype_value = 0.40,
    p_bact_culture_pos_value = 0.80,
    p_nobact_culture_pos_value = 0.05,
    # culpos dirichlet alphas for culture positive branches
    culpos_dir_alphas = c(6.1, 39, 54.9),
    # culneg dirichlet alphas for no-culture
    bact_nocul_dir_alphas = c(24, 60, 16),
    nobact_nocul_dir_alphas = c(7, 37, 56),
    # culneg beta distribution parameters for negative culture
    cefneg_deter_stepup_alpha = 70,
    cefneg_impro_stepdown_alpha = 22,
    cefneg_nocha_stepup_alpha = 56,
    culneg_deter_stepup_alpha = 1.2,
    culneg_impro_stepdown_alpha = 80,
    culneg_nocha_stepup_alpha = 1.3,
    # costs
    cost_henced_microlab_value = 300,
    cost_CefGen_value = 1.76,
    cost_Merope_value = 35.16,
    cost_oralatb_value = 0.14,
    cost_cultpos_value = 120,
    cost_cultneg_value = 20,
    cost_careunchan_value = 158.82,
    cost_caredetrio_value = 203.19,
    cost_careimprov_value = 14.44,
    cost_sdlog = 0.25,
    # model settings
    num_bottle = 2,
    year_lost = 11,
    daly_weight = 0.133,
    length_empiric = 3,
    length_targeted = 7,
    N = 1000
) {
  # infection status and culture results
  if (use_sampling) {
    p.infecttype <- rbeta(1, 40, 60) # proportion of true bacterial infections
    p.bact_culture_pos <- rbeta(1, 80, 20) # proportion of culture positive among infected
    p.nobact_culture_pos <- rbeta(1, 5, 95) # proportion of culture positive among non-infected
  } else {
    p.infecttype <- p_infecttype_value
    p.bact_culture_pos <- p_bact_culture_pos_value
    p.nobact_culture_pos <- p_nobact_culture_pos_value
  }
  
  # postive culture
  if (use_sampling) {
    p.CefGen_culpos_atb <- as.numeric(gtools::rdirichlet(1, culpos_dir_alphas))
  } else {
    p.CefGen_culpos_atb <- culpos_dir_alphas / sum(culpos_dir_alphas)
  }
  
  # No culture
  p.CefGen_bact_nocul_dir <- as.numeric(gtools::rdirichlet(1, bact_nocul_dir_alphas))
  p.CefGen_nobact_nocul_dir <- as.numeric(gtools::rdirichlet(1, nobact_nocul_dir_alphas))
  
  # Bacterial infection and negative cultures: antibiotic changes based on clinical presentation
  # guard against alpha>=100 by using max(1,100-alpha) for beta parameter
  p.CefGen_culneg_deter_stepup <- rbeta(1, cefneg_deter_stepup_alpha, max(1, 100 - cefneg_deter_stepup_alpha))
  p.CefGen_culneg_deter_nochange <- 1 - p.CefGen_culneg_deter_stepup
  p.CefGen_culneg_impro_stepdown <- rbeta(1, cefneg_impro_stepdown_alpha, max(1, 100 - cefneg_impro_stepdown_alpha))
  p.CefGen_culneg_impro_nochange <- 1 - p.CefGen_culneg_impro_stepdown
  p.CefGen_culneg_nocha_stepup <- rbeta(1, cefneg_nocha_stepup_alpha, max(1, 100 - cefneg_nocha_stepup_alpha))
  p.CefGen_culneg_nocha_nochange <- 1 - p.CefGen_culneg_nocha_stepup
  
  # Not bacterial infection and negative culture: antibiotic changes based on clinical presentation
  p.culneg_deter_stepup <- rbeta(1, culneg_deter_stepup_alpha, max(1, 100 - culneg_deter_stepup_alpha))
  p.culneg_deter_nochange <- 1 - p.culneg_deter_stepup
  p.culneg_impro_stepdown <- rbeta(1, culneg_impro_stepdown_alpha, max(1, 100 - culneg_impro_stepdown_alpha))
  p.culneg_impro_nochange <- 1 - p.culneg_impro_stepdown
  p.culneg_nocha_stepup <- rbeta(1, culneg_nocha_stepup_alpha, max(1, 100 - culneg_nocha_stepup_alpha))
  p.culneg_nocha_nochange <- 1 - p.culneg_nocha_stepup
  
  # costs
  if (use_sampling) {
    cost_careunchan <- rlnorm(1, log(cost_careunchan_value), cost_sdlog)
    cost_caredetrio <- rlnorm(1, log(cost_caredetrio_value), cost_sdlog)
    cost_careimprov <- rlnorm(1, log(cost_careimprov_value), cost_sdlog)
  } else {
    cost_careunchan <- cost_careunchan_value
    cost_caredetrio <- cost_caredetrio_value
    cost_careimprov <- cost_careimprov_value
  }
  
  # Compose params list — must include every p.* used by dec_tree()
  params <- list(
    # infection status and culture
    p.bact = p.infecttype,
    p.nobact = 1 - p.infecttype,
    p.bact_culture_pos = p.bact_culture_pos,
    p.nobact_culture_pos = p.nobact_culture_pos,
    
    # culture positive
    p.CefGen_culpos_stepdown = p.CefGen_culpos_atb[1],
    p.CefGen_culpos_stepup   = p.CefGen_culpos_atb[2],
    p.CefGen_culpos_nochange = p.CefGen_culpos_atb[3],
    
    # clinical presentation of culture negative 
    # Here we set bacterial culneg proportions so they sum to 1 across deter/impro/nocha.
    p.CefGen_bact_culneg_deter = 0.048,   
    p.CefGen_bact_culneg_impro = 0.792,
    p.CefGen_bact_culneg_nocha = 0.16,
    p.CefGen_nobact_culneg_deter = 0.07,
    p.CefGen_nobact_culneg_impro = 0.37,
    p.CefGen_nobact_culneg_nocha = 0.56,
    
    # antibiotic changes based on clinical presentation of those with culture negative
    p.CefGen_culneg_deter_stepup = p.CefGen_culneg_deter_stepup,
    p.CefGen_culneg_deter_nochange = p.CefGen_culneg_deter_nochange,
    p.CefGen_culneg_impro_stepdown = p.CefGen_culneg_impro_stepdown,
    p.CefGen_culneg_impro_nochange = p.CefGen_culneg_impro_nochange,
    p.CefGen_culneg_nocha_stepup = p.CefGen_culneg_nocha_stepup,
    p.CefGen_culneg_nocha_nochange = p.CefGen_culneg_nocha_nochange,
    
    # no-culture performed
    p.CefGen_bact_nocul_deter = p.CefGen_bact_nocul_dir[1],
    p.CefGen_bact_nocul_impro = p.CefGen_bact_nocul_dir[2],
    p.CefGen_bact_nocul_nocha = p.CefGen_bact_nocul_dir[3],
    
    p.CefGen_nobact_nocul_deter = p.CefGen_nobact_nocul_dir[1],
    p.CefGen_nobact_nocul_impro = p.CefGen_nobact_nocul_dir[2],
    p.CefGen_nobact_nocul_nocha = p.CefGen_nobact_nocul_dir[3],
    
    # antibiotic changes of those with no-culture performed
    p.CefGen_nocul_deter_stepup = p.culneg_deter_stepup,
    p.CefGen_nocul_deter_nochange = p.culneg_deter_nochange,
    p.CefGen_nocul_impro_stepdown = p.culneg_impro_stepdown,
    p.CefGen_nocul_impro_nochange = p.culneg_impro_nochange,
    p.CefGen_nocul_nocha_stepup = p.culneg_nocha_stepup,
    p.CefGen_nocul_nocha_nochange = p.culneg_nocha_nochange,
    
    # antibiotic changes of those with negative culture 
    p.culneg_deter_stepup = p.culneg_deter_stepup,
    p.culneg_deter_nochange = p.culneg_deter_nochange,
    p.culneg_impro_stepdown = p.culneg_impro_stepdown,
    p.culneg_impro_nochange = p.culneg_impro_nochange,
    p.culneg_nocha_stepup = p.culneg_nocha_stepup,
    p.culneg_nocha_nochange = p.culneg_nocha_nochange,
    
    # costs
    cost_henced_microlab = cost_henced_microlab_value,
    cost_CefGen = cost_CefGen_value,
    cost_Merope = cost_Merope_value,
    cost_oralatb = cost_oralatb_value,
    cost_cultpos = cost_cultpos_value,
    cost_cultneg = cost_cultneg_value,
    cost_careunchan = cost_careunchan,
    cost_caredetrio = cost_caredetrio,
    cost_careimprov = cost_careimprov,
    
    # model settings
    num_bottle = num_bottle,
    length_empiric = length_empiric,
    length_targeted = length_targeted,
    N = N,
    year_lost = year_lost,
    daly_weight = daly_weight
  )
  
  # append vectors
  params$p.CefGen_bact_nocul_dir <- p.CefGen_bact_nocul_dir
  params$p.CefGen_nobact_nocul_dir <- p.CefGen_nobact_nocul_dir
  params$p.CefGen_culpos_atb <- p.CefGen_culpos_atb
  
  # LOS 
  params$los_bact_culpos_stepdown <- rlnorm(1, meanlog = log(7), sdlog = 0.1)
  params$los_bact_culpos_stepup <- params$los_bact_culpos_stepdown * 1.5
  params$los_bact_culpos_nochange <- params$los_bact_culpos_stepdown
  
  params$los_bact_culneg_stepup_deter <- params$los_bact_culpos_stepup * 2
  params$los_bact_culneg_nochange_deter <- params$los_bact_culneg_stepup_deter * 1.5
  params$los_bact_culneg_stepdown_impro <- params$los_bact_culpos_stepdown
  params$los_bact_culneg_nochange_impro <- params$los_bact_culpos_stepdown
  params$los_bact_culneg_stepup_uncha <- params$los_bact_culpos_nochange
  params$los_bact_culneg_nochange_uncha <- params$los_bact_culpos_nochange
  
  params$los_bact_nocul_stepup_deter <- params$los_bact_culneg_stepup_deter * 1.25
  params$los_bact_nocul_nochange_deter <- params$los_bact_culneg_stepup_deter * 2
  params$los_bact_nocul_stepdown_impro <- params$los_bact_culneg_stepup_deter * 1.5
  params$los_bact_nocul_nochange_impro <- params$los_bact_nocul_stepup_deter
  params$los_bact_nocul_stepup_uncha <- params$los_bact_nocul_stepup_deter
  params$los_bact_nocul_nochange_uncha <- params$los_bact_culneg_stepup_deter * 1.5
  
  params$los_nobact_culpos_stepdown <- params$los_bact_culpos_stepdown / 2
  params$los_nobact_culpos_stepup <- params$los_bact_culpos_stepup / 2
  params$los_nobact_culpos_nochange <- params$los_bact_culpos_nochange / 2
  
  params$los_nobact_culneg_stepup_deter <- params$los_bact_culneg_stepup_deter / 2
  params$los_nobact_culneg_nochange_deter <- params$los_bact_culneg_nochange_deter / 2
  params$los_nobact_culneg_stepdown_impro <- params$los_bact_culneg_stepdown_impro / 2
  params$los_nobact_culneg_nochange_impro <- params$los_bact_culneg_nochange_impro / 2
  params$los_nobact_culneg_stepup_uncha <- params$los_bact_culneg_stepup_uncha / 2
  params$los_nobact_culneg_nochange_uncha <- params$los_bact_culneg_nochange_uncha / 2
  
  params$los_nobact_nocul_stepup_deter <- params$los_nobact_culneg_stepup_deter
  params$los_nobact_nocul_nochange_deter <- params$los_nobact_culneg_nochange_deter
  params$los_nobact_nocul_stepdown_impro <- params$los_nobact_culneg_stepdown_impro
  params$los_nobact_nocul_nochange_impro <- params$los_nobact_culneg_nochange_impro
  params$los_nobact_nocul_stepup_uncha <- params$los_nobact_culneg_stepup_uncha
  params$los_nobact_nocul_nochange_uncha <- params$los_nobact_culneg_nochange_uncha
  
  # Mortality
  params$mort_bact_culpos_stepdown <- rbeta(1, 20, 80)
  params$mort_bact_culpos_stepup <- params$mort_bact_culpos_stepdown * 1.25
  params$mort_bact_culpos_nochange <- params$mort_bact_culpos_stepdown * 1.25
  
  params$mort_bact_culneg_stepup_deter <- params$mort_bact_culpos_stepdown * 1.25
  params$mort_bact_culneg_nochange_deter <- params$mort_bact_culneg_stepup_deter * 1.5
  params$mort_bact_culneg_stepdown_impro <- params$mort_bact_culpos_stepdown
  params$mort_bact_culneg_nochange_impro <- params$mort_bact_culpos_stepdown
  params$mort_bact_culneg_stepup_uncha <- params$mort_bact_culpos_nochange
  params$mort_bact_culneg_nochange_uncha <- params$mort_bact_culpos_nochange
  
  params$mort_bact_nocul_stepup_deter <- params$mort_bact_culneg_stepup_deter * 1.25
  params$mort_bact_nocul_nochange_deter <- params$mort_bact_culneg_stepup_deter * 2
  params$mort_bact_nocul_stepdown_impro <- params$mort_bact_culneg_stepup_deter * 1.5
  params$mort_bact_nocul_nochange_impro <- params$mort_bact_nocul_stepup_deter
  params$mort_bact_nocul_stepup_uncha <- params$mort_bact_nocul_stepup_deter
  params$mort_bact_nocul_nochange_uncha <- params$mort_bact_culneg_stepup_deter * 1.5
  
  params$mort_nobact_culpos_stepdown <- params$mort_bact_culpos_stepdown / 2
  params$mort_nobact_culpos_stepup <- params$mort_bact_culpos_stepup / 2
  params$mort_nobact_culpos_nochange <- params$mort_bact_culpos_nochange / 2
  
  params$mort_nobact_culneg_stepup_deter <- params$mort_bact_culneg_stepup_deter / 2
  params$mort_nobact_culneg_nochange_deter <- params$mort_bact_culneg_nochange_deter / 2
  params$mort_nobact_culneg_stepdown_impro <- params$mort_bact_culneg_stepdown_impro / 2
  params$mort_nobact_culneg_nochange_impro <- params$mort_bact_culneg_nochange_impro / 2
  params$mort_nobact_culneg_stepup_uncha <- params$mort_bact_culneg_stepup_uncha / 2
  params$mort_nobact_culneg_nochange_uncha <- params$mort_bact_culneg_nochange_uncha / 2
  
  params$mort_nobact_nocul_stepup_deter <- params$mort_nobact_culneg_stepup_deter
  params$mort_nobact_nocul_nochange_deter <- params$mort_nobact_culneg_nochange_deter
  params$mort_nobact_nocul_stepdown_impro <- params$mort_nobact_culneg_stepdown_impro
  params$mort_nobact_nocul_nochange_impro <- params$mort_nobact_culneg_nochange_impro
  params$mort_nobact_nocul_stepup_uncha <- params$mort_nobact_culneg_stepup_uncha
  params$mort_nobact_nocul_nochange_uncha <- params$mort_nobact_culneg_nochange_uncha
  
  # YLL (mortality * years lost)
  params$yll_bact_culpos_stepdown <- params$mort_bact_culpos_stepdown * params$year_lost
  params$yll_bact_culpos_stepup <- params$mort_bact_culpos_stepup * params$year_lost
  params$yll_bact_culpos_nochange <- params$mort_bact_culpos_nochange * params$year_lost
  
  params$yll_bact_culneg_stepup_deter <- params$mort_bact_culneg_stepup_deter * params$year_lost
  params$yll_bact_culneg_nochange_deter <- params$mort_bact_culneg_nochange_deter * params$year_lost
  params$yll_bact_culneg_stepdown_impro <- params$mort_bact_culneg_stepdown_impro * params$year_lost
  params$yll_bact_culneg_nochange_impro <- params$mort_bact_culneg_nochange_impro * params$year_lost
  params$yll_bact_culneg_stepup_uncha <- params$mort_bact_culneg_stepup_uncha * params$year_lost
  params$yll_bact_culneg_nochange_uncha <- params$mort_bact_culneg_nochange_uncha * params$year_lost
  
  params$yll_bact_nocul_stepup_deter <- params$mort_bact_nocul_stepup_deter * params$year_lost
  params$yll_bact_nocul_nochange_deter <- params$mort_bact_nocul_nochange_deter * params$year_lost
  params$yll_bact_nocul_stepdown_impro <- params$mort_bact_nocul_stepdown_impro * params$year_lost
  params$yll_bact_nocul_nochange_impro <- params$mort_bact_nocul_nochange_impro * params$year_lost
  params$yll_bact_nocul_stepup_uncha <- params$mort_bact_nocul_stepup_uncha * params$year_lost
  params$yll_bact_nocul_nochange_uncha <- params$mort_bact_nocul_nochange_uncha * params$year_lost
  
  params$yll_nobact_culpos_stepdown <- params$mort_nobact_culpos_stepdown * params$year_lost
  params$yll_nobact_culpos_stepup <- params$mort_nobact_culpos_stepup * params$year_lost
  params$yll_nobact_culpos_nochange <- params$mort_nobact_culpos_nochange * params$year_lost
  
  params$yll_nobact_culneg_stepup_deter <- params$mort_nobact_culneg_stepup_deter * params$year_lost
  params$yll_nobact_culneg_nochange_deter <- params$mort_nobact_culneg_nochange_deter * params$year_lost
  params$yll_nobact_culneg_stepdown_impro <- params$mort_nobact_culneg_stepdown_impro * params$year_lost
  params$yll_nobact_culneg_nochange_impro <- params$mort_nobact_culneg_nochange_impro * params$year_lost
  params$yll_nobact_culneg_stepup_uncha <- params$mort_nobact_culneg_stepup_uncha * params$year_lost
  params$yll_nobact_culneg_nochange_uncha <- params$mort_nobact_culneg_nochange_uncha * params$year_lost
  
  # no-culture, no-bacterial-infection YLLs
  params$yll_nobact_nocul_stepup_deter <- params$mort_nobact_nocul_stepup_deter * params$year_lost
  params$yll_nobact_nocul_nochange_deter <- params$mort_nobact_nocul_nochange_deter * params$year_lost
  params$yll_nobact_nocul_stepdown_impro <- params$mort_nobact_nocul_stepdown_impro * params$year_lost
  params$yll_nobact_nocul_nochange_impro <- params$mort_nobact_nocul_nochange_impro * params$year_lost
  params$yll_nobact_nocul_stepup_uncha <- params$mort_nobact_nocul_stepup_uncha * params$year_lost
  params$yll_nobact_nocul_nochange_uncha <- params$mort_nobact_nocul_nochange_uncha * params$year_lost
  
  # YLD (LOS * daly_weight)
  params$yld_bact_culpos_stepdown <- params$los_bact_culpos_stepdown * params$daly_weight
  params$yld_bact_culpos_stepup <- params$los_bact_culpos_stepup * params$daly_weight
  params$yld_bact_culpos_nochange <- params$los_bact_culpos_nochange * params$daly_weight
  params$yld_bact_culneg_stepup_deter <- params$los_bact_culneg_stepup_deter * params$daly_weight
  params$yld_bact_culneg_nochange_deter <- params$los_bact_culneg_nochange_deter * params$daly_weight
  params$yld_bact_culneg_stepdown_impro <- params$los_bact_culneg_stepdown_impro * params$daly_weight
  params$yld_bact_culneg_nochange_impro <- params$los_bact_culneg_nochange_impro * params$daly_weight
  params$yld_bact_culneg_stepup_uncha <- params$los_bact_culneg_stepup_uncha * params$daly_weight
  params$yld_bact_culneg_nochange_uncha <- params$los_bact_culneg_nochange_uncha * params$daly_weight
  
  params$yld_bact_nocul_stepup_deter <- params$los_bact_nocul_stepup_deter * params$daly_weight
  params$yld_bact_nocul_nochange_deter <- params$los_bact_nocul_nochange_deter * params$daly_weight
  params$yld_bact_nocul_stepdown_impro <- params$los_bact_nocul_stepdown_impro * params$daly_weight
  params$yld_bact_nocul_nochange_impro <- params$los_bact_nocul_nochange_impro * params$daly_weight
  params$yld_bact_nocul_stepup_uncha <- params$los_bact_nocul_stepup_uncha * params$daly_weight
  params$yld_bact_nocul_nochange_uncha <- params$los_bact_nocul_nochange_uncha * params$daly_weight
  
  params$yld_nobact_culpos_stepdown <- params$los_nobact_culpos_stepdown * params$daly_weight
  params$yld_nobact_culpos_stepup <- params$los_nobact_culpos_stepup * params$daly_weight
  params$yld_nobact_culpos_nochange <- params$los_nobact_culpos_nochange * params$daly_weight
  
  params$yld_nobact_culneg_stepup_deter <- params$los_nobact_culneg_stepup_deter * params$daly_weight
  params$yld_nobact_culneg_nochange_deter <- params$los_nobact_culneg_nochange_deter * params$daly_weight
  params$yld_nobact_culneg_stepdown_impro <- params$los_nobact_culneg_stepdown_impro * params$daly_weight
  params$yld_nobact_culneg_nochange_impro <- params$los_nobact_culneg_nochange_impro * params$daly_weight
  params$yld_nobact_culneg_stepup_uncha <- params$los_nobact_culneg_stepup_uncha * params$daly_weight
  params$yld_nobact_culneg_nochange_uncha <- params$los_nobact_culneg_nochange_uncha * params$daly_weight
  
  params$yld_nobact_nocul_stepup_deter <- params$los_nobact_nocul_stepup_deter * params$daly_weight
  params$yld_nobact_nocul_nochange_deter <- params$los_nobact_nocul_nochange_deter * params$daly_weight
  params$yld_nobact_nocul_stepdown_impro <- params$los_nobact_nocul_stepdown_impro * params$daly_weight
  params$yld_nobact_nocul_nochange_impro <- params$los_nobact_nocul_nochange_impro * params$daly_weight
  params$yld_nobact_nocul_stepup_uncha <- params$los_nobact_nocul_stepup_uncha * params$daly_weight
  params$yld_nobact_nocul_nochange_uncha <- params$los_nobact_nocul_nochange_uncha * params$daly_weight
  
  return(params)
}

# ---------------------------
# UI — left / center / right layout
# ---------------------------
ui <- fluidPage(
  theme = bs_theme(version = 4, bootswatch = "flatly", primary = "#1565C0"),
  titlePanel("Cost-effectiveness analysis of maintaining active microbiology laboratory services"),
  
  fluidRow(
    # LEFT column: inputs
    column(
      width = 3,
      div(style = "max-height:85vh; overflow-y:auto; padding-right: 12px;",
          wellPanel(
            h4("General parameters", style = "color:#0D47A1"),
            # checkboxInput("use_sampling", "Use sampling (PSA)", value = TRUE),
            numericInput("niter", "Number of iterations", value = 1000, min = 100, step = 100),
            numericInput("N", "Population size (N)", value = 1000, min = 1),
            numericInput("num_bottle", "Number of blood bottles taken per patient", value = 2, min = 1),
            numericInput("length_empiric", "Length of empiric antibiotic treatment (days)", value = 3, min = 1),
            numericInput("length_targeted", "Length of targeted antibiotic treatment (days)", value = 7, min = 1),
            numericInput("year_lost", "Number of years lost", value = 11, min = 1),
            numericInput("daly_weight", "DALY weight", value = 0.133, step = 0.001),
            numericInput("wtp", "Willingness-to-pay (USD per DALY)", value = 500, min = 0, step = 10),
            hr(),
            
            h5("Infection status and culture results", style = "color:#0D47A1"),
            numericInput("p_infecttype_value", "Probability of that the suspected BSI case truly has a bacterial infection", value = 0.40, min = 0, max = 1, step = 0.01),
            sliderInput("p_bact_culture_pos_value", "Probability of culture positive given if the suspected BSI case has a bacterial infection", min = 0, max = 1, value = 0.8, step = 0.01),
            sliderInput("p_nobact_culture_pos_value", "Probability of culture positive given if the suspected BSI case has no bacterial infection", min = 0, max = 1, value = 0.05, step = 0.01),
            hr(),
            
            h5("Patients with culture positive", style = "color:#0D47A1"),
            # p("Edit these Dirichlet alphas"),
            numericInput("culpos_a1", "Parameter to reflect probability of stepping down given if culture result is positive", value = 6.1, min = 0),
            numericInput("culpos_a2", "Parameter to reflect probability of stepping up given if culture result is positive", value = 39, min = 0),
            numericInput("culpos_a3", "Parameter to reflect probability of no change given if culture result is positive", value = 54.9, min = 0),
            hr(),
            
            h5("Health condition given has bacterial infection and no culture done", style = "color:#0D47A1"),
            numericInput("bact_nocul_a1", "Parameter to reflect probability that patient’s health condition deteriorates", value = 24, min = 0),
            numericInput("bact_nocul_a2", "Parameter to reflect probability that patient’s health condition improves", value = 60, min = 0),
            numericInput("bact_nocul_a3", "Parameter to reflect probability that patient’s health condition remains unchanged", value = 16, min = 0),
            hr(),
            
            h5("Health condition given has no bacterial infection and no culture done", style = "color:#0D47A1"),
            numericInput("nobact_nocul_a1", "Parameter to reflect probability that patient’s health condition deteriorates", value = 7, min = 0),
            numericInput("nobact_nocul_a2", "Parameter to reflect probability that patient’s health condition improves", value = 37, min = 0),
            numericInput("nobact_nocul_a3", "Parameter to reflect probability that patient’s health condition remains unchanged", value = 56, min = 0),
            hr(),
          )
      )
    ),
    
    # CENTER column: results (CE-plane + tabs)
    column(
      width = 6,
      wellPanel(
        tabsetPanel(
          tabPanel("CE plane",
                   br(),
                   plotOutput("psa_ceplane", height = "600px"),
                   br(),
                   verbatimTextOutput("psa_summary")
          ),
          tabPanel("Distributions",
                   br(),
                   fluidRow(
                     column(6, plotOutput("psa_cost_hist", height = "360px")),
                     column(6, plotOutput("psa_daly_hist", height = "360px"))
                   ),
                   br(),
                   DTOutput("psa_table")
          )
        )
      )
    ),
    
    # RIGHT column
    column(
      width = 3,
      div(style = "max-height:85vh; overflow-y:auto; padding-left: 12px;",
          wellPanel(
            h4("Changes in antibiotic therapy under different health conditions given culture result is negative", style = "color:#0D47A1"),
            sliderInput("culneg_deter_stepup_alpha", "Probability of stepping up given patient’s health condition deteriorates (%)", min = 1, max = 99, value = 70),
            sliderInput("culneg_impro_stepdown_alpha", "Probability of stepping down given patient’s health condition improves (%)", min = 1, max = 99, value = 80),
            sliderInput("culneg_nocha_stepup_alpha", "Probability of stepping up given patient’s health condition remains unchanged (%)", min = 1, max = 99, value = 10),
            hr(),
            
            h4("Changes in antibiotic therapy under different health conditions given no culture was performed", style = "color:#0D47A1"),
            sliderInput("cefneg_deter_stepup_alpha", "Probability of stepping up given patient’s health condition deteriorates (%)", min = 1, max = 99, value = 90),
            sliderInput("cefneg_impro_stepdown_alpha", "Probability of stepping down given patient’s health condition improves (%)", min = 1, max = 99, value = 20),
            sliderInput("cefneg_nocha_stepup_alpha", "Probability of stepping up given patient’s health condition remains unchanged (%)", min = 1, max = 99, value = 70),
            hr(),
            
            h4("Cost inputs (USD)", style = "color:#0D47A1"),
            sliderInput("cost_henced_microlab", "Maintaining microbiology lab cost per blood bottle", min = 0, max = 2000, value = 300, step = 10, pre = "$"),
            sliderInput("cost_CefGen", "1-day ceftriaxone/gentamicin", min = 0, max = 10, value = 1.76, step = 0.1, pre = "$"),
            sliderInput("cost_Merope", "1-day meropenem/vancomycin", min = 0, max = 200, value = 35.16, step = 0.1, pre = "$"),
            sliderInput("cost_oralatb", "1-day oral amoxicillin", min = 0, max = 20, value = 0.14, step = 0.01, pre = "$"),
            sliderInput("cost_cultpos", "Processing positive cultures", min = 0, max = 500, value = 120, step = 1, pre = "$"),
            sliderInput("cost_cultneg", "Processing negative cultures", min = 0, max = 200, value = 20, step = 1, pre = "$"),
            sliderInput("cost_careunchan", "Per-day care when the condition is unchanged", min = 0, max = 1000, value = 158.82, step = 1, pre = "$"),
            sliderInput("cost_caredetrio", "Per-day care when the condition deteriorated", min = 0, max = 1000, value = 203.19, step = 1, pre = "$"),
            sliderInput("cost_careimprov", "Per-day care when the condition improved", min = 0, max = 500, value = 14.44, step = 1, pre = "$"),
            hr(),
            downloadButton("download_psa", "Download raw PSA results")
          )
      )
    )
  ) # end fluidRow
) # end ui

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {
  
  # Build reactive UI params (only reference inputs that exist)
  ui_params <- reactive({
    list(
      use_sampling = isTRUE(input$use_sampling),
      p_infecttype_value = input$p_infecttype_value,
      p_bact_culture_pos_value = input$p_bact_culture_pos_value,
      p_nobact_culture_pos_value = input$p_nobact_culture_pos_value,
      culpos_dir_alphas = c(input$culpos_a1, input$culpos_a2, input$culpos_a3),
      cefneg_deter_stepup_alpha = input$cefneg_deter_stepup_alpha,
      cefneg_impro_stepdown_alpha = input$cefneg_impro_stepdown_alpha,
      cefneg_nocha_stepup_alpha = input$cefneg_nocha_stepup_alpha,
      culneg_deter_stepup_alpha = input$culneg_deter_stepup_alpha,
      culneg_impro_stepdown_alpha = input$culneg_impro_stepdown_alpha,
      culneg_nocha_stepup_alpha = input$culneg_nocha_stepup_alpha,
      bact_nocul_alphas = c(input$bact_nocul_a1, input$bact_nocul_a2, input$bact_nocul_a3),
      nobact_nocul_alphas = c(input$nobact_nocul_a1, input$nobact_nocul_a2, input$nobact_nocul_a3),
      cost_vals = list(
        cost_henced_microlab = input$cost_henced_microlab,
        cost_CefGen = input$cost_CefGen,
        cost_Merope = input$cost_Merope,
        cost_oralatb = input$cost_oralatb,
        cost_cultpos = input$cost_cultpos,
        cost_cultneg = input$cost_cultneg,
        cost_careunchan = input$cost_careunchan,
        cost_caredetrio = input$cost_caredetrio,
        cost_careimprov = input$cost_careimprov
      ),
      N = input$N,
      num_bottle = input$num_bottle,
      length_empiric = input$length_empiric,
      length_targeted = input$length_targeted,
      year_lost = input$year_lost,
      daly_weight = input$daly_weight,
      niter = input$niter,
      wtp = input$wtp
    )
  })
  
  # Reactive PSA results: re-runs automatically when inputs change
  results_prob <- reactive({
    
    params_ui <- ui_params()
    niter <- params_ui$niter
    
    # For performance: if niter very large, consider reducing or debouncing ui_params
    out_list <- vector("list", niter)
    
    for (i in seq_len(niter)) {
      params <- generate_input_prob(
        use_sampling = params_ui$use_sampling,
        p_infecttype_value = params_ui$p_infecttype_value,
        p_bact_culture_pos_value = params_ui$p_bact_culture_pos_value,
        p_nobact_culture_pos_value = params_ui$p_nobact_culture_pos_value,
        culpos_dir_alphas = params_ui$culpos_dir_alphas,
        bact_nocul_dir_alphas = params_ui$bact_nocul_alphas,
        nobact_nocul_dir_alphas = params_ui$nobact_nocul_alphas,
        cefneg_deter_stepup_alpha = params_ui$cefneg_deter_stepup_alpha,
        cefneg_impro_stepdown_alpha = params_ui$cefneg_impro_stepdown_alpha,
        cefneg_nocha_stepup_alpha = params_ui$cefneg_nocha_stepup_alpha,
        culneg_deter_stepup_alpha = params_ui$culneg_deter_stepup_alpha,
        culneg_impro_stepdown_alpha = params_ui$culneg_impro_stepdown_alpha,
        culneg_nocha_stepup_alpha = params_ui$culneg_nocha_stepup_alpha,
        cost_henced_microlab_value = params_ui$cost_vals$cost_henced_microlab,
        cost_CefGen_value = params_ui$cost_vals$cost_CefGen,
        cost_Merope_value = params_ui$cost_vals$cost_Merope,
        cost_oralatb_value = params_ui$cost_vals$cost_oralatb,
        cost_cultpos_value = params_ui$cost_vals$cost_cultpos,
        cost_cultneg_value = params_ui$cost_vals$cost_cultneg,
        cost_careunchan_value = params_ui$cost_vals$cost_careunchan,
        cost_caredetrio_value = params_ui$cost_vals$cost_caredetrio,
        cost_careimprov_value = params_ui$cost_vals$cost_careimprov,
        num_bottle = params_ui$num_bottle,
        N = params_ui$N,
        year_lost = params_ui$year_lost,
        daly_weight = params_ui$daly_weight,
        length_empiric = params_ui$length_empiric,
        length_targeted = params_ui$length_targeted
      )
      
      temp <- tryCatch({
        as.data.frame(t(dec_tree(params)))
      }, error = function(e) {
        # return a row of NAs if single iteration fails
        warning("dec_tree failed at iteration ", i, " -> returning NA row. Error: ", e$message)
        return(as.data.frame(matrix(NA, nrow = 1, ncol = 10)))
      })
      
      # ensure columns are named consistently
      if (ncol(temp) == 10) {
        colnames(temp) <- c("diff_mort", "diff_cost", "diff_daly", "icer",
                            "mort.noactivelab_CefGen", "mort.activelab_CefGen",
                            "cost.noactivelab_CefGen", "cost.activelab_CefGen",
                            "daly.noactivelab_CefGen", "daly.activelab_CefGen")
      } else {
        tmp <- as.data.frame(matrix(NA, nrow = 1, ncol = 10))
        colnames(tmp) <- c("diff_mort", "diff_cost", "diff_daly", "icer",
                           "mort.noactivelab_CefGen", "mort.activelab_CefGen",
                           "cost.noactivelab_CefGen", "cost.activelab_CefGen",
                           "daly.noactivelab_CefGen", "daly.activelab_CefGen")
        temp <- tmp
      }
      
      out_list[[i]] <- temp
    } # end loop
    
    res <- bind_rows(out_list)
    res
  })
  
  # CE-plane
  # CE-plane with WTP line + shaded CE region + annotated CE probability
  output$psa_ceplane <- renderPlot({
    df <- results_prob()
    validate(need(nrow(df) > 0, "No PSA results yet."))
    
    # get WTP from UI (USD per DALY)
    wtp <- req(input$wtp)
    
    # compute CE probability at this WTP
    ce_bool <- (wtp * df$diff_daly - df$diff_cost) > 0
    ce_prob <- 100 * mean(ce_bool, na.rm = TRUE)
    
    # x and y limits (pad slightly for aesthetics)
    x_min <- min(df$diff_daly, na.rm = TRUE)
    x_max <- max(df$diff_daly, na.rm = TRUE)
    x_pad <- (x_max - x_min) * 0.08
    if (x_pad == 0) x_pad <- max(1, abs(x_max)*0.08)
    x_seq <- seq(x_min - x_pad, x_max + x_pad, length.out = 200)
    
    y_min <- min(df$diff_cost, na.rm = TRUE)
    y_max <- max(df$diff_cost, na.rm = TRUE)
    y_pad <- (y_max - y_min) * 0.12
    if (y_pad == 0) y_pad <- max(1, abs(y_max)*0.12)
    
    # Build a shading polygon for the "cost-effective" region (cost < wtp * daly)
    # shading will go from x_seq left->right, with y = wtp * x, down to y_bottom
    y_bottom <- y_min - y_pad
    shade_df <- data.frame(
      x = c(x_seq, rev(x_seq)),
      y = c(wtp * x_seq, rep(y_bottom, length(x_seq)))
    )
    
    ggplot(df, aes(x = diff_daly, y = diff_cost)) +
      # shaded CE region
      geom_polygon(data = shade_df, aes(x = x, y = y),
                   inherit.aes = FALSE, fill = "lightgreen", alpha = 0.15) +
      # PSA scatter
      geom_point(size = 1.5, alpha = 0.6) +
      # WTP line through origin
      geom_abline(slope = wtp, intercept = 0, linetype = "dashed", size = 0.8, colour = "darkblue") +
      # horizontal & vertical reference lines
      geom_hline(yintercept = 0, colour = "darkgrey") +
      geom_vline(xintercept = 0, colour = "darkgrey") +
      labs(
        x = "DALYs averted (No active lab − Active lab)",
        y = "Incremental cost of maintaining active lab (USD)",
        # title = paste0("Cost-effectiveness plane (PSA) — WTP = ", scales::dollar(wtp), " /DALY"),
        # subtitle = paste0("Prob(cost-effective) = ", sprintf("%.1f%%", ce_prob))
      ) +
      coord_cartesian(xlim = c(x_min - x_pad, x_max + x_pad),
                      ylim = c(y_bottom, y_max + y_pad)) +
      theme_minimal(base_size = 14)
  })
  
  # Cost histogram
  output$psa_cost_hist <- renderPlot({
    df <- results_prob()
    validate(need(nrow(df) > 0, "No PSA results yet."))
    ggplot(df, aes(x = diff_cost)) +
      geom_histogram(bins = 60, alpha = 0.8) +
      labs(x = "Incremental cost: Active lab − No active lab", y = "Frequency", title = "Distribution of incremental cost") +
      theme_minimal(base_size = 13)
  })
  
  # DALY histogram
  output$psa_daly_hist <- renderPlot({
    df <- results_prob()
    validate(need(nrow(df) > 0, "No PSA results yet."))
    ggplot(df, aes(x = diff_daly)) +
      geom_histogram(bins = 60, alpha = 0.8) +
      labs(x = "DALYs averted: No active lab − Active lab", y = "Frequency", title = "Distribution of DALY averted") +
      theme_minimal(base_size = 13)
  })
  
  # quick summary
  output$psa_summary <- renderPrint({
    df <- results_prob()
    req(nrow(df) > 0)
    cat("PSA iterations:", nrow(df), "\n")
    cat("Median incremental cost: ", scales::dollar(median(df$diff_cost, na.rm = TRUE)), "\n")
    cat("Median incremental DALYs (No − Active): ", round(median(df$diff_daly, na.rm = TRUE), 2), "\n")
    # median ICER robust to inf/NA
    med_icer <- median(df$icer[is.finite(df$icer)], na.rm = TRUE)
    cat("Median ICER (cost / DALY): ", ifelse(is.na(med_icer), "NA", scales::dollar(med_icer)), "\n")
    # Probability CE at WTP
    wtp <- ui_params()$wtp
    ce_prob <- mean((wtp * df$diff_daly - df$diff_cost) > 0, na.rm = TRUE)
    cat("Probability active lab is cost-effective at WTP = ", scales::dollar(wtp), "/DALY: ",
        sprintf("%.1f%%", 100 * ce_prob), "\n", sep = "")
  })
  
  # Summary table showing median, Q1, Q3 (publication-style)
  output$psa_table <- renderDT({
    df <- results_prob()
    req(nrow(df) > 0)
    
    metrics <- c("Incremental cost (USD)" = "diff_cost",
                 "Incremental DALYs (No − Active)" = "diff_daly",
                 "Deaths averted (No − Active)" = "diff_im",
                 "ICER (USD per DALY)" = "icer")
    
    summary_df <- lapply(metrics, function(var) {
      x <- df[[var]]
      x <- x[!is.infinite(x)]
      med <- median(x, na.rm = TRUE)
      q25 <- quantile(x, 0.25, na.rm = TRUE)
      q75 <- quantile(x, 0.75, na.rm = TRUE)
      c(median = med, q25 = q25, q75 = q75)
    })
    
    summary_df <- do.call(rbind, summary_df) %>%
      as.data.frame() %>%
      mutate(Metric = rownames(.)) %>%
      select(Metric, median, q25, q75)
    
    # format numbers: currency for cost/ICER, numeric for DALYs/deaths
    summary_df$median <- ifelse(grepl("cost|ICER", summary_df$Metric),
                                scales::dollar(summary_df$median),
                                signif(summary_df$median, 3))
    summary_df$q25 <- ifelse(grepl("cost|ICER", summary_df$Metric),
                             scales::dollar(summary_df$q25),
                             signif(summary_df$q25, 3))
    summary_df$q75 <- ifelse(grepl("cost|ICER", summary_df$Metric),
                             scales::dollar(summary_df$q75),
                             signif(summary_df$q75, 3))
    
    datatable(summary_df, rownames = FALSE, options = list(dom = 't', pageLength = 10))
  })
  
  # Download raw PSA
  output$download_psa <- downloadHandler(
    filename = function() paste0("psa_results_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- results_prob()
      write.csv(df, file, row.names = FALSE)
    }
  )
}

# ---------- launch ----------
shinyApp(ui, server)
