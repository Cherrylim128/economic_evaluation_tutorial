#'Anti-microbial resistance in LMIC- practical applications in R
#'Ben S Cooper, Cherry Lim
#'
#'Script for decision tree model
#'This script contains the R code to build a simple decision tree model
#'that aims to evaluate the cost-effectiveness of maintaining an active microbiology
#'laboratory compared to no microbiology service.
dec_tree <- function(params){
  with(
    as.list(params), 
    {
      #/////// Pathway for active microbiology lab service arm /////// # ----
      ### Expected probabilities for each pathway
      # infection caused by bacteria
      #- culture positive
      ep1 <- p.bact * p.bact_culture_pos * p.CefGen_culpos_stepdown
      ep2 <- p.bact * p.bact_culture_pos * p.CefGen_culpos_stepup
      ep3 <- p.bact * p.bact_culture_pos * p.CefGen_culpos_nochange
      #- culture negative
      ep4 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_deter * p.CefGen_culneg_deter_stepup
      ep5 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_deter * p.CefGen_culneg_deter_nochange
      ep6 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_impro * p.CefGen_culneg_impro_stepdown
      ep7 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_impro * p.CefGen_culneg_impro_nochange
      ep8 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_nocha * p.CefGen_culneg_nocha_stepup
      ep9 <- p.bact * (1-p.bact_culture_pos) * p.CefGen_bact_culneg_nocha * p.CefGen_culneg_nocha_nochange
      checkpoint_bact_cefgen = ep1+ep2+ep3+ep4+ep5+ep6+ep7+ep8+ep9
      
      # not bacterial infection
      #- culture positive
      ep10 <- p.nobact * p.nobact_culture_pos * p.CefGen_culpos_stepdown
      ep11 <- p.nobact * p.nobact_culture_pos * p.CefGen_culpos_stepup
      ep12 <- p.nobact * p.nobact_culture_pos * p.CefGen_culpos_nochange
      #- culture negative
      ep13 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_deter * p.culneg_deter_stepup
      ep14 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_deter * p.culneg_deter_nochange
      ep15 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_impro * p.culneg_impro_stepdown
      ep16 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_impro * p.culneg_impro_nochange
      ep17 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_nocha * p.culneg_nocha_stepup
      ep18 <- p.nobact * (1-p.nobact_culture_pos) * p.CefGen_nobact_culneg_nocha * p.culneg_nocha_nochange
      checkpoint_nobact_cefgen = ep10+ep11+ep12+ep13+ep14+ep15+ep16+ep17+ep18
      
      # To check that all the probabilities should add up to 1
      activelab_CefGen_ptotal <- ep1+ep2+ep3+ep4+ep5+ep6+
        ep7+ep8+ep9+ep10+ep11+ep12+
        ep13+ep14+ep15+ep16+ep17+ep18
      
      #$$ Total costs for each pathway $$# ####
      ### Total costs for active microbiology lab service arm
      #$ ceftriaxone/gentamicin ####
      tc1 <- cost_cultpos + cost_oralatb*length_targeted
      tc2 <- cost_cultpos + cost_Merope*length_targeted                 
      tc3 <- cost_cultpos + cost_CefGen*length_targeted                    
      tc4 <- cost_cultneg + (cost_caredetrio + cost_Merope)*length_targeted                
      tc5 <- cost_cultneg + (cost_caredetrio + cost_CefGen)*length_targeted                 
      tc6 <- cost_cultneg + (cost_careimprov + cost_oralatb)*length_targeted                
      tc7 <- cost_cultneg + (cost_careimprov + cost_CefGen)*length_targeted                
      tc8 <- cost_cultneg + (cost_careunchan + cost_Merope)*length_targeted                
      tc9 <- cost_cultneg + (cost_careunchan + cost_CefGen)*length_targeted                

      #$ Expected Total Costs in active microbiology lab service arm ----
      etc.activelab_CefGen <- ((ep1 * tc1) + (ep2 * tc2) + (ep3 * tc3) + 
                                (ep4 * tc4) + (ep5 * tc5) + 
                                (ep6 * tc6) + (ep7 * tc7) + 
                                (ep8 * tc8) + (ep9 * tc9) + 
                                (ep10 * tc1) + (ep11 * tc2) + (ep12 * tc3) + 
                                (ep13 * tc4) + (ep14 * tc5) + 
                                (ep15 * tc6) + (ep16 * tc7) + 
                                (ep17 * tc8) + (ep18 * tc9) + 
                                (cost_CefGen * length_empiric) + cost_henced_microlab)*N
      
      #^^ Outcome (mortality) //# ----
      mort.activelab_CefGen <- ((ep1 * mort_bact_culpos_stepdown) + (ep2 * mort_bact_culpos_stepup) + (ep3 * mort_bact_culpos_nochange) + 
                                 (ep4 * mort_bact_culneg_stepup_deter) + (ep5 * mort_bact_culneg_nochange_deter) + 
                                 (ep6 * mort_bact_culneg_stepdown_impro) + (ep7 * mort_bact_culneg_nochange_impro) + 
                                 (ep8 * mort_bact_culneg_stepup_uncha) + (ep9 * mort_bact_culneg_nochange_uncha) + 
                                 (ep10 * mort_nobact_culpos_stepdown) + (ep11 * mort_nobact_culpos_stepup) + (ep12 * mort_nobact_culpos_nochange) + 
                                 (ep13 * mort_nobact_culneg_stepup_deter) + (ep14 * mort_nobact_culneg_nochange_deter) + 
                                 (ep15 * mort_nobact_culneg_stepdown_impro) + (ep16 * mort_nobact_culneg_nochange_impro) + 
                                 (ep17 * mort_nobact_culneg_stepup_uncha) + (ep18 * mort_nobact_culneg_nochange_uncha) ) * N
      
      #^ Expected effectiveness (DALY) ----
      daly.activelab_CefGen <- 
        ((ep1 * (yll_bact_culpos_stepdown+ yld_bact_culpos_stepdown)) + (ep2 * (yll_bact_culpos_stepup+ yld_bact_culpos_stepup)) + (ep3 * (yll_bact_culpos_nochange+ yld_bact_culpos_nochange)) + 
           (ep4 * (yll_bact_culneg_stepup_deter+ yld_bact_culneg_stepup_deter)) + (ep5 * (yll_bact_culneg_nochange_deter+ yld_bact_culneg_nochange_deter)) + 
           (ep6 * (yll_bact_culneg_stepdown_impro+ yld_bact_culneg_stepdown_impro)) + (ep7 * (yll_bact_culneg_nochange_impro+ yld_bact_culneg_nochange_impro)) + 
           (ep8 * (yll_bact_culneg_stepup_uncha+ yld_bact_culneg_stepup_uncha)) + (ep9 * (yll_bact_culneg_nochange_uncha+ yld_bact_culneg_nochange_uncha)) + 
           (ep10 * (yll_nobact_culpos_stepdown+ yld_nobact_culpos_stepdown)) + (ep11 * (yll_nobact_culpos_stepup+ yld_nobact_culpos_stepup)) + (ep12 * (yll_nobact_culpos_nochange+ yld_nobact_culpos_nochange)) + 
           (ep13 * (yll_nobact_culneg_stepup_deter+ yld_nobact_culneg_stepup_deter)) + (ep14 * (yll_nobact_culneg_nochange_deter+ yld_nobact_culneg_nochange_deter)) + 
           (ep15 * (yll_nobact_culneg_stepdown_impro+ yld_nobact_culneg_stepdown_impro)) + (ep16 * (yll_nobact_culneg_nochange_impro+ yld_nobact_culneg_nochange_impro)) + 
           (ep17 * (yll_nobact_culneg_stepup_uncha+ yld_nobact_culneg_stepup_uncha)) + (ep18 * (yll_nobact_culneg_nochange_uncha+ yld_nobact_culneg_nochange_uncha)))*N
  
      
      #/////// Pathways for no microbiology laboratory services arm /////// # ----
      #/ ceftriaxone/gentamicin ----
      # Bacterial infection
      np1 <- p.bact * p.CefGen_bact_nocul_deter * p.CefGen_nocul_deter_stepup
      np2 <- p.bact * p.CefGen_bact_nocul_deter * p.CefGen_nocul_deter_nochange
      np3 <- p.bact * p.CefGen_bact_nocul_impro * p.CefGen_nocul_impro_stepdown
      np4 <- p.bact * p.CefGen_bact_nocul_impro * p.CefGen_nocul_impro_nochange
      np5 <- p.bact * p.CefGen_bact_nocul_nocha * p.CefGen_nocul_nocha_stepup
      np6 <- p.bact * p.CefGen_bact_nocul_nocha * p.CefGen_nocul_nocha_nochange
      checkpoint2_bact_cefgen <- np1 + np2 + np3 + np4 + np5 + np6
      # Non-bacterial infection
      np7 <- p.nobact * p.CefGen_nobact_nocul_deter * p.CefGen_nocul_deter_stepup
      np8 <- p.nobact * p.CefGen_nobact_nocul_deter * p.CefGen_nocul_deter_nochange
      np9 <- p.nobact * p.CefGen_nobact_nocul_impro * p.CefGen_nocul_impro_stepdown
      np10 <- p.nobact * p.CefGen_nobact_nocul_impro * p.CefGen_nocul_impro_nochange
      np11 <- p.nobact * p.CefGen_nobact_nocul_nocha * p.CefGen_nocul_nocha_stepup
      np12 <- p.nobact * p.CefGen_nobact_nocul_nocha * p.CefGen_nocul_nocha_nochange
      checkpoint2_nobact_cefgen <- np7 + np8 + np9 + np10 + np11 + np12
      
      # To confirm probabilities add up to 1
      noactivelab_CefGen_ptotal <- np1+np2+np3+np4+np5+np6+
        np7+np8+np9+np10+np11+np12
      
      #$$ Total costs for each pathway $$# ####
      ### Total costs for no microbiology lab arm
      #$ ceftriaxone/gentamicin ----
      nstc1 <- (cost_caredetrio + cost_Merope)*length_targeted            
      nstc2 <- (cost_caredetrio + cost_CefGen)*length_targeted            
      nstc3 <- (cost_careimprov + cost_oralatb)*length_targeted              
      nstc4 <- (cost_careimprov + cost_CefGen)*length_targeted             
      nstc5 <- (cost_careunchan + cost_Merope)*length_targeted                             
      nstc6 <- (cost_careunchan + cost_CefGen)*length_targeted            
      
      #$ Expected Total Costs by empirical treatment choice in enhanced surveillance ----
      etc.noactivelab_CefGen <- 
        ((np1 * nstc1) + (np2 * nstc2) + 
           (np3 * nstc3) + (np4 * nstc4) + 
           (np5 * nstc5) + (np6 * nstc6) + 
           (np7 * nstc1) + (np8 * nstc2) + 
           (np9 * nstc3) + (np10 * nstc4) + 
           (np11 * nstc5) + (np12 * nstc6) + 
           cost_CefGen*length_empiric)*N

      #^^ Outcome (mortality) by empirical treatment choice //# ----
      mort.noactivelab_CefGen <- 
        ((np1 * mort_bact_nocul_stepup_deter) + (np2 * mort_bact_nocul_nochange_deter) + 
           (np3 * mort_bact_nocul_stepdown_impro) + (np4 * mort_bact_nocul_nochange_impro) + 
           (np5 * mort_bact_nocul_stepup_uncha) + (np6 * mort_bact_nocul_nochange_uncha) + 
           (np7 * mort_nobact_nocul_stepup_deter) + (np8 * mort_nobact_nocul_nochange_deter) + 
           (np9 * mort_nobact_nocul_stepdown_impro) + (np10 * mort_nobact_nocul_nochange_impro) + 
           (np11 * mort_nobact_nocul_stepup_uncha) + (np12 * mort_nobact_nocul_nochange_uncha))*N 
      
      #^ Expected effectiveness (DALY) by empirical treatment choice ----
      daly.noactivelab_CefGen <- 
        ((np1 * (yll_bact_nocul_stepup_deter+ yld_bact_nocul_stepup_deter)) + (np2 * (yll_bact_nocul_nochange_deter+ yld_bact_nocul_nochange_deter)) + 
           (np3 * (yll_bact_nocul_stepdown_impro+ yld_bact_nocul_stepdown_impro)) + (np4 * (yll_bact_nocul_nochange_impro+ yld_bact_nocul_nochange_impro)) + 
           (np5 * (yll_bact_nocul_stepup_uncha+ yld_bact_nocul_stepup_uncha)) + (np6 * (yll_bact_nocul_nochange_uncha+ yld_bact_nocul_nochange_uncha)) + 
           (np7 * (yll_nobact_nocul_stepup_deter+ yld_nobact_nocul_stepup_deter)) + (np8 * (yll_nobact_nocul_nochange_deter+ yld_nobact_nocul_nochange_deter)) + 
           (np9 * (yll_nobact_nocul_stepdown_impro+ yld_nobact_nocul_stepdown_impro)) + (np10 * (yll_nobact_nocul_nochange_impro+ yld_nobact_nocul_nochange_impro)) + 
           (np11 * (yll_nobact_nocul_stepup_uncha+ yld_nobact_nocul_stepup_uncha)) + (np12 * (yll_nobact_nocul_nochange_uncha+ yld_nobact_nocul_nochange_uncha)))*N
      
      #*** ICER ***# ----
      ## Overall for the two types of surveillance
      im <- mort.noactivelab_CefGen-mort.activelab_CefGen
      ic <- etc.activelab_CefGen-etc.noactivelab_CefGen
      ie <- daly.noactivelab_CefGen-daly.activelab_CefGen
      icer <- ic/ie
      
      #### Printings ----
      print(c(mort.noactivelab_CefGen, mort.activelab_CefGen))
      print(c(etc.noactivelab_CefGen, etc.activelab_CefGen))
      print(c(daly.noactivelab_CefGen, daly.activelab_CefGen))
      #print(c(checkpoint_bact_cefgen, checkpoint_nobact_cefgen))
      #print(c(checkpoint2_bact_cefgen, checkpoint2_nobact_cefgen))
      ### list to return ----
      return(c(im, ic, ie, icer, 
               mort.noactivelab_CefGen, mort.activelab_CefGen,
               etc.noactivelab_CefGen, etc.activelab_CefGen,
               daly.noactivelab_CefGen, daly.activelab_CefGen))
    }
  )
}
