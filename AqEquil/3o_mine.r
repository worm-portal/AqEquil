library(dplyr)
library(CHNOSZ)
library(stringr)


### helper functions

# trims away leading and trailing spaces and condenses multiple spaces between words
trimspace <- function(str){
    gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", str, perl=TRUE)
}

# isolate a substring by trimming off the portions before and after it
isolate_block <- function(str, begin_str, end_str){
  return(sub(end_str, "", sub(begin_str, "", str)))
}

### main functions

mine_3o <- function(this_file,
                    rxn_table,
                    get_aq_dist=T,
                    get_mass_contribution=T,
                    get_mineral_sat=T,
                    get_redox=T,
                    get_charge_balance=T,
                    get_ion_activity_ratios=T,
                    get_fugacity=T,
                    get_affinity_energy=T,
                    not_limiting=c("H+", "OH-", "H2O"),
                    mass_contribution_other=T,
                    verbose=1){
    
  # set directory to rxn_3o folder where .3o files are kept
  setwd("rxn_3o")
    
  # read .3o file as a string
  fileName <- this_file
  extractme <- readChar(fileName, file.info(fileName)$size) # readChar requires filesize, which is obtained from file.info()$size

  # get sample name
  this_name <- trimspace(isolate_block(extractme, begin_str="^.*\\|Sample:\\s+", end_str="\\|\\n\\|.*$"))

  if(verbose > 1){
    writeLines(paste0("Processing EQ3 output for ", this_name))
  }

  # check if file experienced errors. If so, skip processing the file:
  if (grepl("Normal exit", extractme) == FALSE | grepl("\\* Error", extractme)){
    if(verbose > 0){
      writeLines(paste0("\nSample ", this_name, " experienced errors during speciation:"))

      output_error <- str_extract_all(extractme, regex("\\* Error.*?\n(\n|$)", dotall=T))
        
      output_error <- lapply(output_error, function(x) gsub("\n", " ", x))
      output_error <- lapply(output_error, str_squish)
        
      for(i in 1:length(output_error)){
        writeLines(output_error[[i]])
      }
      
    }

    setwd("../")
    return(list())
  }

  sample_3o <- list()

  sample_3o[["filename"]] <- this_file
  sample_3o[["name"]] <- this_name
  
  ### Begin mining temperature, pressure, water properties
    
  # mine params
  sample_3o[["temperature"]] <- isolate_block(str=extractme, begin_str="^.*Temperature=\\s+", end_str="\\s+.*$")
  sample_3o[["pressure"]] <- isolate_block(str=extractme, begin_str="^.*Pressure=\\s+", end_str="\\s+.*$")
  sample_3o[["logact_H2O"]] <- isolate_block(str=extractme, begin_str="^.*Log activity of water=\\s+", end_str="\\s+.*$")
  sample_3o[["H2O_density"]] <- isolate_block(str=extractme, begin_str="^.*Solution density =\\s+", end_str="\\s+.*$")
  sample_3o[["H2O_molality"]] <- 55.348/as.numeric(sample_3o[["H2O_density"]])
  sample_3o[["H2O_log_molality"]] <- log10(sample_3o[["H2O_molality"]])

  ### Begin extracting 'Distribution of Aqueous Solute Species' ###
  if (get_aq_dist){
    # string to isolate the aqueous species distribution section:
    front_trim <- "^.*\n\n\n\n                --- Distribution of Aqueous Solute Species ---\n\n    Species                  Molality    Log Molality   Log Gamma  Log Activity\n\n\\s+"

    # isolate species distribution block
    species_block <- isolate_block(str=extractme, begin_str=front_trim, end_str="\n\n.*$")

    # split into substrings, each representing a separate row in the table
    species_block <- strsplit(species_block, "\n")

    #create an empty data frame to store results
    df <- data.frame(species = character(0),
                      molality = character(0),
                      log_molality = character(0),
                      log_gamma = character(0),
                      log_activity = character(0),
                      stringsAsFactors = FALSE)

    # convert into dataframe
    for(this_row in species_block[[1]]){

      # mine row data
      this_row <- trimspace(this_row)
      this_row_data <- strsplit(this_row, " ")[[1]]

      # create a dataframe with results
      this_df <- data.frame(species=this_row_data[1],
                            molality=this_row_data[2],
                            log_molality=this_row_data[3],
                            log_gamma=this_row_data[4],
                            log_activity=this_row_data[5],
                            stringsAsFactors=FALSE)

      # bind distribution of aqueous species to dataframe
      df <- rbind(df, this_df)

    }

    if(!("H2O" %in% df$species)){
      # add a row for water
      df <- rbind(df, data.frame(species="H2O",
                      molality=sample_3o[["H2O_molality"]],
                      log_molality=sample_3o[["H2O_log_molality"]],
                      log_gamma=1,
                      log_activity=sample_3o[["logact_H2O"]],
                      stringsAsFactors=FALSE))
    }

    # set rownames of aqueous species block as species names
    rownames(df) <- df$species
    df$species <- NULL

    # add aqueous block to this sample data
    sample_3o[["aq_distribution"]] <- df
      
  } # end of 'aqueous distribution' extraction
    

  if(get_mass_contribution){
    ### begin extracting 'Major Species by Contribution to Aqueous Mass Balances' ###
    
    # string to isolate the species saturation section:
    front_trim <- "^.*\n\n\n      --- Major Species by Contribution to Aqueous Mass Balances ---\n\n\n"

    # isolate contribution block
    contrib_block <- isolate_block(str=extractme, begin_str=front_trim, end_str="\n\n\n\n.*$")
    
    # split into substrings, each representing a separate row in the table
    contrib_block <- strsplit(contrib_block, "\n")
    # remove blank linkes
    contrib_block[[1]] <- contrib_block[[1]][contrib_block[[1]] != ""]
  
    # loop through rows in this block and mine contributions
    mine_vals <- FALSE
    for(this_row in contrib_block[[1]]){
      if(grepl("Accounting for", this_row)){
        # get basis species for this block
        this_basis <- sub(" Species Accounting for 99% or More of Aqueous ", "", this_row)
      } else if (grepl("Per Cent", this_row)){
        # get ready to mine data for this basis species
        mine_vals <- TRUE
        df_basis <- data.frame(species=character(0),
                               molality=character(0),
                               factor=character(0),
                               percent=character(0),
                               stringsAsFactors = FALSE)
      } else if (mine_vals && !grepl(" - - - - - - - - -", this_row)){
        # mine data from this row
        row_data <- trimspace(this_row)
        row_data <- strsplit(row_data, " ")[[1]]
        df_basis <- rbind(df_basis, data.frame(species=row_data[1], factor=row_data[2],
                        molality=row_data[3], percent=row_data[4], stringsAsFactors = FALSE))
      } else if (grepl(" - - - - - - - - -", this_row)){
        # stop mining for this basis species
        mine_vals <- FALSE
        # specify rownames for this contribution block
        rownames(df_basis) <- df_basis$species
        df_basis$species <- NULL
        # add contribution data to list of sample data
        sample_3o[["mass_contribution"]][[this_basis]] <- df_basis
      }
    }
  } # end 'aqueous contribution' extraction


  ### Begin mining mineral saturation section  
  if(get_mineral_sat){
      
    # string to isolate the mineral saturation section:
    front_trim <- "^.*\n\n\n\n           --- Saturation States of Pure Solids ---\n\n       Phase                      Log Q/K    Affinity, kcal\n\n\\s+"

    # isolate mineral block
    mineral_block <- isolate_block(str=extractme, begin_str=front_trim, end_str="\n\n.*$")

    # split into substrings, each representing a separate row in the table
    mineral_block <- strsplit(mineral_block, "\n")

    # create an empty data frame to store results
    df <- data.frame(mineral = character(0),
                     logQoverK = character(0),
                     affinity = character(0),
                     stringsAsFactors = FALSE)

    # convert into dataframe
    for(this_row in mineral_block[[1]]){
      # get row data
      this_row_data <- strsplit(trimspace(this_row), " ")[[1]]
      
      # create a dataframe with mined data
      this_df <- data.frame(mineral = this_row_data[1],
                            logQoverK = this_row_data[2],
                            affinity = this_row_data[3],
                            stringsAsFactors = FALSE)
      
      # bind results to dataframe
      df <- rbind(df, this_df)
    }
    rownames(df) <- df$mineral
    df$mineral <- NULL
      
    # add mineral saturation block to this sample data
    sample_3o[["mineral_sat"]] <- df
  } # end 'mineral saturation affinity' extraction
    

  ### Begin mining redox data
  if(get_redox){
    # string to isolate the redox section:
    front_trim <- "^.*\n\n\n\n                --- Aqueous Redox Reactions ---\n\n   Couple                           Eh, volts      pe-      log fO2   Ah, kcal\n\n\\s+"

    # isolate redox block
    redox_block <- isolate_block(str=extractme, begin_str=front_trim, end_str="\n\n.*$")

    # split into substrings, each representing a separate row in the table
    redox_block <- strsplit(redox_block, "\n")

    #create an empty data frame to store results
    df <- data.frame(couple = character(0),
                     Eh = character(0),
                     pe = character(0),
                     logfO2 = character(0),
                     Ah = character(0),
                     stringsAsFactors = FALSE)

    # convert into dataframe
    for(this_row in redox_block[[1]]){
      
      # get row data
      this_row_data <- strsplit(trimspace(this_row), " ")[[1]]

      # create a dataframe with results
      this_df <- data.frame(couple = this_row_data[1],
                            Eh = this_row_data[2],
                            pe = this_row_data[3],
                            logfO2 = this_row_data[4],
                            Ah = this_row_data[5],
                            stringsAsFactors = FALSE)
      
      # bind results to dataframe
      df <- rbind(df, this_df)
    }
    rownames(df) <- df$couple
    df$couple <- NULL
    
    # add mineral saturation block to this sample data
    sample_3o[["redox"]] <- df
  } # end redox extraction
    

  ### begin mining charge balance data
  if(get_charge_balance){
    # string to isolate ionic strength:
    front_trim <- "^.*Ionic strength \\(I\\)=\\s+"

    # isolate ionic strength
    IS <- isolate_block(str=extractme, begin_str=front_trim, end_str="\\s+.*$")
    names(IS) <- "ionic strength"

    # string to isolate stoichiometric ionic strength:
    front_trim <- "^.*Stoichiometric ionic strength=\\s+"

    IS_stoich <- isolate_block(str=extractme, begin_str=front_trim, end_str="\\s+.*$")
    names(IS_stoich) <- "stoichiometric ionic strength"

    # string to isolate the electrical balance section:
    front_trim <- "^.*Sigma\\(mz\\) cations=\\s+"

    elec_block <- isolate_block(str=extractme, begin_str=front_trim, end_str="\n\n.*$")

    # split electrical block into strings and numerics
    elec_block <- strsplit(elec_block, "=\\s+|\n\\s+")[[1]]

    elec_block <- c("sigma(mz) cations"=elec_block[1],
                    "sigma(mz) anions"=elec_block[3],
                    "total charge"=elec_block[5],
                    "mean charge"=elec_block[7],
                    "charge imbalance"=elec_block[9])

    # string to isolate charge balance:
    front_trim <- "^.*The electrical imbalance is:\n\n\\s+"

    cbal_bal <- isolate_block(str=extractme, begin_str=front_trim, end_str="\n\n.*$")

    # split electrical block into strings and numerics
    cbal_block <- strsplit(cbal_bal, " per cent|\n\\s+")[[1]]

    cbal_block <- c("charge imbalance % of total charge"=cbal_block[1],
                    "charge imbalance % of mean charge"=cbal_block[3])
  
    sample_3o[["charge_balance"]] <- c(IS, IS_stoich, elec_block, cbal_block)
  } # end charge balance extraction


  if(get_ion_activity_ratios){
    ion_ratio_block <- isolate_block(extractme, "^.*--- Ion-H\\+ Activity Ratios ---\n\n", "\n\n.*$")
    ion_ratio_block_split <- strsplit(ion_ratio_block, "\n")[[1]]
    ion_ratio_block_split <- strsplit(ion_ratio_block_split, "=")
    ion_ratio_logs <- trimspace(lapply(ion_ratio_block_split, `[[`, 1))
    ion_ratio_values <- suppressWarnings(as.numeric(lapply(ion_ratio_block_split, `[[`, 2)))
      
    which_to_divide <- grepl("/", ion_ratio_logs) # which of these ratios divide by H+ (instead of multiply?)
    hydrogen_exponents <- unlist(lapply(lapply(strsplit(gsub("^Log \\( a\\(", "", gsub(" \\)$", "", ion_ratio_logs)), "\\)xx "), `[[`, 2), as.numeric))
    ion <- unlist(lapply(strsplit(gsub("^Log \\( a\\(", "", gsub(" \\)$", "", ion_ratio_logs)), "\\) [x|/] a\\("), `[[`, 1))
    ion_times <- lapply(strsplit(gsub("^Log \\( a\\(", "", gsub(" \\)$", "", ion_ratio_logs)), "\\) x a\\(")[!which_to_divide], `[[`, 1)
    ion_divide <- lapply(strsplit(gsub("^Log \\( a\\(", "", gsub(" \\)$", "", ion_ratio_logs)), "\\) / a\\(")[which_to_divide], `[[`, 1)
  
    names(ion_ratio_values)[!which_to_divide] <- paste0("Log(", ion_times, " x H+**", hydrogen_exponents[!which_to_divide], ")")
    names(ion_ratio_values)[which_to_divide]  <- paste0("Log(", ion_divide, " / H+**", hydrogen_exponents[which_to_divide], ")")
    
    ion_ratio_values[is.na(ion_ratio_values)] <- "NA"

    df <- data.frame("values"=ion_ratio_values,
               "H_exponent"=hydrogen_exponents,
               "divide"=which_to_divide,
               "ion"=ion)

    df <- transform(transform(df, values = as.character(values)), values=as.numeric(values))
    df <- transform(df, ion = as.character(ion))
    sample_3o[["ion_activity_ratios"]] <- df
  }

  ### begin fugacity mining
  if(get_fugacity){
    fugacity_block <- isolate_block(extractme, "^.*--- Fugacities ---\n\n", "\n\n\n.*$")
    str <- str_squish(strsplit(fugacity_block, "\n")[[1]])
    str <- str[3:length(str)]
    split_str <- strsplit(str, " ")
    df <- data.frame(gas=unlist(lapply(split_str, `[[`, 1)),
                     log_fugacity=as.numeric(unlist(lapply(split_str, `[[`, 2))),
                     stringsAsFactors=F, row.names=1)
    sample_3o[["fugacity"]] <- df
  }
    
  ### begin energy mining
  if(get_affinity_energy){
    this_temp <- as.numeric(sample_3o[["temperature"]])
    this_pres <- as.numeric(sample_3o[["pressure"]])
    this_logact_H2O <- as.numeric(sample_3o[["logact_H2O"]])
    this_H2O_log_molality <- as.numeric(sample_3o[["H2O_log_molality"]])

    df <- sample_3o[["aq_distribution"]]
    
    # get a list of mineral names in CHNOSZ
    CHNOSZ_cr_names <- unlist(thermo()$OBIGT %>% filter(state == "cr") %>% select(name))

    # clear any molal values upon moving to this sample
    remaining_react_molal_with_product_molal <- c()
    other_reactants_and_prod_names <- c()
    other_reactants_and_prod <- c()
    
    # create a dataframe for storing results
    df_rxn <- data.frame(rxn = character(0),
                         affinity = numeric(0),
                         energy_supply = numeric(0),
                         mol_rxn = numeric(0),
                         electrons= numeric(0),
                         reaction = character(0),
                         limiting = character(0),
                         stringsAsFactors = FALSE)
      
    for(rxn in rxn_table){
      rxn_split <- unlist(strsplit(rxn, "\t"))
      rxn_name  <- rxn_split[1]
      electrons <- as.numeric(rxn_split[2])
      full_rxn  <- rxn_split[3:length(rxn_split)]
      if(length(full_rxn) %% 2 != 0){
        stop(paste("Error: Number of reaction coefficients and species do not match in reaction", rxn_name))
      }
        
      stoichs <- as.numeric(full_rxn[c(TRUE, FALSE)])
      species_EQ3 <- full_rxn[c(FALSE, TRUE)]
#       species_CHNOSZ <- gsub(",AQ", "", species_EQ3) # this won't work for mineral names, gases, etc.!
#       species_CHNOSZ <- gsub("METHANE", "CH4", species_CHNOSZ) # temporary fix for methane
#       species_CHNOSZ <- gsub("SULFUR", "S", species_CHNOSZ) # temporary fix for sulfur
#       species_CHNOSZ <- gsub("Ca(CO3)", "CaCO3", species_CHNOSZ) # temporary fix for calcium carbonate
#       species_CHNOSZ <- gsub("Ca(CO3),AQ", "CaCO3", species_CHNOSZ) # temporary fix for calcium carbonate

      species_CHNOSZ <- species_EQ3 # assuming no differences between EQ3 and CHNOSZ naming
        
      # substitute CHNOSZ's lowercase mineral names
      for(species in species_CHNOSZ){
          
        lowercase_species <- tolower(species)
        if(lowercase_species %in% CHNOSZ_cr_names){
          #species_CHNOSZ <- gsub(species, lowercase_species, species_CHNOSZ)
          # add the mineral to master_df and master_df_mol dataframes assuming an
          # activity of 1 (log activity 0)
          df[species, "log_activity"] <- 0
          not_limiting <- c(not_limiting, species)
        }
      }
      
      ### calculate Q using EQ3-speciated activities
      # get speciated activities from master_df
      activities <- c()
      molalities <- c()

      for(species in species_EQ3){
        
        if(species %in% rownames(df)){
        
          activities <- c(activities, 10^as.numeric(df[species, "log_activity"]))

          if(!grepl("sub$", rxn_name)){
            molalities <- c(molalities, as.numeric(df[species, "molality"]))
          } else {
            molalities <- c(molalities, remaining_react_molal_with_product_molal[species])
          }

        } else {
          activities <- c(activities, NA)
          molalities <- c(molalities, NA)
        }
      }
      
      
      if(!(NA %in% activities) & !is.null(molalities) & !is.null(activities)){
        
        names(activities) <- species_EQ3
        names(molalities) <- species_EQ3
        
        if(grepl("sub$", rxn_name)){
          other_reactants_and_prod_names <- setdiff(names(remaining_react_molal_with_product_molal), names(molalities))
          other_reactants_and_prod <- remaining_react_molal_with_product_molal[other_reactants_and_prod_names]
        }

        # calculate Q
        reactant_stoich <- stoichs[which(stoichs < 0)]
        reactant_activities <- activities[which(stoichs < 0)]
        reactant_molalities <- molalities[which(stoichs < 0)]
        product_stoich <- stoichs[which(stoichs > 0)]
        product_activities <- activities[which(stoichs > 0)]
        product_molalities <- molalities[which(stoichs > 0)]
        this_logQ <- sum(abs(product_stoich)*log10(product_activities)) - sum(abs(reactant_stoich)*log10(reactant_activities))

        if(!is.na(this_logQ)){
          ### calculate K using subcrt() function in CHNOSZ
          this_logK <- suppressMessages(subcrt(species=species_CHNOSZ,
                                               coeff=stoichs,
                                               #state=phase,
                                               T=this_temp,
                                               P=this_pres)$out$logK)
        }else{
          this_logK <- NA
        }

        ### calculate activity of limiting reactant
        reactant_names <- species_EQ3[which(stoichs < 0)]
        product_names <- species_EQ3[which(stoichs > 0)]
        not_lim_index <- which(reactant_names %in% not_limiting)

        names(reactant_molalities) <- reactant_names

        if(length(not_lim_index) == 0){
          molality_div_stoich <- reactant_molalities/abs(reactant_stoich)
        }else if(all(reactant_names %in% not_limiting)){
          # if there are no limiting reactants (e.g., reactants are all minerals)
          # then append NAs and continue...
          df_rxn <- rbind(df_rxn, data.frame(rxn=rxn_name,
                                             affinity=affinity,
                                             energy_supply=NA,
                                             mol_rxn=NA,
                                             electrons=electrons,
                                             reaction=paste(full_rxn,collapse=" "),
                                             limiting=NA,
                                             stringsAsFactors=FALSE))
          next
        }else{
          molality_div_stoich <- reactant_molalities[-not_lim_index]/abs(reactant_stoich[-not_lim_index])
        }
        limiting_reactant <- min(molality_div_stoich)


        which_limiting <- which(reactant_molalities/abs(reactant_stoich) == limiting_reactant)

        # create a string of all limiting reactants
        limiting_reactants <- paste(reactant_names[which_limiting], collapse=" ")

        ### If there is more than one limiting reactant, pick the first one.
        # Prevents warnings when calculating reactant molalities when the limiting
        # reactant runs out. The math should work out the same regardless of which limiting
        # reactant is chosen.
        which_limiting <- which_limiting[1]

        # calculate moles of rxn before limiting reactant runs out
        mol_rxn <- reactant_molalities[which_limiting] / abs(reactant_stoich[which_limiting])
        
        # calculate reactant molalities that remain when the limiting reactant runs out
        remaining_reactant_molalities <- reactant_molalities - abs(reactant_stoich) * mol_rxn

        
        # if 'nonlimiting' species are specified by user, restore their molalities back to their original values.
        if(length(not_lim_index) > 0){
          remaining_reactant_molalities[not_lim_index] <- reactant_molalities[not_lim_index]
        }


        remaining_react_molal_with_product_molal <- c(remaining_reactant_molalities, product_molalities)
        names(remaining_react_molal_with_product_molal) <- c(reactant_names, product_names)

        # attach molalities of other reactants and products that are in the "mother" reaction
        if(grepl("sub$", rxn_name)){
          remaining_react_molal_with_product_molal <- c(remaining_react_molal_with_product_molal, other_reactants_and_prod)
        }

        ### calculate affinity, A
        this_A <- 0.008314*(this_temp+273.15)*2.302585*(this_logK-this_logQ) # in kJ/mol, A=RT*ln(K/Q)=RT*2.302585*(logK-logQ)

        ### calculate 'energy' in kJ/kg H2O by multiplying affinity (kJ/mol) by activity (mol/kg) of limiting reactant
        this_energy <- this_A * limiting_reactant


      } else { # if there is an NA in one of the activities
        this_A <- NA
        this_energy <- NA
        mol_rxn <- NA
        limiting_reactants <- NA
      }
        
      # unit conversion
      affinity <- ((this_A*1000)/4.184)/electrons # in cal/mol e-
      energy_supply <- (this_energy*1000)/4.184 # in cal/kg
        
      # append results
      df_rxn <- rbind(df_rxn, data.frame(rxn=rxn_name,
                                         affinity=affinity,
                                         energy_supply=energy_supply,
                                         mol_rxn=mol_rxn,
                                         electrons=electrons,
                                         reaction=paste(full_rxn,collapse=" "),
                                         limiting=limiting_reactants,
                                         stringsAsFactors=FALSE))
    } # end rxn loop
      
    rownames(df_rxn) <- df_rxn$rxn
    df_rxn$rxn <- NULL
      
    # create a dataframe for storing results
    df_rxn_sum <- data.frame(affinity = numeric(0),
                         energy_supply = numeric(0),
                         electrons = numeric(0),
                         reaction = character(0),
                         stringsAsFactors = FALSE)
      
    ### sum reaction clusters (a reaction and its subreactions)
    if(sum(grepl("_sub$", rownames(df_rxn))) > 0){
      # perform this chunk of code if "_sub" rxns are present
      rxn_list_sub <- c() # initialize vector of sub-reaction names
      df_rxn[, "mol_rxn_perc"] <- NA # add a new column to df_rxn to store percent mol rxn
      
      for(rxn in c(rownames(df_rxn), "final_energy")){
        # loop through columns (plus a dummy "final_energy" column)
          
        if(grepl("_sub", rxn)){
          # if the rxn represents a sub-reaction in the cluster, add column name to vector
          rxn_list_sub <- c(rxn_list_sub, rxn)
          
        } else {
          # if the rxn does not represent a sub-reaction...
          if(rxn != "final_energy"){
            reaction <- df_rxn[rxn, "reaction"]
          }
          
          if(length(rxn_list_sub) != 0){
            # sum the previous reaction cluster and append to dataframe
            rxn_sum_sub_E <- colSums(df_rxn[rxn_list_sub, "energy_supply", drop=FALSE])
            rxn_sum_sub_mol <- colSums(df_rxn[rxn_list_sub, "mol_rxn", drop=FALSE])
            rxn_clust_sub_A <- df_rxn[rxn_list_sub, , drop=FALSE] %>% mutate(A_weighted=affinity*(mol_rxn/rxn_sum_sub_mol))
            rxn_clust_sub_mol_perc <- df_rxn[rxn_list_sub, , drop=FALSE] %>% mutate(mol_perc=round(100*(mol_rxn/rxn_sum_sub_mol), 1))
            df_rxn[rxn_list_sub, "mol_rxn_perc"] <- rxn_clust_sub_mol_perc[, "mol_perc"]
            rxn_sum_sub_A <- colSums(rxn_clust_sub_A[, "A_weighted", drop=FALSE])
            rxn_sum_row <- data.frame(row.names=rxn_name, "affinity"=rxn_sum_sub_A, "energy_supply"=rxn_sum_sub_E, electrons=electrons, reaction=reaction)
            df_rxn_sum <- rbind(df_rxn_sum, rxn_sum_row)
          }
        
          # re-initialize vector for new reaction or reaction cluster
          rxn_name <- rxn
          rxn_list_sub <- c(rxn)
        
        }
      }
    }else{
      df_rxn_sum <- df_rxn
    }
    
    # append affinity and energy results to this sample's data
    sample_3o[["affinity_energy_raw"]] <- df_rxn
    sample_3o[["affinity_energy"]] <- df_rxn_sum
  } # end calculation of affinity and energy supply

  setwd("../")

  return(sample_3o)

}


# function to melt aqueous contribution data from multiple samples into
# a single dataframe and then return it.
melt_mass_contribution <- function(batch_3o, other=F, verbose=1){

  # initialize empty dataframe
  df_aq_cont <- data.frame(sample=character(0),
                           basis=character(0),
                           species=character(0),
                           factor=character(0),
                           molality=character(0),
                           percent=character(0),
                           stringsAsFactors=FALSE)

  # get all aqueous contribution data
  mass_contributions <- lapply(batch_3o[["sample_data"]], `[[`, 'mass_contribution')  

  # loop through each sample and basis species
  for(sample in names(mass_contributions)){
    if(verbose > 1){
      writeLines(paste0("Processing mass contribution of basis species in ", sample, "..."))
    }
    for(basis in names(mass_contributions[[sample]])){
      df <- mass_contributions[[sample]][[basis]]
      df[, "basis"] <- basis
      df[, "sample"] <- sample
      df[, "species"] <- rownames(df)
      rownames(df) <- NULL
        
      if(other){
        percent <- round(100-sum(as.numeric(df[, "percent"])), 2)
        df <- rbind(df, data.frame(sample=sample,
                                   basis=basis,
                                   species="Other",
                                   factor=NA,
                                   molality=NA,
                                   percent=toString(percent),
                                   stringsAsFactors=FALSE))
      }
        
      df_aq_cont <- rbind(df_aq_cont, df)
    }
  }
  
  df_aq_cont <- df_aq_cont[, c("sample", "basis", "species", "factor", "molality", "percent")]
    
  return(df_aq_cont)

}


# function to create report versions of data categories (aq distributions, etc.)
create_report_df <- function(data, category, out_type){

  df_cat <- lapply(data, `[[`, category)

  all_species <- unique(unlist(lapply(lapply(df_cat, FUN=t), FUN=colnames)))

  df <- read.csv(text=paste(all_species, collapse="\t"), check.names=FALSE, sep="\t", stringsAsFactors=FALSE)

  for(i in 1:length(df_cat)){
    row <- as.data.frame(t(df_cat[[i]])[out_type, , drop=FALSE], stringsAsFactors=FALSE)
    df <- bind_rows(mutate_all(df, as.character), mutate_all(row, as.character))
  }
  df <- df[, order(colnames(df))]
  df_cat <- cbind.data.frame(sample=names(df_cat), df,  stringsAsFactors = FALSE)
  
  return(df_cat)
    
}

# function to compile a report
compile_report <- function(data, csv_filename, aq_dist_type, mineral_sat_type,
                           redox_type, get_aq_dist, get_mineral_sat, get_redox,
                           get_charge_balance, get_ion_activity_ratios, get_fugacity,
                           get_affinity_energy, input_processed_df,
                           df_input_processed_names){
    
  report_list <- list()
    
  # open processed input file and initialize report with it
  report <- input_processed_df
  names(report) <- df_input_processed_names
  #report <- read.csv(csv_filename, check.names=FALSE, stringsAsFactors=FALSE)
   
  report_list[["divs"]][["input"]] <- names(report)[2:length(report)] # start at 2 to exclude "sample" column
    
  # create report versions of EQ3 output blocks
  if(get_aq_dist){
    aq_distribution <- create_report_df(data=data, category='aq_distribution', out_type=aq_dist_type)
    report_list[["divs"]][["aq_distribution"]] <- names(aq_distribution)[2:length(aq_distribution)] # start at 2 to exclude "sample" column
    report <- report %>% inner_join(aq_distribution, by=c("Sample"="sample"))
  }
  
  if(get_mineral_sat){
    mineral_sat <- create_report_df(data=data, category='mineral_sat', out_type=mineral_sat_type)
    report_list[["divs"]][["mineral_sat"]] <- names(mineral_sat)[2:length(mineral_sat)] # start at 2 to exclude "sample" column
    report <- report %>% inner_join(mineral_sat, by=c("Sample"="sample"))
  }
  
  if(get_redox){
    redox <- create_report_df(data=data, category='redox', out_type=redox_type)
    report_list[["divs"]][["redox"]] <- names(redox)[2:length(redox)] # start at 2 to exclude "sample" column
    report <- report %>% inner_join(redox, by=c("Sample"="sample"))
  }
  
  if(get_charge_balance){
    charge_balance <- create_report_df(data=data, category='charge_balance', out_type=1)
    report_list[["divs"]][["charge_balance"]] <- names(charge_balance)[2:length(charge_balance)] # start at 2 to exclude "sample" column
    report <- report %>% inner_join(charge_balance, by=c("Sample"="sample"))
  }

  if(get_ion_activity_ratios){
    ion_activity_ratios <- create_report_df(data=data, category='ion_activity_ratios', out_type=1)
    report_list[["divs"]][["ion_activity_ratios"]] <- names(ion_activity_ratios)[2:length(ion_activity_ratios)] # start at 2 to exclude "sample" column
    report <- report %>% inner_join(ion_activity_ratios, by=c("Sample"="sample"))
  }

  if(get_fugacity){
    fugacities <- create_report_df(data=data, category='fugacity', out_type=1)
    report_list[["divs"]][["fugacity"]] <- names(fugacities)[2:length(fugacities)] # start at 2 to exclude "sample" column
    report <- report %>% inner_join(fugacities, by=c("Sample"="sample"))
  }

  if(get_affinity_energy){
    affinity <- create_report_df(data=data, category='affinity_energy', out_type=1)
    names(affinity)[2:length(names(affinity))] <- paste0(names(affinity)[2:length(names(affinity))], "_affinity")
    report_list[["divs"]][["affinity"]] <- names(affinity)[2:length(affinity)] # start at 2 to exclude "sample" column
    energy <- create_report_df(data=data, category='affinity_energy', out_type=2)
    names(energy)[2:length(names(energy))] <- paste0(names(energy)[2:length(names(energy))], "_energy")
    report_list[["divs"]][["energy"]] <- names(energy)[2:length(energy)] # start at 2 to exclude "sample" column
    report <- report %>%
      inner_join(affinity, by=c("Sample"="sample")) %>%
      inner_join(energy, by=c("Sample"="sample"))
  
  }

  rownames(report) <- report$Sample
  report$Sample <- NULL

  report_list[["report"]] <- report
    
  return(report_list)

}


### main
main_3o_mine <- function(files_3o,
                         rxn_filename,
                         get_aq_dist,
                         get_mass_contribution,
                         get_mineral_sat,
                         get_redox,
                         get_charge_balance,
                         get_ion_activity_ratios,
                         get_fugacity,
                         get_affinity_energy,
                         not_limiting,
                         mass_contribution_other,
                         csv_filename,
                         aq_dist_type,
                         mineral_sat_type,
                         redox_type,
                         input_filename,
                         batch_3o_filename,
                         df_input_processed,
                         df_input_processed_names,
                         custom_obigt,
                         verbose){
    
  start_time <- Sys.time()

  # allow user to add their custom data as an OBIGT
  if(!is.null(custom_obigt)){
    custom_obigt <- read.csv(custom_obigt, stringsAsFactors=F)
    custom_obigt <- custom_obigt[, c("name", "abbrv", "formula", "state", "ref1", "ref2", "date", "E_units", "G", "H", "S", "Cp", "V", "a1.a", "a2.b", "a3.c", "a4.d", "c1.e", "c2.f", "omega.lambda", "z.T")]
    suppressMessages(mod.OBIGT(custom_obigt, replace=T))
  }

  rxn_table <- NULL
  if(get_affinity_energy){
    # read table of reactions
    rxn_table <- readLines(rxn_filename)
  }

  # instantiate an empty object to store data from all 3o files
  batch_3o <- list()
  
  if(verbose > 1){
    writeLines("Now processing EQ3 output files...")
  }

  # process each .3o file
  for(file in files_3o){
    # add this sample's aqueous data to list of all sample data
    sample_3o <- mine_3o(file,
                         rxn_table=rxn_table,
                         get_aq_dist=get_aq_dist,
                         get_mass_contribution=get_mass_contribution,
                         get_mineral_sat=get_mineral_sat,
                         get_redox=get_redox,
                         get_charge_balance=get_charge_balance,
                         get_ion_activity_ratios=get_ion_activity_ratios,
                         get_fugacity=get_fugacity,
                         get_affinity_energy=get_affinity_energy,
                         not_limiting=not_limiting,
                         mass_contribution_other=mass_contribution_other,
                         verbose=verbose)

    # if this file could be processed, add its data to the batch_3o object
    if(length(sample_3o)>1){
      batch_3o[["sample_data"]][[sample_3o[["name"]]]] <- sample_3o
    }
  }
    
  if(verbose > 1){
    writeLines("Finished processing EQ3 output files...")
  }
  
  # compile aqueous contribution data into a single melted dataframe and
  # append it to the batch_3o object.
  if(get_mass_contribution && length(batch_3o)>0){
    if(verbose > 1){
      writeLines("Now processing mass contribution data...")
    }
    batch_3o[["mass_contribution"]] <- melt_mass_contribution(batch_3o=batch_3o,
                                                              other=mass_contribution_other,
                                                              verbose=verbose)
    if(verbose > 1){
      writeLines("Finished processing mass contribution data...")
    }
  }
    
  if(length(batch_3o)>0){
    # create a report summarizing 3o data from all samples
    report_list <- compile_report(data=batch_3o[["sample_data"]],
                                  csv_filename=csv_filename,
                                  aq_dist_type,
                                  mineral_sat_type,
                                  redox_type,
                                  get_aq_dist,
                                  get_mineral_sat,
                                  get_redox,
                                  get_charge_balance,
                                  get_ion_activity_ratios,
                                  get_fugacity,
                                  get_affinity_energy,
                                  df_input_processed,
                                  df_input_processed_names)
    

    # add the report to the batch_3o object
    report <- report_list[["report"]]
    batch_3o[["report"]] <- report
    batch_3o[["report_divs"]] <- report_list[["divs"]]
  }else{
    return(list())
  }

  # store user input file data
  batch_3o[["input"]] <- read.csv(input_filename, check.names=FALSE, stringsAsFactors=FALSE)

  # save the batch_3o object as an rds file
  if(!is.null(batch_3o_filename)){
    saveRDS(batch_3o, file=batch_3o_filename)
  }
      
  time_elapsed <- Sys.time() - start_time
  if(verbose > 1){
    writeLines(paste("Finished mining .3o files. Time elapsed:", round(time_elapsed, 2), "seconds"))
  }

  return(batch_3o)
}