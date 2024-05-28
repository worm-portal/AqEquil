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
                    this_pressure,
                    get_aq_dist=T,
                    get_mass_contribution=T,
                    get_mineral_sat=T,
                    get_redox=T,
                    get_charge_balance=T,
                    get_ion_activity_ratios=T,
                    get_fugacity=T,
                    get_basis_totals=T,
                    get_solid_solutions=T,
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
#     if(verbose > 0){
#       writeLines(paste0("\nSample ", this_name, " experienced errors during speciation:"))
#       output_error <- str_extract_all(extractme, regex("\\* Error.*?\n(\n|$)", dotall=T))
#       output_error <- lapply(output_error, function(x) gsub("\n", " ", x))
#       output_error <- lapply(output_error, str_squish)
        
#       for(i in 1:length(output_error)){
#         writeLines(output_error[[i]])
#       }
      
#     }

    setwd("../")
    return(list())
  }

  sample_3o <- list()

  sample_3o[["filename"]] <- this_file
  sample_3o[["name"]] <- this_name
  
  ### Begin mining temperature, pressure, water properties
    
  # mine params
  sample_3o[["temperature"]] <- isolate_block(str=extractme, begin_str="^.*Temperature=\\s+", end_str="\\s+.*$")
  sample_3o[["pressure"]] <- this_pressure #isolate_block(str=extractme, begin_str="^.*Pressure=\\s+", end_str="\\s+.*$")
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

    if(get_solid_solutions){
      if(grepl("--- Saturation States of Hypothetical Solid Solutions ---", extractme)){
        # string to isolate the solid solution saturation section:
        front_trim <- "^.*\n\n\n                --- Saturation States of Hypothetical Solid Solutions ---\n\n"

        # isolate solid solution block
        ss_block <- isolate_block(str=extractme, begin_str=front_trim, end_str="\n\n                     --- Fugacities ---.*$") # end_str might be different if fugacity section is not present

        if(ss_block != " None"){
          
        # split into substrings, each representing a separate solid solution
        ss_block <- strsplit(ss_block, "\n\n\n                --- ")[[1]]

        
          
        ss_entries <- list()
        for(ss_entry in ss_block){
          ss_entry <- strsplit(ss_entry, " ---\n\n   ")[[1]]
          ss_name <- ss_entry[1]
          ss_name <- sub("\n                --- ", "", ss_name) # needed to clean up the first ss_entry name
          ss_data <- ss_entry[2]
          ss_data <- sub("Ideal solution\n\n    Component                    x           Log x   Log lambda  Log activity\n\n", "", ss_data)[[1]]
          ss_entry <- strsplit(ss_data, "\n\n\n    Mineral                       Log Q/K         Aff, kcal    State\n\n")



          ss_entry <- lapply(ss_entry, FUN=strsplit, "\n")[[1]]
          ss_entry <- lapply(ss_entry, FUN=trimws, "l")
          ss_entry <- lapply(ss_entry, FUN=strsplit, "[ ]{2,}", perl=TRUE)
            
          names(ss_entry) <- c("ideal solution", "mineral")


          ideal_sol_df <- data.frame(component=character(),
                                   x=numeric(),
                                   `Log x`=numeric(),
                                   `Log lambda`=numeric(),
                                   `Log activity`=numeric(),
                                   stringsAsFactors = FALSE)
          for(row in ss_entry[["ideal solution"]]){

            x <- suppressWarnings(as.numeric(row[2]))
            if(is.na(x)){
              x <- 0
            }

            ideal_sol_df <- rbind(ideal_sol_df, data.frame(component=as.character(row[1]),
                                                         x = x,
                                                         `Log x` = as.numeric(row[3]),
                                                         `Log lambda` = as.numeric(row[4]),
                                                         `Log activity` = as.numeric(row[5])),
                                                         stringsAsFactors=FALSE)
          }
          names(ideal_sol_df) <- c("component", "x", "Log x", "Log lambda", "Log activity") # rename columns because check.names doesn't work when creating the dataframes
          ss_entry[["ideal solution"]] <- ideal_sol_df

          mineral_ss_df <- data.frame(mineral=character(),                                   
                                    `Log Q/K`=numeric(),
                                    `Aff, kcal`=numeric(),
                                    `State`=character(),
                                    stringsAsFactors = FALSE)
          for(row in ss_entry[["mineral"]]){


            if(length(row) == 3){
              row <- c(row, "")
            }
            mineral_ss_df <- rbind(mineral_ss_df, data.frame(mineral=as.character(row[1]),
                                                             `Log Q/K`=as.numeric(row[2]),
                                                             `Aff, kcal`=as.numeric(row[3]),
                                                             `State`=row[4]),
                                                             stringsAsFactors = FALSE)
          }
          names(mineral_ss_df) <- c("mineral", "Log Q/K", "Aff, kcal", "State") # rename columns because check.names doesn't work when creating the dataframes
          ss_entry[["mineral"]] <- mineral_ss_df

          ss_entries[[ss_name]] <- ss_entry
        }
        sample_3o[["solid_solutions"]] <- ss_entries
      }else{
        sample_3o[["solid_solutions"]] <- NA
        }
      }
    }
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
      
    if (!identical(ion_ratio_block_split, character(0))){
      
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
      
    # fix an annoying behavior with R dataframes ignoring rownames when there is
    # only one row
    if(nrow(df) == 1){
      rownames(df) <- unlist(lapply(split_str, `[[`, 1))[1]
      df <- df[ , !(names(df)=="gas"), drop=FALSE]
    }
      
    sample_3o[["fugacity"]] <- df
      
  }

                             
  ### begin sensible composition mining ("basis totals")
  if(get_basis_totals){
    sc_block <- isolate_block(extractme, "^.*--- Sensible Composition of the Aqueous Solution ---\n\n", "\n\n   The above data have.*$")
    str <- str_squish(strsplit(sc_block, "\n")[[1]])
    str <- str[3:length(str)]
    split_str <- strsplit(str, " ")
    sc_names <- unlist(lapply(lapply(split_str, `[[`, 1), paste0, "_total"))
    df <- data.frame(species=sc_names,
                     `mg/L`=as.numeric(unlist(lapply(split_str, `[[`, 2))),
                     `mg/kg.sol`=as.numeric(unlist(lapply(split_str, `[[`, 3))),
                     `molarity`=as.numeric(unlist(lapply(split_str, `[[`, 4))),
                     `molality`=as.numeric(unlist(lapply(split_str, `[[`, 5))),
                     stringsAsFactors=F, row.names=1)
    
    # fix an annoying behavior with R dataframes ignoring rownames when there is
    # only one row
    if(nrow(df) == 1){
      rownames(df) <- sc_names[1]
      df <- df[ , !(names(df)=="species"), drop=FALSE]
    }
      
    sample_3o[["basis_totals"]] <- df

  }
                             
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
                           get_basis_totals, input_processed_df,
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
    if('ion_activity_ratios' %in% names(data)){
    ion_activity_ratios <- create_report_df(data=data, category='ion_activity_ratios', out_type=1)
    report_list[["divs"]][["ion_activity_ratios"]] <- names(ion_activity_ratios)[2:length(ion_activity_ratios)] # start at 2 to exclude "sample" column
    report <- report %>% inner_join(ion_activity_ratios, by=c("Sample"="sample"))
    }
  }
    
  if(get_fugacity){
    fugacities <- create_report_df(data=data, category='fugacity', out_type=1)
    report_list[["divs"]][["fugacity"]] <- names(fugacities)[2:length(fugacities)] # start at 2 to exclude "sample" column
    report <- report %>% inner_join(fugacities, by=c("Sample"="sample"))
  }
    
  if(get_basis_totals){
    sc <- create_report_df(data=data, category='basis_totals', out_type=4) # 4 is the molality column
    report_list[["divs"]][["basis_totals"]] <- names(sc)[2:length(sc)] # start at 2 to exclude "sample" column
    report <- report %>% inner_join(sc, by=c("Sample"="sample"))
  }
    
  rownames(report) <- report$Sample
  report$Sample <- NULL

  report_list[["report"]] <- report
    
  return(report_list)

}


### main
main_3o_mine <- function(files_3o,
                         get_aq_dist,
                         get_mass_contribution,
                         get_mineral_sat,
                         get_redox,
                         get_charge_balance,
                         get_ion_activity_ratios,
                         get_fugacity,
                         get_basis_totals,
                         get_solid_solutions,
                         mass_contribution_other,
                         csv_filename,
                         aq_dist_type,
                         mineral_sat_type,
                         redox_type,
                         input_filename,
                         input_pressures,
                         batch_3o_filename,
                         df_input_processed,
                         df_input_processed_names,
                         verbose){
    
  start_time <- Sys.time()

  # instantiate an empty object to store data from all 3o files
  batch_3o <- list()
  
  if(verbose > 1){
    writeLines("Now processing EQ3 output files...")
  }
    
  names(input_pressures) <- files_3o
    
  # process each .3o file
  for(file in files_3o){
      
    # add this sample's aqueous data to list of all sample data
    sample_3o <- mine_3o(file,
                         this_pressure=input_pressures[file],
                         get_aq_dist=get_aq_dist,
                         get_mass_contribution=get_mass_contribution,
                         get_mineral_sat=get_mineral_sat,
                         get_redox=get_redox,
                         get_charge_balance=get_charge_balance,
                         get_ion_activity_ratios=get_ion_activity_ratios,
                         get_fugacity=get_fugacity,
                         get_basis_totals=get_basis_totals,
                         get_solid_solutions=get_solid_solutions,
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
                                  get_basis_totals,
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