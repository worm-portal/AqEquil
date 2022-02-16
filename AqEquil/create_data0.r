# load packages
suppressMessages({
  library(CHNOSZ)
  library(dplyr)
  library(comprehenr)
  library(stringr)
})

############ Helper functions

# print messages if 'verbose' setting >= vlevel of message.
vmessage <- function(m, vlevel, verbose){
  if(verbose >= vlevel){
    print(m)
  }
}


# function to create a string of a certain length by adding additional spaces.
# e.g., "H2O" becomes "H2O    " if nspaces=7.
# Spaces can be added before the string by specifying spaces_after=FALSE
fillspace <- function(str, nspaces, spaces_after=TRUE){
    
  ifelse(spaces_after,
         paste0(str, paste(rep(" ", nspaces-nchar(str)), collapse="")),
         paste0(paste(rep(" ", nspaces-nchar(str)), collapse=""), str))
}


# function to specify how many decimals are printed
# e.g. 12.433 becomes "12.4330" if k=4
s_d <- function(x, k) trimws(format(round(x, k), nsmall=k, scientific=F))


# main function
create_data0 <- function(thermo_df,
                         filename_ss=NULL,
                         db,
                         water_model,
                         template,
                         dissrxns,
                         basis_pref=c(),
                         exceed_Ttr=FALSE,
                         fixed_species=c("H2O", "H+", "O2(g)", "water", "Cl-", "e-"),
                         verbose){
    
  # set water model
  water(water_model)
    
  # initialize lists and vectors
  azero_vec <- c()
  neutral_ion_type_vec <- c()
  dissociation_list <- list()
  tag_vec <- c()
    
  for(i in 1:nrow(thermo_df)){
      # for each row in the data file...

      # look up azero bdot param
      azero_temp <- thermo_df[i, "azero"]
      names(azero_temp) <- thermo_df[i, "name"]
      azero_vec <- c(azero_vec, azero_temp)

      # look up neutral ion type
      nit_temp <- thermo_df[i, "neutral_ion_type"]
      names(nit_temp) <- thermo_df[i, "name"]
      neutral_ion_type_vec <- c(neutral_ion_type_vec, nit_temp)
      
      # look up dissociation reaction
      species_name <- thermo_df[i, "name"]
      
      if(thermo_df[i, "dissrxn"] != "" && is.null(dissrxns[[species_name]])){
        dissrxn <- thermo_df[i, "dissrxn"]
      }else if(!(is.null(dissrxns[[species_name]]))){
        # use auto-balanced dissociation reaction if no dissociation reaction is found for non-basis
        dissrxn <- dissrxns[[species_name]]
      }else if(thermo_df[i, "tag"] == "basis"){
        # pass
      }else{
        message("Error: dissociation reaction could not be generated for:")
        message(thermo_df[i, "name"])
      }
      
      if(thermo_df[i, "tag"] != "basis"){
        dissrxn <- strsplit(dissrxn, " ")[[1]] # split the rxn into coeffs and species
        #dissrxn_temp <- dissrxn[3:length(dissrxn)] # ignore the and diss. spec. name and coeff
        dissrxn_names <- dissrxn[c(FALSE, TRUE)] # get names of reactants and products
        dissrxn_coefs <- dissrxn[c(TRUE, FALSE)] # get coeffs of reactants and products
        dissrxn_coefs <- as.numeric(dissrxn_coefs) # convert coeffs from str to numeric
        names(dissrxn_coefs) <- dissrxn_names # set vector names to react and prod names
        dissociation_list[[species_name]] <- dissrxn_coefs # assign to dissociation list
      }
    
      # look up tag
      tag_temp <- thermo_df[i, "tag"]
      names(tag_temp) <- thermo_df[i, "name"]
      tag_vec <- c(tag_vec, tag_temp)
  }
    

  if(!is.null(filename_ss)){
    ss_params <- read.csv(filename_ss, stringsAsFactors=FALSE)
  }else{
    ss_params <- data.frame()
  }

  # initialize vector of name differences between OBIGT and SLOP
  CHNOSZ_data0_name_diff <- c()

  add_obigt_df <- thermo_df

  # remove raw line indicators '\r' in data0.min template
  data0_template <- gsub("\r", "", template)
      
  # initialize a vector to store names of species that must be skipped
  # due to one or more NA in its dissrxn logK grid
  skipped_species <- c()
    
  # loop through species in OBIGT file
  for(idx in 1:nrow(thermo_df)){

    entry <- thermo_df[idx, ]
    name <- entry$name

    if (name %in% c("O2(g)", "H2O", "H+")){
      vmessage(paste(name, "is included as a basis species by default. Moving to the next species..."), 2, verbose)
      next
    }
    if (entry$state %in% paste0("cr", 2:20)){
      next
    }
      
    date <- entry$date
    date <- strtrim(date, 9) # truncates date if greater than 9 letters
    elem <- makeup(thermo_df[idx, "formula_modded"]) # get elemental composition
      
    aux_basis <- FALSE
    if(name %in% lapply(dissrxns[["basis_list"]], `[[`, 1)){
      # if this species is marked as a preferred basis species, move to the next species
      vmessage(paste0("'", name, "' (basis species) processed successfully."), 2, verbose)
      next
    }else if(thermo_df[thermo_df[, "name"]==name, "tag"] == "basis"){
      # if this is marked as a basis in the data0 supplemental file, move to the next species
      vmessage(paste0("'", name, "' (basis species) processed successfully."), 2, verbose)
      next
    }else if(thermo_df[thermo_df[, "name"]==name, "tag"] == "aux"){
      # if this species is an auxiliary basis species, flag and continue with aqueous formatting
      aux_basis <- TRUE
    }

    # format charge for this species' data0 entry
    if("Z" %in% names(elem)){
      charge <- elem["Z"]
      elem <- elem[which(names(elem) != "Z")]
      formatted_charge <- format(round(charge, 1), nsmall = 1, scientific=F)
    } else {
      formatted_charge <- "0.0"
    }
      
    # format the element block of this species' data0 entry
    elem_list <- c()
        
    for(i in 1:length(names(elem))){

      # get element value and name from makeup
      elem_val <- s_d(elem[i], 4)
      elem_name <- names(elem)[i]

      # conditional formatting based on position
      if(i == 1 | i %% 4 == 0){ # first entry of a line
        max_length <- 8
        end_char <- ""
        if(nchar(elem_name) != 2){elem_name <- paste0(elem_name, " ")}
      }else if(i %% 3 == 0 && i != length(names(elem))){ # last entry of a line
        max_length <- 15
        end_char <- "\n"
      } else {
        max_length <- 15
        end_char <- ""
        if(nchar(elem_name) != 2){elem_name <- paste0(elem_name, " ")}
      }

      # paste together value and element name
      pasted_entry <- paste(elem_val, elem_name)

      # get decimal position and format spaces accordingly
      decimal_position <- gregexpr(pattern ='\\.', pasted_entry)[[1]]
      pasted_entry <- paste0(paste(rep(" ", max_length-decimal_position), collapse=""), pasted_entry, end_char)

      # add entry to element list
      elem_list <- c(elem_list, pasted_entry)
    }
      
    n_elements <- as.character(length(elem_list))
    element_list <- paste(elem_list, collapse="")
      
    # format the dissociation reaction block of this species' data0 entry
    if(name %in% names(dissociation_list)){
      species_name_list <- names(dissociation_list[[name]])
      species_val_list <- unname(dissociation_list[[name]])
    }

    n_species <- length(species_name_list)

    nchar <- 25 # spaces taken by entries in dissociation reaction (25 char)

    # convert species names to their data0 counterparts
    for(species in species_name_list){
      if(species %in% names(CHNOSZ_data0_name_diff)){
        species_name_list[which(species_name_list == species)] <- CHNOSZ_data0_name_diff[species]
      }
    }
    
    logK_list <- paste0("logK_grid_", name, "")
      
    # loop through species in dissociation reaction and format for data0
    spec_list <- c()
    for(i in 1:n_species){
      # get species value and name
      species_val <- format(round(species_val_list[i], 4), nsmall=4, scientific=F)
      species_name <- species_name_list[i]
      # conditional formatting based on position
      if(i == 1 | i %% 2 != 0){ # first entry of a line
        max_length <- 7
        end_char <- ""
        species_name <- fillspace(species_name, nchar)
      }else if(i %% 2 == 0 && i != n_species){ # last entry of a line
        max_length <- 8
        end_char <- "\n"
      } else {
        max_length <- 8
        end_char <- ""
        species_name <- fillspace(species_name, nchar)
      }
      # paste together coeff and element name
      pasted_entry <- paste0(species_val, "  ", species_name)
      # get decimal position and format spaces accordingly
      decimal_position <- gregexpr(pattern ='\\.', pasted_entry)[[1]]
      pasted_entry <- paste0(paste(rep(" ", max_length-decimal_position), collapse=""), pasted_entry, end_char)
      # add entry to element list
      spec_list <- c(spec_list, pasted_entry)
    }
    
    # instantiate template and begin formatting aq, cr, gas, liq entries
    species_list <- paste(spec_list, collapse="")
    n_species <- as.character(n_species) # convert to string
    template <- "+--------------------------------------------------------------------\n%s\n    date last revised = %s\n%s\n     charge  =   %s\n%s\n     %s element(s):\n%s\n     %s species in aqueous dissociation reaction:      \n%s\n**** logK grid [T, P @ Miscellaneous parameters]=%s\n%s"
    if(entry$state == "aq"){
      formatted_name = name
      volume = "*"
      if(aux_basis){
        keys = " keys   = aux              active"
        insertline_regex <- "\\+-+\naqueous species"
        insertline <- "+--------------------------------------------------------------------\naqueous species"
      } else {
        keys = " keys   = aqueous          active"
        insertline_regex <- "\\+-+\nsolids"
        insertline <- "+--------------------------------------------------------------------\nsolids"
      }
    } else if (entry$state == "cr"){
      formatted_name = paste0(fillspace(name, 24), entry$formula)
      tag = tag_vec[entry$name] # for solids, this is a tag like "refsate", "polymorph", or "idealized"
      formatted_tag = fillspace(tag, 17)
      keys = sprintf(" keys   = solid            %sactive", formatted_tag)
      formatted_V0PrTr = fillspace(entry$V, 9, spaces_after=FALSE)
      volume = sprintf("       V0PrTr = %s cm**3/mol", formatted_V0PrTr)
      insertline_regex <- "\\+-+\nliquids"
      insertline <- "+--------------------------------------------------------------------\nliquids"
    } else if (entry$state == "gas"){
      formatted_name = name
      tag = tag_vec[entry$name] # for gases, this is a tag like "refsate"
      formatted_tag = fillspace(tag, 17)
      keys = sprintf(" keys   = gas              %sactive", formatted_tag)
      volume = "       V0PrTr = 24465.000 cm**3/mol  (source = ideal gas law)"
      insertline_regex <- "\\+-+\nsolid solutions"
      insertline <- "+--------------------------------------------------------------------\nsolid solutions"
    } else if (entry$state == "liq"){
      formatted_name = paste0(fillspace(name, 24), entry$formula)
      tag = tag_vec[entry$name] # for liquids, this is a tag like "refsate"
      formatted_tag = fillspace(tag, 17)
      keys = sprintf(" keys   = liquid           %sactive", formatted_tag)
      formatted_V0PrTr = fillspace(entry$V, 9, spaces_after=FALSE)
      volume = sprintf("       V0PrTr = %s cm**3/mol", formatted_V0PrTr)
      insertline_regex <- "\\+-+\ngases"
      insertline <- "+--------------------------------------------------------------------\ngases"
    } else {
      # throw an error if entry is not aq, gas, cr, liq, or ss.
      stop(paste("Error: in", entry$name, "...", entry$state,
                 "is not a recognized state. Must be aq, gas, cr, liq, or ss."))
    }
      
    # append to aq, solid, gas, or liq portion of data0.min template
    output <- sprintf(template, formatted_name, date, keys, formatted_charge, volume, n_elements, element_list, n_species, species_list, name, logK_list)
    
    if(name == "O2"){
      O2_entry <- "\\+-+\nO2\n.*=O2\n         2.6560    3.0310    3.1080    3.0350\n         2.8740    2.6490    2.3540    1.8830"
      data0_template <- sub(O2_entry, paste0(output), data0_template)
    }else if(name == "OH-"){
      OH_entry <- "\\+-+\nOH-\n.*=OH-\n        14.9400   13.2710   12.2550   11.6310\n        11.2840   11.1670   11.3000   11.8280"
      data0_template <- sub(OH_entry, paste0(output), data0_template)
    }else if(name == "H2O(g)"){
      steam_entry <- "\\+-+\nH2O\\(g\\)\n.*=H2O\\(g\\)\n         2.2990    0.9950   -0.0060   -0.6630\n        -1.1560   -1.5340   -1.8290   -2.0630"
      data0_template <- sub(steam_entry, paste0(output), data0_template)
    }else{
      data0_template <- sub(insertline_regex, paste0(output, "\n", insertline), data0_template)
    }
      
    vmessage(paste0("'", name, "' processed successfully."), 2, verbose)
  }
    
  # handle basis species
  basis_entry_template <- "+--------------------------------------------------------------------\n%s\n    date last revised =  %s\n keys   = basis            active\n     charge  =   %s\n     %s element(s):\n%s"
    
  for(basis in lapply(dissrxns[["basis_list"]], `[[`, 1)){

    # go to the next basis species if it is among these hard-coded by EQ3:
    # (these are already in the data0.min template used to build all data0 files)
    if(basis %in% fixed_species){
      next
    }

    # get the date on the basis species entry
    suppressMessages({
      date_basis <- info(info(basis))$date
    })

    # get the basis species formula
    basis_formula <- add_obigt_df[add_obigt_df[, "name"] == basis, "formula_modded"][1]

    # get the elemental makeup of the basis species
    elem_list_basis <- makeup(basis_formula)

    # extract charge from elemental composition of basis species
    if("Z" %in% names(elem_list_basis)){
      charge_basis <- elem_list_basis["Z"]
      elem_list_basis <- elem_list_basis[which(names(elem_list_basis) != "Z")]
      formatted_charge_basis <- format(round(charge_basis, 1), nsmall = 1, scientific=F)
    } else {
      formatted_charge_basis <- "0.0"
    }
      
    # begin formatting for data0 entry
    n_elem_basis <- length(elem_list_basis)
    elem_list_basis_names <- names(elem_list_basis)
    elem_list_basis_coefs <- as.character(format(round(elem_list_basis, 4), nsmall = 4, scientific=F))
    elem_list_basis_formatted <- ""
              
    for(i in 1:length(elem_list_basis_names)){
        if(i == 1){
            sp <- "      "
        }else{
            sp <- "              "
        }
        elem_list_basis_formatted <- paste0(elem_list_basis_formatted, sp, elem_list_basis_coefs[i], " ", elem_list_basis_names[i])
    }

    basis_entry <- sprintf(basis_entry_template, basis, date_basis, formatted_charge_basis, n_elem_basis, elem_list_basis_formatted)

    # add basis entry to data0.min template
    basis_insertline_regex <- "\\+--------------------------------------------------------------------\nO2\\(g\\)"
    basis_insertline <- "+--------------------------------------------------------------------\nO2\\(g\\)"

    data0_template <- sub(basis_insertline_regex, paste0(basis_entry, "\n", basis_insertline), data0_template)
  }

  # handle elements

  # check for elements not already in data0.min template and append
  elem_template <- sub("^.*\nelements\n\\+-+", "", data0_template) # strip away string before element section in template
  elem_template <- sub("\n\\+-+\nbasis species.*$", "", elem_template) # strip away string after element section
  elem_temp_lines <- strsplit(elem_template, "\n")[[1]] # split element section by line
  elem_temp_names <- sub("\\s+.*$", "", elem_temp_lines) # isolate element names by stripping everything after them
  elem_addme <- setdiff(names(dissrxns[["basis_list"]]), elem_temp_names) # check which elements to add to template

  # REQUIRED BY EQ3: ensure that elements appear in the same order as the basis species representing those elements.
  elem_addme <- elem_addme[order(match(elem_addme,names(dissrxns[["basis_list"]])))]

  # loop through elements that need to be added to the data0 template
  for(elem in elem_addme){
      
    weight <- trimws(format(round(mass(elem), 5), nsmall=5, scientific=F))
      
    # format a line for the new element
    elem_temp_lines <- c(elem_temp_lines,
                         paste(elem,
                               paste(rep(" ", 16-(nchar(elem) + nchar(weight))), collapse=''),
                               weight,
                               collapse=""))
  }
    
  # join the elements back together into a full element block
  elem_block <- paste(elem_temp_lines, collapse="\n")

  # insert the full element block into the data0 template
  data0_template <- sub("elements\n\\+-+\n.*\n\\+-+\nbasis species",
                        paste0("elements\n+--------------------------------------------------------------------", elem_block, "\n+--------------------------------------------------------------------\nbasis species"),
                        data0_template)
    

  vmessage("Handling solid solutions...", 2, verbose)

  # handle solid solutions
  if(!is.null(filename_ss)){
    for(i in 1:nrow(ss_params)){

      entry <- ss_params[i, ]
      name <- entry$name

      vmessage(paste("Processing", name), 2, verbose)

      date <- entry$date
      species_name_list <- strsplit(entry$components, " ")[[1]]
      nextflag <- FALSE

      for(species in species_name_list){
        if(!(species %in% add_obigt_df$name)){
          vmessage(paste0("Error when the solid solution '", name, "': '", species, "' was not found in the data file as a pure mineral. Skipping it..."), 1, verbose)
          nextflag <- TRUE
          break
        } else if (species %in% skipped_species){
          vmessage(paste0("Error when processing the solid solution '", name, "': the dissociation reaction for '", species, "' contained one or more NAs in its logK grid. Skipping it..."), 1, verbose)
          nextflag <- TRUE
          break
        }
      }
      if(nextflag){next}

      n_species <- length(species_name_list)
      species_val_list <- rep(1, n_species)
      nchar <- 23 # spaces taken by entries in ss endmember list (23 char)

      # loop through species and format for data0
      spec_list <- c()
      for(i in 1:n_species){

        # get species value and name
        species_val <- format(round(species_val_list[i], 4), nsmall = 4, scientific=F)
        species_name <- species_name_list[i]

        # conditional formatting based on position
        if(i == 1 | i %% 2 != 0){ # first entry of a line
          max_length <- 7
          end_char <- ""
          species_name <- fillspace(species_name, nchar)
        }else if(i %% 2 == 0 && i != n_species){ # last entry of a line
          max_length <- 8
          end_char <- "\n"
        } else {
          max_length <- 8
          end_char <- ""
          species_name <- fillspace(species_name, nchar)
        }

        # paste together value and element name
        pasted_entry <- paste0(species_val, "  ", species_name)

        # get decimal position and format spaces accordingly
        decimal_position <- gregexpr(pattern ='\\.', pasted_entry)[[1]]
        pasted_entry <- paste0(paste(rep(" ", max_length-decimal_position), collapse=""), pasted_entry, end_char)

        # add entry to element list
        spec_list <- c(spec_list, pasted_entry)
      }
      species_list <- paste(spec_list, collapse="")

      n_species <- as.character(n_species) # convert to string

      ss_template <- "+--------------------------------------------------------------------\n%s\n    date last revised = %s\n%s\n  %s components\n%s\n type   = %s\n%s\n%s"

      model_param_list <- c(entry$mp1, entry$mp2, entry$mp3, entry$mp4, entry$mp5, entry$mp6)
      site_param_list  <- c(entry$sp1, entry$sp2, entry$sp3, entry$sp4,  entry$sp5, entry$sp6)

      model_param_list <- paste(s_d(model_param_list, 3), collapse=" ")
      site_param_list  <- paste(s_d(site_param_list, 3), collapse=" ")

      if(as.numeric(entry$n_model_params) > 0){
        model_param_template <- sprintf("  %s model parameters\n%s", entry$n_model_params, model_param_list)
      }else{
        model_param_template <- "  0 model parameters"
      }

      if(as.numeric(entry$n_site_params) > 0){
        site_param_template <- sprintf("  %s site parameters\n%s", entry$n_site_params, site_param_list)
      }else{
        site_param_template <- "  0 site parameters"
      }

      # ss:
      formatted_name = paste0(fillspace(name, 24), entry$formula)

      tag = entry$tag # for solid solutions, this is a tag like "cubic maclaurin", "ideal", "regular"
      formatted_tag = fillspace(tag, 17)
      keys = sprintf(" keys   = ss               %sactive", formatted_tag)

      insertline_regex <- "\\+-+\nreferences"
      insertline <- "+--------------------------------------------------------------------\nreferences"

      # create output for solid solution entries

      output <- sprintf(ss_template,
                        formatted_name,
                        date,
                        keys,
                        n_species,
                        species_list,
                        entry$type,
                        model_param_template,
                        site_param_template)

      data0_template <- sub(insertline_regex, paste0(output, "\n", insertline), data0_template)

    }
  } else {
    vmessage("No solid solutions supplied. Moving on...", 2, verbose)
  }

  # format basis and non-basis species for bdot parameter section
  bdot_formatted <- c()
  for(i in 1:length(azero_vec)){
      if(add_obigt_df[i, "name"] %in% c("Cl-", "O2", "OH-", "H2O", "H+")){
        next
      }
      
      if(add_obigt_df[i, "state"] == "aq"){
        spec_name <- names(azero_vec)[i]
        spec_azero <- as.character(format(round(azero_vec[i], 4), nsmall = 4, scientific=F))
        neutral_ion_type <- neutral_ion_type_vec[i]
        bdot_entry <- paste(spec_name,
              paste(rep(" ", 36-(nchar(spec_name) + nchar(spec_azero))), collapse=''),
              spec_azero,
              "  ",
              fillspace(neutral_ion_type, 2, spaces_after=FALSE),
              collapse="")

        # add bdot entry to data0.min template
        bdot_insertline_regex <- "\\+--------------------------------------------------------------------\nelements"
        bdot_insertline <- "+--------------------------------------------------------------------\nelements"
        data0_template <- sub(bdot_insertline_regex, paste0(bdot_entry, "\n", bdot_insertline), data0_template)

      }
  }

  return(data0_template)
    
}