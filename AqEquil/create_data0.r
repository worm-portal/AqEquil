# load packages
suppressMessages({
  library(CHNOSZ)
  library(dplyr)
  library(comprehenr)
})

known_oxstates <- c("H"="+", "O"="-2", "F"="-",
                    "Cl"="-", "Na"="+", "Li"="+",
                    "K"="+", "Ca"="+2", "Mg"="+2", "Z"="+")#,
#                   "Si"="+4", "Al"="+3") # useful for aluminosilicates


############ Helper functions

# helper function from Jeff Dick's CHNOSZ
# to retrieve the coefficients of reactions to form the species from the basis species
species.basis <- function(species=get("thermo", CHNOSZ)$species$ispecies) {
  # current basis elements
  bmat <- basis.elements()
  tbmat <- t(bmat)
  # what are the elements?
  belem <- rownames(tbmat)
  # get the species makeup into a matrix
  mkp <- as.matrix(sapply(makeup(species, count.zero=TRUE), c))
  # the positions of the species elements in the basis elements
  ielem <- match(rownames(mkp), belem)
  # the elements of the species must be contained by the basis species
  if(any(is.na(ielem))) stop(paste("element(s) not in the basis:", 
    paste(rownames(mkp)[is.na(ielem)], collapse=" ")))
  # the positions of the basis elements in the species elements
  jelem <- match(belem, rownames(mkp))
  # keep track of which ones are NA's; 
  # index them as one here but turn them to zero later
  ina <- is.na(jelem)
  jelem[ina] <- 1
  # now put the species matrix into the same order as the basis
  mkp <- mkp[jelem, , drop=FALSE]
  # fill zeros for any basis element not in the species
  mkp[ina, ] <- 0
  # solve for the basis coefficients and transpose
  nbasis <- t(apply(mkp, 2, function(x) solve(tbmat, x)))
  # very small numbers are probably a floating point artifact
  # can cause problems in situations where zeros are needed
  # (manifests as issue in longex("phosphate"), where which.balance()
  #  identifies H2O as conserved component)
  # 20140201 set digits (to R default) becuase getOption("digits") is changed in knitr
  out <- zapsmall(nbasis, digits=7)
  # add names of species and basis species
  colnames(out) <- colnames(tbmat)
  # add names of species only if it was a character argument
  if(all(is.character(species))) rownames(out) <- species
  return(out)
} 

                    
# modified version of subcrt() that only performs balancing
subcrt_bal <- function(ispecies, coeff){
    
  thermo <- get("thermo", CHNOSZ)
  species <- as.character(thermo$OBIGT$name[ispecies])
  state <- as.character(thermo$OBIGT$state[ispecies])
    
  # the mass balance; should be zero for a balanced reaction
  mss <- makeup(ispecies, coeff, sum=TRUE)
  # take out very small numbers
  mss[abs(mss) < 1e-7] <- 0
  # report and try to fix any non-zero mass balance
  if(any(mss!=0)) {
    # the missing composition: the negative of the mass balance
    miss <- -mss
    # drop elements that are zero
    miss <- miss[miss!=0]
    # look for basis species that have our compositoin
    tb <- thermo$basis
    if(all(names(miss) %in% colnames(tb)[1:nrow(tb)])) {
      # the missing composition as formula
      ft <- as.chemical.formula(miss)
      # the basis species needed to supply it
      bc <- species.basis(ft)
      # drop zeroes
      bc.new <- bc[,(bc[1,]!=0),drop=FALSE]
      bc <- bc.new
      ispecies.new <- tb$ispecies[match(colnames(bc),rownames(tb))]
      b.species <- thermo$OBIGT$formula[ispecies.new]
      newspecies <- c(ispecies, tb$ispecies[match(colnames(bc), rownames(tb))])
      newcoeff <- c(coeff, as.numeric(bc[1, ]))
      return(list('newspecies'=newspecies, 'newcoeff'=newcoeff))
    }
  } else {
    return(list('newspecies'=ispecies, 'newcoeff'=as.numeric(coeff)))
  }
}


# print messages if 'verbose' setting >= vlevel of message.
vmessage <- function(m, vlevel, verbose){
  if(verbose >= vlevel){
    message(m)
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
sd <- function(x, k) trimws(format(round(x, k), nsmall=k))


# trims away leading and trailing spaces and condenses multiple spaces between words
trimspace <- function(str){
  gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", str, perl=TRUE)
}


# function to return 0 if a value is NA
check_get <- function(val){
  ifelse(is.na(val), 0, val)
}


get_oxstate <- function(makeup_list){
    
  ifelse(length(setdiff(names(makeup_list), names(known_oxstates))) == 1,
  {
  ox_state <- c(- check_get(makeup_list["F"])*(-1) -
              check_get(makeup_list["Cl"])*(-1) -
              check_get(makeup_list["H"]) -
              check_get(makeup_list["Li"]) -
              check_get(makeup_list["K"]) -
              check_get(makeup_list["Na"]) -
              check_get(makeup_list["Ca"])*(2) -
              check_get(makeup_list["Mg"])*(2) -
#               check_get(makeup_list["Si"])*(4) - # useful for aluminosilicates
#               check_get(makeup_list["Al"])*(3) - # useful for aluminosilicates
              check_get(makeup_list["O"])*(-2) +
              check_get(makeup_list["Z"]))
  ox_state <- ox_state/check_get(makeup_list[setdiff(names(makeup_list), names(known_oxstates))])
  ox_state <- ifelse(ox_state %% 1 == 0, ox_state, NA) # filter out results with non-integer oxidation states
  names(ox_state) <- setdiff(names(makeup_list), names(known_oxstates))
  return(ox_state)
  },
  {
   return(NA)
  })
}


# function to print charge (e.g. a charge of +1 is displayed as "+", a charge of -2 is displayed as "-2")
format_charge <- function(charge){
  sign <- ifelse(charge >= 0, "+", "-")
  mag <- ifelse(abs(charge) == 1, "", as.character(abs(charge)))

  sign <- ifelse(charge == 0, "", sign)
  mag <- ifelse(charge == 0, "", mag)

  return(paste0(sign, mag))
}


errcheck <- function(expr, me="", mw="", verbose){
  tryCatch(expr,
    error = function(e){
      if(me != ""){
        vmessage(paste0("An error occurred:\n", me, "Additional detail:\n", e), 1, verbose)
        stop()
      }else{
        vmessage(paste0("An error occurred:\n", e), 1, verbose)
        stop()
      }
            
      },
      warning = function(w){
        if(mw != ""){
          vmessage(paste0("A warning occurred:\n", mw, "Additional detail:\n", w), 1, verbose)
        }else{
          vmessage(paste0("A warning occurred:\n", w), 1, verbose)
        }
      }
  )
}


match_basis_comp <- function(sp_elems, elem){
  if(elem %in% sp_elems){
    match <- all(sp_elems %in% unique(c(elem, "H", "O", "Z")))
  }else{
    match <- FALSE
  }
  return(match)
}


calc_bdot <- function(T){
  # Calculate bdot parameter at a given temperature.

  # GB notes:
  # The equation used by dbcreate to approximate the curve in Fig 3
  # of Helgeson 1969 results in numbers that are close to, but not
  # exactly the same as those in data0.jus:

  # Bdot parameter grid:
  #  0.0376   0.0443   0.0505   0.0529   0.0479   0.0322   0.0000   0.0000  # from dbcreate
  #  0.0374   0.0430   0.0460   0.0470   0.0470   0.0340   0.0000   0.0000  # from data0.jus

  # Close but not exact! data0.jus is closer to what is depicted in Fig 3 of Helgeson 1969.
  # Not sure what other equation to use, though. Will keep the dbcreate equation for now.
  # TODO: look into alternative equations.

  b1 <-  0.0374
  b2 <-  1.3569e-4
  b3 <-  2.6411e-7
  b4 <- -4.6103e-9  

  result <- b1 + b2*T + b3*(T-25.0)^2 + b4*(T-25.0)^3

  return(ifelse(T >= 300, 0, result))

}


################################ get_dissrxn

### INPUT EXAMPLE
# sp_name <- "tyrosine"
# # Preferred basis species for elements
# # (Not all elements have to be specified)
# basis_pref <- c("C"="HCO3-", "N"="NH4+")
# HOZ_balancers <- c("H+", "O2", "H2O") # might be data0-specific (e.g., "O2(g)")

get_dissrxn <- function(sp_name, basis_pref=c(), HOZ_balancers=c("H+", "O2", "H2O"), db_idx=c()){
    
  # get a vector of elements that make up the species
  basis_elem <- (function (x) setdiff(names(unlist(makeup(info(info(x), check.it=F)$formula))), c("H", "O", "Z"))) (sp_name)
  
  if (length(db_idx) == 0){
    sp <- thermo()$OBIGT
  } else {
    sp <- thermo()$OBIGT[db_idx, ]
  }
  sp_formula <- sp$formula
  sp_formula_makeup <- makeup(sp_formula)
  names(sp_formula_makeup) <- sp$name

  # loop through each element and assign preferred basis species (if supplied)
  # or automatically choose compatible basis species
  basis_list <- list()
  for(elem in basis_elem){

    # if a preferred basis species exists for this element, assign it and move to next element
    if(elem %in% names(basis_pref)){
      basis_list[[elem]] <- c(info(info(basis_pref[elem]), check.it=F)$name)
      names(basis_list[[elem]]) <- basis_pref[elem]
      next
    }

    # get index of species that fulfill requirements of EQ3's definition of a basis species
    idx <- unlist(lapply(lapply(sp_formula_makeup, names), match_basis_comp, elem=elem))

    # provide a vector of potential basis species for this element
    basis_list[[elem]] <- sp$name[idx]
    names(basis_list[[elem]]) <- sp$name[idx]

    # get formula makeup of potential basis species
    this <- lapply(unlist(lapply(lapply(lapply(basis_list[[elem]], FUN=info), FUN=info), `[`, "formula")), makeup)
        
    # if possible, narrow down this list to basis species with a matching element abundance of 1
    abund1 <- unlist(lapply(this, function(x) x[elem]==1))
    if(sum(abund1 != 0)){
      # assign narrowed list of matches to basis list
      basis_list[[elem]] <-  unlist(basis_list[[elem]][abund1])
    }

    # filter out species with "[" in the name (e.g. "[-CH3]"). These kinds of groups are in CHNOSZ's OBIGT.
    basis_list[[elem]] <-  unlist(basis_list[[elem]][!grepl("\\[", names(basis_list[[elem]]))])

    # error handling in case a suitable basis species cannot be found
    errcheck(length(basis_list[[elem]]) == 0, me=paste0("A suitable basis species could be found to represent the element ", elem, "."))

  }

  # further narrow down the list of potential basis species by grabbing the basis species
  # with the fewest elements+charge from available basis species. In the end there should
  # only be one basis species assigned to each element.

  smallest_elem_vec <- c()
  for(elem in basis_elem){
    form <- lapply(info(info(basis_list[[elem]]))$formula, makeup)
    sums <- unlist(lapply(lapply(form, abs), sum))
    idx <- which(sums==min(sums))
    if(length(idx) > 1){
      idx <- idx[1]
    }
    smallest_elem_vec[elem] <- basis_list[[elem]][idx]
  }

  dissrxns <- lapply(sp_name, spec_diss, smallest_elem_vec, sp_formula_makeup, HOZ_balancers)
  names(dissrxns) <- sp_name
  dissrxns[["basis_list"]] <- smallest_elem_vec
  return(dissrxns)

}

# supply a species name, get a dissrxn.
# (meant to lapply this over all species names)
spec_diss <- function(sp, smallest_elem_vec, sp_formula_makeup, HOZ_balancers){
  # determine a balanced dissociation reaction into basis species

  # get the elemental composition
  comp <- makeup(info(info(sp), check.it=F)$formula)

  # get the elements that are not H, O, Z (charge)
  basis_elem <- setdiff(names(comp), c("H", "O", "Z"))

  # set basis species in CHNOSZ so that a dissociation reaction can be balanced by subcrt()
  basis(unique(c(smallest_elem_vec[basis_elem], HOZ_balancers)))
      
  # loop through each element in the composition of this species. Divide by the coefficient of that element in its basis species.
  # e.g. If the species is S3-2 and the basis species for S is S2-, ensure the coefficient for S2- in the dissrxn is 3/2.
  for(elem in names(comp[basis_elem])){
    if(!(elem %in% c("Cl", "H", "O"))){
      n <- sp_formula_makeup[[smallest_elem_vec[[elem]]]][elem]
    }else{
      n <- 1
    }
    comp[elem] <- comp[elem]/n
  }

###### OPTION A. this portion works, but is slow. Uncomment A and comment B.
      
#     rxn_block <- subcrt(c(this_species, smallest_elem_vec[basis_elem]), c(-1, comp[basis_elem]))$reaction
#     spec_names <-  unlist(rxn_block$name)
#     spec_names[spec_names == "water"] <- "H2O" # prevent CHNOSZ from renaming 'H2O' (name req. by EQ3) to 'water'
#     sp_dissrxn <-  unlist(rxn_block$coeff)
#     names(sp_dissrxn) <- spec_names
      
###### END OF OPTION A
      
###### OPTION B. This chunk is fast but needs testing.
      
  if(identical(basis_elem, character(0))){
    calculate <- TRUE
  }else{
    calculate <- sp != smallest_elem_vec[basis_elem]
  }
  if (calculate){
    if(identical(basis_elem, character(0))){
      isp <- suppressMessages(info(sp))
      rxn_coeffs <- c(-1)
    }else{
      isp <- info(c(sp, smallest_elem_vec[basis_elem]))
      rxn_coeffs <- c(-1, comp[basis_elem])
    }
      
    autobalance <- subcrt_bal(isp, rxn_coeffs)
    spec_names <- thermo()$OBIGT[autobalance$newspecies, "name"]
    spec_names[spec_names == "water"] <- "H2O" # prevent CHNOSZ from renaming 'H2O' (name req. by EQ3) to 'water'
    sp_dissrxn <- autobalance$newcoeff
    names(sp_dissrxn) <- spec_names
  }else{
    sp_dissrxn <- c(-1, 1)
    names(sp_dissrxn) <- c(sp, smallest_elem_vec[basis_elem])
  }
      
###### END OF OPTION B
    
  sp_dissrxn <- paste(c(rbind(sprintf("%.4f", round(sp_dissrxn, 4)), names(sp_dissrxn))), collapse=" ")

  return(sp_dissrxn)
}

######################################## create_data0 function

# main function
create_data0 <- function(supp_file,
                         supp_file_ss=NULL,
                         grid_temps,
                         grid_press,
                         db,
                         template,
                         dissrxns,
                         modified_custom_obigt,
                         db_idx,
                         basis_pref=list(),
                         exceed_Ttr=FALSE,
                         verbose){

  # read the supplementary data file
  data0_params <- read.csv(supp_file, stringsAsFactors=FALSE)

  # initialize lists and vectors
  azero_vec <- c()
  neutral_ion_type_vec <- c()
  dissociation_list <- list()
  tag_vec <- c()
      
  for(i in 1:nrow(data0_params)){
      # for each row in the supplementary file...

      # look up azero bdot param
      azero_temp <- data0_params[i, "azero"]
      names(azero_temp) <- data0_params[i, "name"]
      azero_vec <- c(azero_vec, azero_temp)

      # look up neutral ion type
      nit_temp <- data0_params[i, "neutral_ion_type"]
      names(nit_temp) <- data0_params[i, "name"]
      neutral_ion_type_vec <- c(neutral_ion_type_vec, nit_temp)

      # get dissociation reaction
      if(data0_params[i, "tag"] != "basis"){
          dissrxn <- dissrxns[[data0_params[i, "name"]]]
          dissrxn <- strsplit(dissrxn, " ")[[1]] # split the rxn into coeffs and species
          species_name <- data0_params[i, "name"] # get the name of the dissociating species
          dissrxn_temp <- dissrxn[3:length(dissrxn)] # ignore the "-1" and diss. spec. name
          dissrxn_names <- dissrxn_temp[c(FALSE, TRUE)] # get names of reactants and products
          dissrxn_coefs <- dissrxn_temp[c(TRUE, FALSE)] # get coeffs of reactants and products
          dissrxn_coefs <- as.numeric(dissrxn_coefs) # convert coeffs from str to numeric
          names(dissrxn_coefs) <- dissrxn_names # set vector names to react and prod names
          dissociation_list[[species_name]] <- dissrxn_coefs # assign to dissociation list
      }

      # look up tag
      tag_temp <- data0_params[i, "tag"]
      names(tag_temp) <- data0_params[i, "name"]
      tag_vec <- c(tag_vec, tag_temp)
  }

  

  if(!is.null(supp_file_ss)){
    ss_params <- read.csv(supp_file_ss, stringsAsFactors=FALSE)
  }else{
    ss_params <- data.frame()
  }

  # initialize vector of name differences between OBIGT and SLOP
  CHNOSZ_data0_name_diff <- c()

  add_obigt_df <- modified_custom_obigt

  # check if temperature grid has the necessary 8 entries
  if(length(grid_temps) != 8){
      vmessage("Error: temperature grid does not have 8 entries.", 1, verbose)
      grid_temps <- c(0.0100, 50.0000, 100.0000, 150.0000, 200.0000, 250.0000, 300.0000, 350.0000)
      vmessage("Resorting to using a temperature grid of:", 1, verbose)
      vmessage(paste0(grid_temps), 1, verbose)
      grid_press <- "Psat" # "Psat" for liquid-vapor saturation curve from temperature grid
  }

  # round grid temperatures to four decimal places
  grid_temps <- round(grid_temps, 4)

  # calculate PSAT pressure if specified by user or if pressure grid
  # has a number of values that does not equal temperature grid length.
  if(tolower(grid_press) == "psat" | length(grid_press) != length(grid_temps)){
      vmessage("Calculating pressure grid along liquid-vapor saturation curve...", 2, verbose)
      grid_press <- water("Psat", T=grid_temps + 273.15)[[1]]
  }

  # remove raw line indicators '\r' in data0.min template
  data0_template <- gsub("\r", "", template)
      
  # initialize a vector to store names of species that must be skipped
  # due to one or more NA in its dissrxn logK grid
  skipped_species <- c()

  # loop through species in OBIGT file
  for(idx in db_idx){

    entry <- suppressMessages(info(idx))
    name <- entry$name
    date <- entry$date
    date <- strtrim(date, 9) # truncates date if greater than 9 letters
    elem <- makeup(idx) # get elemental composition

    if (name == "O2(g)"){
      vmessage("O2(g) is included as a basis species by default. Moving to the next species...", 2, verbose)
      next
    }

    aux_basis <- FALSE
    if(name %in% lapply(dissrxns[["basis_list"]], `[[`, 1)){
      # if this species is marked as a preferred basis species, move to the next species
      vmessage(paste0("'", name, "' (basis species) processed successfully."), 2, verbose)
      next
    }else if(data0_params[which(data0_params["name"]==name), "tag"] == "basis"){
      # if this is marked as a basis in the data0 supplemental file, move to the next species
      vmessage(paste0("'", name, "' (basis species) processed successfully."), 2, verbose)
      next
    }else if(data0_params[which(data0_params["name"]==name), "tag"] == "aux"){
      # if this species is an auxiliary basis species, flag and continue with aqueous formatting
      aux_basis <- TRUE
    }

    # format charge for this species' data0 entry
    if("Z" %in% names(elem)){
      charge <- elem["Z"]
      elem <- elem[which(names(elem) != "Z")]
      formatted_charge <- format(round(charge, 1), nsmall = 1)
    } else {
      formatted_charge <- "0.0"
    }
      
    # format the element block of this species' data0 entry
    elem_list <- c()
    for(i in 1:length(names(elem))){

      # get element value and name from makeup
      elem_val <- sd(elem[i], 4)
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
    species_name_list <- c(name) # reactant and product names
    species_val_list <- c(1) # stoiciometric coeffs
    for(i in 1:length(species_name_list)){
      # replace non-basis species in dissociation rxn with basis species
      species_name <- species_name_list[i]
      species_val <- species_val_list[i]
      if(species_name %in% names(dissociation_list)){
        species_name_list <- species_name_list[-i]
        species_val_list <- species_val_list[-i]
        species_name_list <- c(species_name_list, names(dissociation_list[[species_name]]))
        species_val_list <- c(species_val_list, unname(dissociation_list[[species_name]])*species_val)
      }
    }
    species_name_list <- c(name, species_name_list)
    species_val_list <- c(-1, species_val_list)
    n_species <- length(species_name_list)
 
    # format the logK reaction block of this species' data0 entry
    # This is done within a tryCatch() in case this fails.
    tryCatch({
      # the subcrt() calculation for each P-T in the grid
      if(entry$state == "cr"){
        # temporary fix for CHNOSZ not recognizing uppercase mineral names
        # when considering its thermo properties above its transition state.
        # Long term fix: have lowercase names be the WORM default.
        to_check <- suppressMessages(info(tolower(species_name_list[1])))
        if(!is.na(to_check)){
          logK_grid <- suppressMessages(subcrt(c(tolower(species_name_list[1]), species_name_list[2:length(species_name_list)]), species_val_list, T=grid_temps, P=grid_press, exceed.Ttr=exceed_Ttr)$out$logK)
        }
      } else {
          logK_grid <- suppressMessages(subcrt(species_name_list, species_val_list, T=grid_temps, P=grid_press, exceed.Ttr=exceed_Ttr)$out$logK)
      }
      # if CHNOSZ can't perform a calculation, assign a logK grid of zeros
      }, error=function(e){
        vmessage(paste0("Warning: CHNOSZ is unable to calculate a logK grid ",
                        "for the formation reaction of ", name,
                        ". A logK grid of zeros will be output."), 1, verbose)
        logK_grid <<- rep(0, length(grid_temps)) # assign global variable with <<- because this is within the error function
    })
    nchar <- 25 # spaces taken by entries in dissociation reaction (25 char)
    if(TRUE %in% is.na(logK_grid)){
      skipped_species <- c(skipped_species, species_name_list[1])
      vmessage(paste0("WARNING: One or more missing values are present in the logK grid calculated for ", species_name_list[1], ". This species will be skipped."), 1, verbose)
      next
    }
    # now that CHNOSZ has calculated logKs, convert species names to their data0 counterparts
    for(species in species_name_list){
      if(species %in% names(CHNOSZ_data0_name_diff)){
        species_name_list[which(species_name_list == species)] <- CHNOSZ_data0_name_diff[species]
      }
    }
    # loop through logK values and format for data0
    logK_list <- c()
    for(i in 1:length(logK_grid)){
      logK_val <- sd(logK_grid[i], 4)
      # conditional formatting based on position
      if(i == 1 | i %% 5 == 0){ # first entry of a line
        max_length <- 11
        end_char <- ""
      }else if(i %% 4 == 0 && i != length(logK_grid)){ # last entry of a line
        max_length <- 6
        end_char <- "\n"
      } else {
        max_length <- 6
        end_char <- ""
      }
      # get decimal position and format spaces accordingly
      decimal_position <- gregexpr(pattern ='\\.', logK_val)[[1]]
      logK_val <- paste0(paste(rep(" ", max_length-decimal_position), collapse=""), logK_val, end_char)
      # append to logk list
      logK_list <- c(logK_list, logK_val)
    }
    logK_list <- paste(logK_list, collapse="")
      
    # loop through species in dissociation reaction and format for data0
    spec_list <- c()
    for(i in 1:n_species){
      # get species value and name
      species_val <- format(round(species_val_list[i], 4), nsmall=4)
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
    data0_template <- sub(insertline_regex, paste0(output, "\n", insertline), data0_template)
    vmessage(paste0("'", name, "' processed successfully."), 2, verbose)
  }
    
  # handle basis species
  basis_entry_template <- "+--------------------------------------------------------------------\n%s\n    date last revised =  %s\n keys   = basis            active\n     charge  =   %s\n     %s element(s):\n%s"

  for(basis in lapply(dissrxns[["basis_list"]], `[[`, 1)){

    # go to the next basis species if it is among these hard-coded by EQ3:
    # (these are already in the data0.min template used to build all data0 files)
    if(basis == "O2(g)" | basis == "H2O" | basis == "Cl-" | basis == "H+"){
      next
    }

    # get the date on the basis species entry
    suppressMessages({
      date_basis <- info(info(basis))$date
    })

    # get the basis species formula
    basis_formula <- add_obigt_df[which(add_obigt_df[, "name"] == basis), "formula"]

    # get the elemental makeup of the basis species
    elem_list_basis <- makeup(basis_formula)

    # extract charge from elemental composition of basis species
    if("Z" %in% names(elem_list_basis)){
      charge_basis <- elem_list_basis["Z"]
      elem_list_basis <- elem_list_basis[which(names(elem_list_basis) != "Z")]
      formatted_charge_basis <- format(round(charge_basis, 1), nsmall = 1)
    } else {
      formatted_charge_basis <- "0.0"
    }

    # begin formatting for data0 entry
    n_elem_basis <- length(elem_list_basis)
    elem_list_basis_names <- names(elem_list_basis)
    elem_list_basis_coefs <- as.character(format(round(elem_list_basis, 4), nsmall = 4))
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
    weight <- trimws(format(round(mass(elem), 5), nsmall=5))

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
  if(!is.null(supp_file_ss)){
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
        species_val <- format(round(species_val_list[i], 4), nsmall = 4)
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

      model_param_list <- paste(sd(model_param_list, 3), collapse=" ")
      site_param_list  <- paste(sd(site_param_list, 3), collapse=" ")

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

      if(add_obigt_df[i, "state"] == "aq"){
        spec_name <- names(azero_vec)[i]
        spec_azero <- as.character(format(round(azero_vec[i], 4), nsmall = 4))
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


  # calculate debye huckel a and b parameters for the grid
  A_DH_grid <- unlist(water("A_DH", T=273.15+grid_temps, P=grid_press))
  B_DH_grid <- unlist(water("B_DH", T=273.15+grid_temps, P=grid_press)*10^-8)

  # format grid values
  grid_temps_f <- as.character(format(round(grid_temps, 4), nsmall = 4))
  grid_press_f <- as.character(format(round(grid_press, 4), nsmall = 4))
  A_DH_grid_f <- as.character(format(round(A_DH_grid, 4), nsmall = 4))
  B_DH_grid_f <- as.character(format(round(B_DH_grid, 4), nsmall = 4))

  bdot_grid_f <- as.character(format(round(calc_bdot(grid_temps), 4), nsmall = 4))


  # cco2 (coefficients for the drummond (1981) polynomial)
  # GB note: might not change with T or P?
  # Examination of various data0 files seems to indicate that DBcreate does not change these values.

  # Calculate the "log k for eh reaction" grid.
  # From eq. 9 in EQPT (version 7.0) user manual, part 2, by Wolery:
  suppressMessages({
      logK_Eh_vals <- subcrt(c("H2O", "O2", "e-", "H+"),
                            c(-2, 1, 4, 4),
                            c("liq", "gas", "aq", "aq"),
                            T=grid_temps,
                            P=grid_press)$out$logK
  })

  logk_grid_f <- as.character(format(round(logK_Eh_vals, 4), nsmall = 4))


  tempgrid <- c("     ")
  presgrid <- c("     ")
  A_DHgrid <- c("     ")
  B_DHgrid <- c("     ")
  bdotgrid <- c("     ")
  logkgrid <- c("     ")
  for(i in 1:8){
      if(i == 5){
          tempgrid <- c(tempgrid, "\n     ")
          presgrid <- c(presgrid, "\n     ")
          A_DHgrid <- c(A_DHgrid, "\n     ")
          B_DHgrid <- c(B_DHgrid, "\n     ")
          bdotgrid <- c(bdotgrid, "\n     ")
          logkgrid <- c(logkgrid, "\n     ")
      }
      tempgrid <- c(tempgrid, paste0(paste(rep(" ", 10-nchar(grid_temps_f[i])), collapse=""), grid_temps_f[i]))
      presgrid <- c(presgrid, paste0(paste(rep(" ", 10-nchar(grid_press_f[i])), collapse=""), grid_press_f[i]))
      A_DHgrid <- c(A_DHgrid, paste0(paste(rep(" ", 10-nchar(A_DH_grid_f[i])), collapse=""), A_DH_grid_f[i]))
      B_DHgrid <- c(B_DHgrid, paste0(paste(rep(" ", 10-nchar(B_DH_grid_f[i])), collapse=""), B_DH_grid_f[i]))
      bdotgrid <- c(bdotgrid, paste0(paste(rep(" ", 10-nchar(bdot_grid_f[i])), collapse=""), bdot_grid_f[i]))
      logkgrid <- c(logkgrid, paste0(paste(rep(" ", 10-nchar(logk_grid_f[i])), collapse=""), logk_grid_f[i]))

  }
  tempgrid <- paste(tempgrid, collapse="")
  presgrid <- paste(presgrid, collapse="")
  A_DHgrid <- paste(A_DHgrid, collapse="")
  B_DHgrid <- paste(B_DHgrid, collapse="")
  bdotgrid <- paste(bdotgrid, collapse="")
  logkgrid <- paste(logkgrid, collapse="")

  # insert temperature grid values into data0 template
  tempgrid_insertlines <- "\ntemperatures\n.*\npressures\n"
  data0_template <- sub(tempgrid_insertlines, paste0("\ntemperatures\n", tempgrid, "\npressures\n"), data0_template)

  # insert temperature grid values into data0 template
  presgrid_insertlines <- "\npressures\n.*\ndebye huckel a \\(adh\\)\n"
  data0_template <- sub(presgrid_insertlines, paste0("\npressures\n", presgrid, "\ndebye huckel a (adh)\n"), data0_template)

  # insert Debeye Huckel A and B parameter values into data0 template
  A_DHgrid_insertlines <- "\ndebye huckel a \\(adh\\)\n.*\ndebye huckel b \\(bdh\\)\n"
  data0_template <- sub(A_DHgrid_insertlines, paste0("\ndebye huckel a (adh)\n", A_DHgrid, "\ndebye huckel b (bdh)\n"), data0_template)
  B_DHgrid_insertlines <- "\ndebye huckel b \\(bdh\\)\n.*\nbdot\n"
  data0_template <- sub(B_DHgrid_insertlines, paste0("\ndebye huckel b (bdh)\n", B_DHgrid, "\nbdot\n"), data0_template)

  # insert bdot grid values into data0 template
  bdotgrid_insertlines <- "\nbdot\n.*\ncco2"
  data0_template <- sub(bdotgrid_insertlines, paste0("\nbdot\n", bdotgrid, "\ncco2"), data0_template)

  # insert logk (eh) grid values into data0 template
  logkgrid_insertlines <- "\nlog k for eh reaction\n.*\n\\+-+\nbdot parameters"
  logkgrid_end_insert <- "\n+--------------------------------------------------------------------\nbdot parameters"
  data0_template <- sub(logkgrid_insertlines, paste0("\nlog k for eh reaction\n", logkgrid, logkgrid_end_insert), data0_template)

  # modify the data0 header lines
  desc <- "data0.%s"
  min_desc <- "data0.min\nminimal working data0 file"
  data0_template <- sub(min_desc, sprintf(desc, db), data0_template)
  vmessage("Finished.", 2, verbose)
  write(data0_template, paste0("data0.", db))
}


######################################## Main function
                            
main_create_data0 <- function(filename,
                              supp_file,
                              supp_file_ss,
                              grid_temps,
                              grid_press,
                              db,
                              template,
                              exceed_Ttr,
                              data0_formula_ox_name,
                              suppress_redox,
                              infer_formula_ox,
                              generate_template,
                              template_name,
                              verbose,
                              basis_prefs,
                              basis_pref_names){

  # load thermodynamic data
  data0_params <- read.csv(supp_file, stringsAsFactors=F)
  rownames(data0_params) <- data0_params$name


  # get all unique oxidation states of elements in the entire dataset
  elem_ox <- c()
  for(entry in data0_params$formula_ox){
    entry <- strsplit(entry, " ")[[1]]
    for(elem in entry){
      e <- unlist(strsplit(elem, "[0-9](?=[A-Z])", perl=T))
      e <- e[length(e)]
      elem_ox <- c(elem_ox, e)
    }
  }
  elem_ox <- unique(elem_ox)

  # generate a formatted list of elements, oxidation states, and their EQ3 element names
  redox_elem_states <- list()
  for(elem in suppress_redox){
    # check for multiple oxidation states for this element in data0_params$formula_ox
    elem_ox_match <- unlist(lapply(lapply(makeup(elem_ox), FUN=names), `[`, 1)) # remove 'Z' in makeup()
    redox_entry <- c()

    matches <- elem_ox[elem_ox_match == elem]

    for(match in matches){
      if(grepl("\\+", match)){
        mag <- tolower(as.roman(as.numeric(strsplit(match, "\\+")[[1]][2])))
        if(is.na(mag)){
          mag <- "i"
        }
        redox_entry[match] <- paste0(elem, mag)
      }else if(grepl("\\-", match)){
        mag <- tolower(as.roman(as.numeric(strsplit(match, "\\-")[[1]][2])))
        if(is.na(mag)){
          mag <- "i"
        }
        redox_entry[match] <- paste0(elem, mag, "n") # n for negative charge
      }else{
        redox_entry[match] <- paste0(elem, "z") # z for zero charge
      }
    }
    redox_elem_states[[elem]] <- redox_entry

  }

  # Assign Americium a mass of 223 in the CHNOSZ element database.
  e <- thermo()$element
  e[e[,"element"] == "Am", "mass"] <- 223
  thermo(element = e)

  # Add oxidation-separated elements to CHNOSZ's database of elements
  # element, state, source, mass, s, n
  for(elem in names(redox_elem_states)){
    for(elem_ox in names(redox_elem_states[[elem]])){
      old <- thermo()$element
      elem_entry <- filter(thermo()$element, element==elem)[1, ] # todo: more elegant solution for getting one matching entry, e.g. in the case of matching Cl
      warning(paste(elem_entry))
      Xprops <- data.frame(element=unname(redox_elem_states[[elem]][elem_ox]),
                          state=elem_entry[, "state"],
                          source=paste("redox suppression workaround for", elem_ox),
                          mass=elem_entry[, "mass"],
                          s=elem_entry[, "s"],
                          n=elem_entry[, "n"])
      new <- rbind(old, Xprops)
      thermo(element=new)
    }
  }

  # load thermodynamic data
  custom_obigt <- read.csv(filename, stringsAsFactors=F)
  rownames(custom_obigt) <- custom_obigt$name


  # check that data0 supp and obigt file have same number of rows and same order.
  stopifnot(nrow(data0_params)==nrow(custom_obigt) && data0_params["name"]==custom_obigt["name"])


  known_species <- to_vec(for(i in 1:length(known_oxstates)) paste0(names(known_oxstates)[i], known_oxstates[i]))


  if(infer_formula_ox){

    makeup_list <- makeup(custom_obigt$formula)
    names(makeup_list) <- custom_obigt$name

    result <- lapply(makeup_list, FUN=get_oxstate)

    inferred_formula_ox <- c()

    for(entry in names(result)){
      if(!is.na(result[[entry]])){
        ox <- result[[entry]]
        elem <- names(ox)
        n_elem <- makeup_list[[entry]][elem]
        if(n_elem == 1){
          n_elem <- ""
        }
        out_elem <- paste0(n_elem, elem, format_charge(as.numeric(ox)))
      
        elem_other <- names(makeup_list[[entry]])[!(names(makeup_list[[entry]]) %in% c("Z", elem))]
        out <- to_vec(for(x in elem_other) if(makeup_list[[entry]][x] != 1) paste0(makeup_list[[entry]][x], x, known_oxstates[x]) else(paste0(x, known_oxstates[x])))
        inferred_formula_ox[[entry]] <- paste(out_elem, paste(out, collapse=" "))

        data0_params[entry, "formula_ox"] <- inferred_formula_ox[[entry]]
      } else if (entry %in% known_species){
        data0_params[entry, "formula_ox"] <- entry
      }
    }
    write.csv(data0_params, data0_formula_ox_name, row.names=F, na="")
  }

  modified_custom_obigt <- custom_obigt

  if(length(suppress_redox) > 0){

    # loop through species
    for(idx in 1:nrow(custom_obigt)){

      species_name <- custom_obigt[idx, "name"]

      # get makeup and charge of this species
      this_makeup <- makeup(custom_obigt[idx, "formula"])
      elems <- names(this_makeup)
      charge <- this_makeup["Z"]
      if(is.na(charge)){
        charge <- 0
      }
        
      # for each element in this species, see if it matches an entry in redox_elem_states
      this_formula <- c()
      for(elem in elems){
        
        if(elem %in% names(redox_elem_states)){

          # get oxidation state info from obigt data0 add-on file
          formula_ox <- strsplit(data0_params[species_name, "formula_ox"], " ")
          names(formula_ox) <- species_name

          # assign abundances to each ox elem
          for(elem_ox in formula_ox[[species_name]]){

            # determine the formula coefficient
            this_coeff <- as.numeric(gsub("([0-9]?)[A-Z].*", "\\1", elem_ox, perl=T))
            if(is.na(this_coeff)){
              this_coeff <- 1
            }

            # determine the element
            e <- unlist(strsplit(elem_ox, "[0-9](?=[A-Z])", perl=T))
            e <- e[length(e)]

            this_formula[e] <- this_coeff
          }
        }
      }

      modified_formula <- c()
      for(name in names(this_formula)){
        elem <- names(which(sapply(lapply(redox_elem_states, FUN=names), FUN=function(x) name %in% x)))
        if(!identical(elem, character(0))){
          modified_formula[redox_elem_states[[elem]][name]] <- this_formula[name]
        }else{
          e_makeup <- names(makeup(name))
          modified_formula[e_makeup[!(e_makeup %in% c("Z"))]] <- this_formula[name]
        }
      }

      formula_vec <- c(rbind(names(modified_formula), modified_formula))
      formula <- paste(formula_vec[formula_vec!="1"], collapse="")

      if(formula==""){
        formula <- custom_obigt[idx, "formula"]
      }else{
        formula <- paste0(formula, format_charge(charge))
      }
      modified_custom_obigt[species_name, "formula"] <- formula
    }
  }
                                   
  modified_custom_obigt %>% mutate_if(is.factor, as.character) -> modified_custom_obigt

  # include modified file in CHNOSZ database
  suppressMessages({
    db_idx <- mod.OBIGT(modified_custom_obigt) # produces a message
    j_df <- modified_custom_obigt %>% left_join(data0_params, join_by=name) # produces a message
  })

  j_df_basis <- filter(j_df, tag=="basis")

  basis_elem <- unlist(lapply(lapply(lapply(j_df_basis$formula, makeup), names), function(x) x[!(x %in% c("Z", "O", "H"))]))

  basis_pref <- unlist(j_df_basis["name"])
  names(basis_pref) <- basis_elem


  # EQ3 has Cl-, H2O, and O2(g) hard-coded as basis species for the
  # elements Cl, H, and O, respectively.
  basis_pref["Cl"] <- "Cl-"
  basis_pref["H"] <- "H2O"
  basis_pref["O"] <- "O2(g)"
  
  # handle the rest of the preferred basis species, e.g.,
  # define basis species for redox pseudoelements
  # basis_pref["Feiii"] <- "Fe+3"
  # basis_pref["Feii"] <- "Fe+2"
  # basis_pref["Fez"] <- "iron"
  # basis_pref["Siin"] <- "HS-"
  # basis_pref["Svi"] <- "SO4-2"
  for(i in 1:length(basis_pref_names)){
    basis_pref[basis_pref_names[i]] <- basis_prefs[i]
  }

  # specify molecules to balance H, O, and charge (Z)
  HOZ_balancers <- c("H+", "O2(g)", "H2O") # might be data0-specific (e.g., "O2(g)")

  # generate dissociation reactions
  dissrxns <- suppressMessages(get_dissrxn(sp_name=unlist(j_df["name"]),
                                           basis_pref=basis_pref,
                                           HOZ_balancers=HOZ_balancers,
                                           db_idx=db_idx))

  create_data0(supp_file,
               supp_file_ss,
               grid_temps,
               grid_press,
               db,
               template,
               dissrxns,
               modified_custom_obigt,
               db_idx,
               basis_pref,
               exceed_Ttr,
               verbose)

  if(generate_template){
    # create a template for sample input
    input_template <- data.frame(Sample=c("id"), `H+`=c("pH"), Temperature=c("degC"), O2=c("Molality"), check.names=F)
    for(basis in dissrxns[["basis_list"]]){
      input_template[[basis]] <- c("Molality")
    }
    write.csv(input_template, template_name, row.names=F)
  }
}