# load packages
suppressMessages({
  library(CHNOSZ)
  library(dplyr)
  library(comprehenr)
  library(stringr)
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

# function to round up the last digit of a four-digit number
# e.g., 52.6112 becomes 52.6113.
# does not work if the last digit is 9... TODO
# This function is useful for reporting rounded PSAT pressures that keeps water in the liquid phase
roundup <- function(x, n){
  endval <- as.numeric(str_sub(format(round(x, n), nsmall = n, scientific=F),start=-1))+1
  firstvals <- substr(x, start = 1, stop = nchar(format(round(x, n), nsmall = n, scientific=F))-1)
  return(paste0(firstvals, endval))
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
sd <- function(x, k) trimws(format(round(x, k), nsmall=k, scientific=F))


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
# thermo_df has to be the entire thermo datasheet.
#   Needed to determine whether aux species are better than strict basis depending on ox state.
#   If NULL, just uses strict basis species, or what it determines to be strict basis.

get_dissrxn <- function(sp_name, redox_elem_states, basis_pref=c(), aux_pref=c(),
                        HOZ_balancers=c("H+", "O2", "H2O"),
                        thermo_df=NULL, verbose=2){
  
  if(length(sp_name) > 0){
    # get a vector of elements that make up the (non-basis) species
    basis_elem <- (function (x) setdiff(names(unlist(makeup(info(info(x), check.it=F)$formula))), c("H", "O", "Z"))) (rownames(thermo_df))
    basis_elem <- c(basis_elem, names(basis_pref))
  }else{
    # if there are no species, define basis list in dissrxns and exit the function
    dissrxns <- list(basis_list = basis_pref)
    return(dissrxns)
  }
  basis_elem <- unique(basis_elem)
                   
  sp <- thermo_df
  sp_formula <- sp$formula_modded
  sp_formula_makeup <- makeup(sp_formula)
  names(sp_formula_makeup) <- sp$name

  # loop through each element and assign preferred basis species (if supplied)
  # or automatically choose compatible basis species
  basis_list <- list()

  # check that every basis element has a basis pref
  if(!setequal(basis_elem, names(basis_pref))){
    missing_basis <- basis_elem[!(basis_elem %in% names(basis_pref))]
    stop(paste("Error: the element(s)", paste(missing_basis, collapse=","), "require strict basis species in the database."))
  }

  for(elem in basis_elem){
    # if a preferred basis species exists for this element, assign it and move to next element
    if(elem %in% names(basis_pref)){
      basis_list[[elem]] <- c(info(info(basis_pref[elem]), check.it=F)$name)
      names(basis_list[[elem]]) <- basis_pref[elem]
      next
    }

    # get index of species that fulfill requirements of EQ3's definition of a basis species for this element
    idx <- unlist(lapply(lapply(sp_formula_makeup, names), match_basis_comp, elem=elem))

    # provide a vector of potential basis species for this element
    basis_list[[elem]] <- sp$name[idx]
    names(basis_list[[elem]]) <- sp$name[idx]

    # get formula makeup of potential basis species
    frm_mkup <- lapply(unlist(lapply(lapply(lapply(basis_list[[elem]], FUN=info), FUN=info), `[`, "formula")), makeup)
        
    # if possible, narrow down this list to basis species with a matching element abundance of 1
    abund1 <- unlist(lapply(frm_mkup, function(x) x[elem]==1))
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
  simplest_basis <- c()
  for(elem in basis_elem){
    form <- lapply(info(info(basis_list[[elem]]))$formula, makeup)
    sums <- unlist(lapply(lapply(form, abs), sum))
    idx <- which(sums==min(sums))
    if(length(idx) > 1){
      idx <- idx[1]
    }
    simplest_basis[elem] <- basis_list[[elem]][idx]
  }

  dissrxns <- lapply(sp_name, spec_diss, simplest_basis, sp_formula_makeup,
                     HOZ_balancers, redox_elem_states, aux_pref, thermo_df, verbose)
  names(dissrxns) <- sp_name
  dissrxns[["basis_list"]] <- simplest_basis
  return(dissrxns)

}
                            
# supply a species name, get a dissrxn.
# (meant to lapply this over all species names)
spec_diss <- function(sp, simplest_basis, sp_formula_makeup, HOZ_balancers,
                      redox_elem_states, aux_pref, thermo_df=NULL, verbose=2){
 
  # determine a balanced dissociation reaction into basis species

  # get the elemental composition
  comp <- makeup(thermo_df[thermo_df["name"]==sp, "formula_modded"])
    
  # get the elements that are not H, O, Z (charge)
  basis_elem <- setdiff(names(comp), c("H", "O", "Z"))
    
  # Instantiate a list of chosen basis species and their coeffs in the dissrxns.
  # Start with strict basis species only. E.g., if SO4-2 is a strict basis species:
  # Pb(HS)3-, chosen_basis$Pb is Pb+2:1 and chosen_basis$S is SO4-2:3
  # In the code below, auxiliary basis species can be chosen instead:
  # Pb(HS)3-, chosen_basis$Pb is Pb+2:1 and chosen_basis$S is HS-:3
  chosen_basis <- list()
    
  # Check whether an auxiliary basis species given in thermo_df would work
  # better than a strict basis species.
    
  if(!is.null(thermo_df)){
      
    sp_formula_ox = thermo_df[thermo_df["name"]==sp, "formula_ox_modded"]
      
    for(elem in basis_elem){
        
      sp_formula_ox_split <- strsplit(sp_formula_ox, " ")
    
#       print(paste("species:", sp))
#       print(paste("elem:", elem))
#       print(paste("formula ox:", sp_formula_ox))
#       print("formula ox split:")
        
      # get the element plus oxidation state from formula ox
      # e.g., turn "2Al+3" into "Al+3"
      # or "0.165Mg+2" to "Mg+2"
      step1 <- sp_formula_ox_split[[1]]
      step2 <- step1[grepl(paste0(elem,"([^a-z]|$)"), step1)]
      step3 <- gsub("^[[:digit:]]+(\\.[[:digit:]]+)?", "", step2)
      sp_ox_elem <- makeup(step3)
        
      sp_num_elem <- as.numeric(gsub("([0-9]?)[A-Z].*", "\\1", sp_formula_ox_split[[1]][grepl(paste0(elem,"([^a-z]|$)"), sp_formula_ox_split[[1]], perl=T)], perl=T))
      sp_num_elem <- ifelse(is.na(sp_num_elem), 1, sp_num_elem)
      names(sp_num_elem) <- rep(simplest_basis[[elem]], length(sp_num_elem))
      n_elem_in_basis <- sp_formula_makeup[[simplest_basis[elem]]][elem]
      chosen_basis[[elem]] <- sp_num_elem/n_elem_in_basis
        
      if(elem %in% names(simplest_basis) && elem %in% names(aux_pref)){
        # get oxidation states of this element available in strict and aux basis species
        strict_df = thermo_df %>%
          filter(name == simplest_basis[elem]) %>%
          select(formula_ox_modded)
        aux_df = thermo_df %>%
          filter(name %in% aux_pref[names(aux_pref)==elem]) %>%
          filter(name != sp) %>% # prevents an aux species from dissociating into itself
          select(formula_ox_modded)
        basis_df <- rbind(strict_df, aux_df)
        basis_df["ave_ox_state_of_elem"] <- 0
        basis_formula_ox_split <- lapply(basis_df["formula_ox_modded"], strsplit, " ")$formula_ox_modded
          
        idx <- 1
        for(line in basis_formula_ox_split){
          basis_num_elem <- as.numeric(gsub("([0-9]?)[A-Z].*", "\\1", line[grepl(paste0(elem,"([^a-z]|$)"), line, perl=T)], perl=T))
          basis_num_elem <- ifelse(is.na(basis_num_elem), 1, basis_num_elem)
          basis_ox_elem <- makeup(gsub("^[[:digit:]]+(\\.[[:digit:]]+)?", "", line[grepl(paste0(elem,"([^a-z]|$)"), line, perl=T)]))
          if(length(basis_num_elem) > 1){
            basis_oxstates <- c()
            basis_total_atoms <- c()
            for(i in 1:length(basis_num_elem)){
              if(!("Z" %in% names(basis_ox_elem[[i]]))){
                basis_ox_elem[[i]]["Z"] <- 0
              }
              basis_oxstates <- c(basis_oxstates, basis_ox_elem[[i]]["Z"]*basis_num_elem[i])
              basis_total_atoms <- c(basis_total_atoms, basis_num_elem[i])
            }
            basis_ave_ox_state_of_elem <- sum(basis_oxstates)/sum(basis_total_atoms)
          }else{
            if(!("Z" %in% names(basis_ox_elem))){
              basis_ave_ox_state_of_elem <- 0
            }else{
              basis_ave_ox_state_of_elem <- basis_ox_elem["Z"]
            }

          }
          basis_df[idx, "ave_ox_state_of_elem"] <- basis_ave_ox_state_of_elem
          idx <- idx + 1
        }
          
        # get the oxidation state of the element in the non-basis species
        sp_formula_ox = thermo_df[thermo_df["name"]==sp, "formula_ox_modded"]
        
        # For each oxidation state of this element in the non-basis species,
        # pick the closest strict or aux basis species.
        if(is.list(sp_ox_elem)){
          # When there is more than one instance of this element in the non-basis species,
          # find the basis species with the closest oxidation states and assign.
          chosen_basis_species_vec <- c()
          basis_ave_ox_states <- unlist(basis_df["ave_ox_state_of_elem"]) 
        
          names(basis_ave_ox_states) <- row.names(basis_df)
          for(i in 1:length(sp_ox_elem)){
            if(!("Z" %in% names(sp_ox_elem[[i]]))){
              sp_ox_elem[[i]]["Z"] <- 0
            }
            chosen_basis_species_vec <- c(chosen_basis_species_vec, names(which.min(abs(basis_ave_ox_states-sp_ox_elem[[i]]["Z"])))[1])
          }
          
          basis_num_elems <- c()
          for(basis in chosen_basis_species_vec){
            n_elem_in_basis <- sp_formula_makeup[[basis]][elem]
            basis_num_elems <- c(basis_num_elems, n_elem_in_basis)
          }
          names(sp_num_elem) <- chosen_basis_species_vec
          chosen_basis[[elem]] <- sp_num_elem/basis_num_elems
          
            
        }else{
          # Find the basis species with the closest oxidation state to the
          # element inside the non-basis species and assign.
          if(!("Z" %in% names(sp_ox_elem))){
            sp_ox_elem["Z"] <- 0
          }
          basis_ave_ox_states <- unlist(basis_df["ave_ox_state_of_elem"])
          names(basis_ave_ox_states) <- row.names(basis_df)
          
          # if there is more than one basis species matching the oxidation state of
          # this element in the non-basis species, just pick the first one, thus the [1].
          # Maybe there is a better way to decide but I'm going with this for now.
          chosen_basis_name <- names(which.min(abs(basis_ave_ox_states-sp_ox_elem["Z"])))[1]
          
          n_elem_in_basis <- sp_formula_makeup[[chosen_basis_name]][elem]
            
          names(sp_num_elem) <- chosen_basis_name
            
          chosen_basis[[elem]] <- sp_num_elem/n_elem_in_basis
            
        }
        
      }

      # sum together basis species coeffs sharing a name (e.g., c(HS-:1, HS-:1) is summed to c(HS-:2))
      chosen_basis[[elem]] <- tapply(unlist(chosen_basis[[elem]]), names(unlist(chosen_basis[[elem]])), sum)
    }
  }
   
#   print("Chosen basis")
#   print(chosen_basis)
    
  # set basis species in CHNOSZ so that a dissociation reaction can be autobalanced
  basis(unique(c(simplest_basis[basis_elem], HOZ_balancers)))

  basis_to_include <- unlist(lapply(chosen_basis, names))
  coeffs_to_include <- unlist(lapply(chosen_basis, `[`))

  isp <- suppressMessages(info(c(sp, basis_to_include)))

  autobalance <- subcrt_bal(isp, c(-1, coeffs_to_include))
  spec_names <- thermo()$OBIGT[autobalance$newspecies, "name"]
  spec_names[spec_names == "water"] <- "H2O" # prevent CHNOSZ from renaming 'H2O' (name req. by EQ3) to 'water'
  sp_dissrxn <- autobalance$newcoeff
  names(sp_dissrxn) <- spec_names
    
  sp_dissrxn <- paste(c(rbind(sprintf("%.4f", round(sp_dissrxn, 4)), names(sp_dissrxn))), collapse=" ")
    
  return(sp_dissrxn)
}

######################################## create_data0 function

# main function
create_data0 <- function(thermo_df,
                         filename_ss=NULL,
                         grid_temps,
                         grid_press,
                         db,
                         water_model,
                         template,
                         dissrxns,
                         db_idx,
                         basis_pref=c(),
                         exceed_Ttr=FALSE,
                         verbose){

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
      dissrxn_in_datafile <- thermo_df[i, "dissrxn"]
      species_name <- thermo_df[i, "name"]
      
      if(dissrxn_in_datafile != "" && is.null(dissrxns[[species_name]])){
        dissrxn <- dissrxn_in_datafile
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
        dissrxn_temp <- dissrxn[3:length(dissrxn)] # ignore the "-1" and diss. spec. name
        dissrxn_names <- dissrxn_temp[c(FALSE, TRUE)] # get names of reactants and products
        dissrxn_coefs <- dissrxn_temp[c(TRUE, FALSE)] # get coeffs of reactants and products
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
  for(idx in db_idx){

    entry <- suppressMessages(info(idx))
    name <- entry$name

    if (name == "O2(g)"){
      vmessage("O2(g) is included as a basis species by default. Moving to the next species...", 2, verbose)
      next
    }
      
    date <- entry$date
    date <- strtrim(date, 9) # truncates date if greater than 9 letters
    elem <- makeup(idx) # get elemental composition

    aux_basis <- FALSE
    if(name %in% lapply(dissrxns[["basis_list"]], `[[`, 1)){
      # if this species is marked as a preferred basis species, move to the next species
      vmessage(paste0("'", name, "' (basis species) processed successfully."), 2, verbose)
      next
    }else if(thermo_df[which(thermo_df["name"]==name), "tag"] == "basis"){
      # if this is marked as a basis in the data0 supplemental file, move to the next species
      vmessage(paste0("'", name, "' (basis species) processed successfully."), 2, verbose)
      next
    }else if(thermo_df[which(thermo_df["name"]==name), "tag"] == "aux"){
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
    logK_grid <- rep(0, 8)
    tryCatch({
      # the subcrt() calculation for each P-T in the grid
      if(thermo_df[entry$name, "tag"] != "basis"){
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
    if(basis == "O2(g)" | basis == "H2O" | basis == "Cl-" | basis == "H+" | basis == "water"){
      next
    }

    # get the date on the basis species entry
    suppressMessages({
      date_basis <- info(info(basis))$date
    })

    # get the basis species formula
    basis_formula <- add_obigt_df[which(add_obigt_df[, "name"] == basis), "formula_modded"]

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
      if(add_obigt_df[i, "name"] == "Cl-" | add_obigt_df[i, "name"] == "O2" | add_obigt_df[i, "name"] == "OH-"){
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


  # calculate debye huckel a and b parameters for the grid
  A_DH_grid <- unlist(water("A_DH", T=273.15+grid_temps, P=grid_press))
  B_DH_grid <- unlist(water("B_DH", T=273.15+grid_temps, P=grid_press)*10^-8)

  # format grid values
  grid_temps_f <- as.character(format(round(grid_temps, 4), nsmall = 4, scientific=F))
  grid_press_f <- as.character(format(round(grid_press, 4), nsmall = 4, scientific=F))
  A_DH_grid_f <- as.character(format(round(A_DH_grid, 4), nsmall = 4, scientific=F))
  B_DH_grid_f <- as.character(format(round(B_DH_grid, 4), nsmall = 4, scientific=F))

  bdot_grid_f <- as.character(format(round(calc_bdot(grid_temps), 4), nsmall = 4, scientific=F))


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

  # insert minimum and maximum temperature values into data0 template
  temp_min_max_insertlines <- "\nTemperature limits \\(degC\\)\n.*\ntemperatures\n"
  t_min <- min(grid_temps)
  t_max <- max(grid_temps)
  t_min_f <- as.character(format(round(t_min, 4), nsmall = 4, scientific=F))
  t_max_f <- as.character(format(round(t_max, 4), nsmall = 4, scientific=F))
  t_min_max <- paste0(paste(rep(" ", 10-nchar(t_min_f)), collapse=""), t_min_f)
  t_min_max <- paste0(t_min_max, paste(rep(" ", 10-nchar(t_max_f)), collapse=""), t_max_f)
  t_min_max <- paste0("     ", t_min_max)
  data0_template <- sub(temp_min_max_insertlines, paste0("\nTemperature limits (degC)\n", t_min_max, "\ntemperatures\n"), data0_template)

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
  desc <- "data0.%s\nWater model: %s"
  min_desc <- "data0.min\nminimal working data0 file"
  data0_template <- sub(min_desc, sprintf(desc, db, water_model), data0_template)
  vmessage("Finished.", 2, verbose)
  write(data0_template, paste0("data0.", db))
}


######################################## Main function
                            
main_create_data0 <- function(filename,
                              filename_ss,
                              grid_temps,
                              grid_press,
                              db,
                              water_model,
                              template,
                              exceed_Ttr,
                              data0_formula_ox_name,
                              suppress_redox,
                              infer_formula_ox,
                              generate_template,
                              template_name,
                              template_type,
                              verbose){
  
  # set water model
  suppressMessages(water(water_model))
  
  # round grid temperatures to four decimal places
  grid_temps <- round(grid_temps, 4)
    
  # check that water is a liquid at each T-P point
  if(length(grid_press) > 1){
    TP_grid_errors <- c()
    for(i in 1:length(grid_temps)){
      tryCatch({
        psat_press <<- suppressMessages(water("Psat", T=grid_temps[i]+273.15)[[1]])
      }, error=function(e){
        psat_press <<- NA
      })
      if(is.na(psat_press)){
        #pass
      }else if(grid_press[i] < psat_press){
        TP_grid_errors <- c(TP_grid_errors,
                            paste("\n", grid_press[i], "bar is below liquid-vapor",
                                  "saturation pressure", roundup(psat_press, 4),
                                  "bar at", grid_temps[i], "degrees C."))
      }
    }
    if(length(TP_grid_errors) > 0){
      stop(paste(paste(TP_grid_errors, collapse="\n"),
                 "\n\nIncrease the pressure at these temperature points in 'grid_P'",
                 "to keep water in a liquid state."))
    }
  }
    
  # calculate PSAT pressure if specified by user or if pressure grid
  # has a number of values that does not equal temperature grid length.
  if(tolower(grid_press) == "psat"){
      vmessage("Calculating pressure grid along liquid-vapor saturation curve...", 2, verbose)
      grid_press <- water("Psat", T=grid_temps+273.15)[[1]]
  }
    
  # check that pressure polynomials can be calculated with temperature grid
  grid_temps <- as.numeric(grid_temps)
  grid_press <- as.numeric(grid_press)

  # third order polynomial for the first T-P range
  poly_coeffs_1 <- lm(grid_press[1:4] ~ poly(grid_temps[1:4], 3, raw=T))$coefficients
  
  # fourth order polynomial for the second T-P range
  poly_coeffs_2 <- lm(grid_press[4:8] ~ poly(grid_temps[4:8], 4, raw=T))$coefficients
  
  for(T in grid_temps[1:4]){
    if(is.na(poly_coeffs_1[1] + poly_coeffs_1[2]*T + poly_coeffs_1[3]*T^2 + poly_coeffs_1[4]*T^3)){
      stop(paste0("Error: Could not compute the coefficients of an interpolating polynomial
                   for the first four values of the temperature grid: [", paste(grid_temps[1:4], collapse=", "), "]."))
    }
  }
  for(T in grid_temps[4:8]){
    if(is.na(poly_coeffs_2[1] + poly_coeffs_2[2]*T + poly_coeffs_2[3]*T^2 + poly_coeffs_2[4]*T^3)){
      stop(paste0("Error: Could not compute the coefficients of an interpolating polynomial
                   for the last five values of the temperature grid: [", paste(grid_temps[4:8], collapse=", "), "]."))
    }
  }
    
  # specify molecules to balance H, O, and charge (Z)
  HOZ_balancers <- c("H+", "O2(g)", "H2O") # might be dataset-specific (e.g., "O2(g)")
  
  # load thermodynamic data
  thermo_df <- read.csv(filename, stringsAsFactors=F)
  rownames(thermo_df) <- thermo_df$name
    
  # TODO: ensure that Cl- (and perhaps other hard-coded species) are in thermo_df
  # and return an error if not.
    
  # get all unique oxidation states of elements in the entire dataset
  elem_ox <- c()
  for(entry in thermo_df$formula_ox){
    entry <- strsplit(entry, " ")[[1]]
    for(elem in entry){
      e <- unlist(strsplit(elem, "[0-9](?=[A-Z])", perl=T))
      e <- e[length(e)]
      elem_ox <- c(elem_ox, e)
    }
  }
  elem_ox <- unique(elem_ox)

  # create pseudoelement names
  redox_elem_states <- list()
  for(elem in suppress_redox){
    # check for multiple oxidation states for this element in thermo_df$formula_ox
    elem_ox_match <- unlist(lapply(lapply(makeup(elem_ox), FUN=names), `[`, 1)) # remove 'Z' in makeup()
    redox_entry <- c()

    matches <- elem_ox[elem_ox_match == elem]

    for(match in matches){
      if(grepl("\\+", match)){
        mag <- tolower(as.roman(as.numeric(strsplit(match, "\\+")[[1]][2])))
        if(is.na(mag)){
          mag <- "i"
        }
        redox_entry[match] <- paste0(elem, mag, "p")
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
      Xprops <- data.frame(element=unname(redox_elem_states[[elem]][elem_ox]),
                          state=elem_entry[, "state"],
                          source=paste("redox suppression workaround for", elem_ox),
                          mass=elem_entry[, "mass"],
                          s=elem_entry[, "s"],
                          n=elem_entry[, "n"])
      new <- rbind(old, Xprops)
      suppressMessages(thermo(element=new))
    }
  }

  known_species <- to_vec(for(i in 1:length(known_oxstates)) paste0(names(known_oxstates)[i], known_oxstates[i]))

  if(infer_formula_ox){
      
    makeup_list <- makeup(thermo_df$formula)
    names(makeup_list) <- thermo_df$name

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

        thermo_df[entry, "formula_ox"] <- inferred_formula_ox[[entry]]
      } else if (entry %in% known_species){
        thermo_df[entry, "formula_ox"] <- entry
      }
    }
    write.csv(thermo_df, data0_formula_ox_name, row.names=F, na="")
  }
                      
  thermo_df["formula_modded"] <- thermo_df["formula"]
  thermo_df["formula_ox_modded"] <- thermo_df["formula_ox"]
                      
  # create pseudoelements and assign to molecular formulae
  if(length(suppress_redox) > 0){

    # loop through species
    for(idx in 1:nrow(thermo_df)){

      species_name <- thermo_df[idx, "name"]

      # get makeup and charge of this species
      this_makeup <- makeup(thermo_df[idx, "formula"])
      elems <- names(this_makeup)
      charge <- this_makeup["Z"]
      if(is.na(charge)){
        charge <- 0
      }
        
      # for each element in this species, see if it matches an entry in redox_elem_states
      this_formula <- c()
      for(elem in elems){
        
        if(elem %in% names(redox_elem_states)){
            
          sp_formula_ox <- thermo_df[species_name, "formula_ox"]

          # get oxidation state info
          formula_ox <- strsplit(sp_formula_ox, " ")
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
        
          # create a modified formula_ox column in the modified_thermo_df if there
          # is redox suppression, such that "Fe+2 2Cl-" becomes "Feii+2 2Cl-"
          this_formula_modded <- this_formula
          for(ox_elem in names(redox_elem_states)){
            for(ox_state in names(this_formula)){
              if(grepl(paste0(ox_elem,"([^a-z]|$)"), ox_state)){
                  
                # handle positive oxidation states
                splitname <- strsplit(ox_state, "\\+")[[1]]
                if(length(splitname) > 1){
                  ox_charge <- "+"
                }else{
                  # handle negative oxidation states
                  splitname <- strsplit(ox_state, "\\-")[[1]]
                  if(length(splitname) > 1){
                    ox_charge <- "-"
                  }else{
                    ox_charge <- ""
                  }
                }
                
                if(length(splitname) != 1){
                  ox_charge_magnitude <- as.numeric(splitname[2])
                }else{
                  ox_charge_magnitude <- 0
                }

                pseudoelem_ox_charge <- redox_elem_states[[ox_elem]][ox_state]
                if(ox_charge_magnitude !=0){
                  pseudoelem_ox_charge <- paste0(pseudoelem_ox_charge, ox_charge, ox_charge_magnitude)
                }
                
                this_formula_names <- names(this_formula_modded)
                this_formula_names[this_formula_names == ox_state] <- pseudoelem_ox_charge
                names(this_formula_modded) <- this_formula_names
                
              }
            }
          }
            
          to_collapse <- c()
          for(x in names(this_formula_modded)){
            formula_ox_vec <- c(rbind(this_formula_modded[x], x))
            this_formula_ox_modded <- paste(formula_ox_vec[formula_ox_vec!="1"], collapse="")
            to_collapse <- c(to_collapse, this_formula_ox_modded)
          }
          this_formula_ox_modded <- paste(to_collapse, collapse=" ")

          thermo_df[species_name, "formula_ox_modded"] <- this_formula_ox_modded
            
        }
      }
        
      modified_formula <- c()
      previous_elems <- c()
      for(name in names(this_formula)){
        elem <- names(which(sapply(lapply(redox_elem_states, FUN=names), FUN=function(x) name %in% x)))
        if(!identical(elem, character(0))){
          modified_formula[redox_elem_states[[elem]][name]] <- this_formula[name]
        }else{
          e_makeup <- names(makeup(name))
          elem <- e_makeup[!(e_makeup %in% c("Z"))]
          if(elem %in% previous_elems){
            modified_formula[elem] <- modified_formula[elem] + this_formula[name]
          }else{
            modified_formula[elem] <- this_formula[name]
          }
        }
        previous_elems <- c(previous_elems, elem)
      }
                                   
      formula_vec <- c(rbind(names(modified_formula), modified_formula))
      formula <- paste(formula_vec[formula_vec!="1"], collapse="")
                                   
      if(formula==""){
        formula <- thermo_df[idx, "formula_modded"]
      }else{
        formula <- paste0(formula, format_charge(charge))
      }
                                   
      thermo_df[species_name, "formula_modded"] <- formula
    }
  }
                                   
  thermo_df %>% mutate_if(is.factor, as.character) -> thermo_df
                                   
  # include modified file in CHNOSZ database
  to_mod_OBIGT <- thermo_df[c("name", "abbrv", "formula_modded",
                              "state", "ref1", "ref2", "date",
                              "E_units", "G", "H", "S", "Cp",
                              "V", "a1.a", "a2.b", "a3.c",
                              "a4.d", "c1.e", "c2.f",
                              "omega.lambda", "z.T")]
                                   
  names(to_mod_OBIGT) <- c("name", "abbrv", "formula",
                            "state", "ref1", "ref2", "date",
                            "E_units", "G", "H", "S", "Cp",
                            "V", "a1.a", "a2.b", "a3.c",
                            "a4.d", "c1.e", "c2.f",
                            "omega.lambda", "z.T")
                                   
  suppressMessages({
    db_idx <- mod.OBIGT(to_mod_OBIGT, replace=TRUE) # produces a message
  })

                                   
  # begin handling basis preferences
  basis_df <- thermo_df %>%
    filter(tag=="basis")
  aux_df <- thermo_df %>%
    filter(tag=="aux")

  basis_pref <- basis_df[, "name"]
  names(basis_pref) <- lapply(lapply(lapply(basis_df[,"formula_modded"], makeup), names), setdiff, c("Z", "O", "H"))
  aux_pref <- unlist(aux_df[, "name"])
  aux_pref_names <- lapply(lapply(lapply(aux_df[,"formula_modded"], makeup), names), setdiff, c("Z", "O", "H"))
  
  # Remove aux basis species from the preferred list if they contain more than
  # one element besides O and H.
  #     E.g., remove CN-, OCN-, etc.
  keep <- c()
  for(i in 1:length(aux_pref)){
    if(length(aux_pref_names[[i]])==1){
      keep <- c(keep, i)
    }
  }
  aux_pref <- aux_pref[keep]
  aux_pref_names <- aux_pref_names[keep]
    
  # Also remove aux basis species from the preferred list if they have more than
  # one atom of the same element.
  #     E.g., remove S2-2, S2O3-2, etc.
  keep <- c()
  for(i in 1:length(aux_pref)){
    if(makeup(thermo_df[thermo_df["name"]==aux_pref[[i]], "formula_modded"])[aux_pref_names[[i]]] == 1){
      keep <- c(keep, i)
    }
  }
  aux_pref <- aux_pref[keep]
  aux_pref_names <- aux_pref_names[keep]
    
  names(aux_pref) <- aux_pref_names
                                   
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
                                   
  # Ensure that each element has a representative strict basis species.
  # If not, assign one from aux species.
  # If that fails, assign one from a valid non-basis species.
  # First check that each element has only one valid basis species.
  #   1. get a vector of ALL elements except for O and H
  all_elements <- lapply(lapply(lapply(thermo_df[,"formula_modded"], makeup), names), setdiff, c("Z", "O", "H"))
  all_elements <- unique(unlist(all_elements))
  #   2. check that all elements are represented by basis species
  elem_need_basis <- setdiff(all_elements, names(basis_pref))
  #   3. check auxiliary basis species to see if any can be used as strict basis for unrepresented elements
  #      and assign to basis_prefs. Remove from aux basis prefs.
  if(length(elem_need_basis) > 0){
    for(elem in elem_need_basis){
      if(elem %in% names(aux_pref)){
        valid_aux <- aux_pref[names(aux_pref) == elem]
        
        basis_pref <- c(basis_pref, valid_aux[1]) # add valid aux to basis prefs. The [1] just selects the first valid aux species if there are multiple. Selection method can probably be improved.
        aux_pref <- aux_pref[aux_pref != valid_aux]
        
        thermo_df[valid_aux, "tag"] <- "basis"

      }else{
        # Pick a basis species from valid non-basis...
        # The dissociation reaction generator will do this automatically.
      }
    }
  }

  # get names of species that need dissrxns:
  #  1. non-basis species in the datafile lacking a dissociation reaction
  #  2. non-basis species with incorrect or unbalanced dissociation reactions.
  #     Note: unbalanced dissociation reactions can occur when suppressing redox on element(s).
  #     (e.g., chalcopyrite, CuFeS2, which contains Fe+3, might be written to dissociate into Fe+2.
  #     This creates an unbalanced reaction if iron oxidation states are all different pseudoelements)
  
  # for each dissrxn, determine if unbalanced. If so, flag a thermo_df 'need_dissrxn' column.
  thermo_df[, "regenerate_dissrxn"] <- F
  for(species in rownames(thermo_df)){
    
    tag <- thermo_df[species, "tag"]
    if(tag != "basis"){
      dissrxn <- thermo_df[species, "dissrxn"]
      dissrxn <- strsplit(dissrxn, " ")[[1]] # split the rxn into coeffs and species
      dissrxn_names <- dissrxn[c(FALSE, TRUE)] # get names of reactants and products
      dissrxn_ispecies <- suppressMessages(info(dissrxn_names))
      dissrxn_coefs <- dissrxn[c(TRUE, FALSE)] # get coeffs of reactants and products
      dissrxn_coefs <- as.numeric(dissrxn_coefs) # convert coeffs from str to numeric
      tryCatch({
        subcrt_bal(dissrxn_ispecies, dissrxn_coefs)
      }, error=function(e){
        thermo_df[species, "regenerate_dissrxn"] <<- T # assign global variable with <<- because this is within the error function
      })
    }
  }
                                   
  df_needs_dissrxns <- thermo_df %>%
    filter(tag != "basis") %>%
    filter(regenerate_dissrxn == T)
      
  if(nrow(df_needs_dissrxns) > 0){
    if(length(suppress_redox) == 0){
      needs_dissrxns_message <- paste("Balanced dissociation reactions are missing for species:",
                                      paste(unlist(df_needs_dissrxns["name"]), collapse=", "))
    }else{
      needs_dissrxns_message <- paste("Balanced and/or redox-suppressed dissociation reactions are missing for species:",
                                      paste(unlist(df_needs_dissrxns["name"]), collapse=", "))
    }
      
    vmessage(needs_dissrxns_message, 1, verbose)
    vmessage("Generating dissociation reactions for these species using strict and auxiliary basis species containing a maximum of one atom of one element besides O and H...", 1, verbose)
  }
                                   
  # generate dissociation reactions
  dissrxns <- suppressMessages(get_dissrxn(sp_name=unlist(df_needs_dissrxns["name"]),
                                           basis_pref=basis_pref,
                                           aux_pref=aux_pref,
                                           HOZ_balancers=HOZ_balancers,
                                           thermo_df=thermo_df,
                                           verbose=verbose,
                                           redox_elem_states=redox_elem_states))
                                   
  # Produce a warning message about which dissrxns were (re)generated and what they are.
  if(nrow(df_needs_dissrxns) > 0){
    names <- df_needs_dissrxns[["name"]]
    generated_dissrxns <- dissrxns[names]
    nonbasis_idx <- unlist(lapply(lapply(lapply(generated_dissrxns, strsplit, " "), `[[`, 1), FUN=function(x) length(x)!=4))
    basis_idx <- !nonbasis_idx
    nonbasis_names <- names(nonbasis_idx)[nonbasis_idx]
    basis_names <- names(basis_idx)[basis_idx]
    vmessage(paste(nonbasis_names, ":", generated_dissrxns[nonbasis_names], "\n"), 1, verbose)
    vmessage(paste("Species that have been converted into strict basis:", paste(basis_names, collapse=", ")), 1, verbose)
  }

  create_data0(thermo_df,
               filename_ss,
               grid_temps,
               grid_press,
               db,
               water_model,
               template,
               dissrxns,
               db_idx,
               basis_pref,
               exceed_Ttr,
               verbose)
                                  
  if(generate_template){
    # create a template for sample input
      
    if(template_type == 'strict'){
      species_names <- thermo_df %>% filter(tag == "basis")
    }else if(template_type == 'all basis'){
      species_names <- thermo_df %>% filter(tag == "basis" | tag == "aux")
    }else{
      species_names <- thermo_df
    }

    species_names <- sort(unlist(species_names[, "name"]))
    input_template <- data.frame(Sample=c("id"), `H+`=c("pH"), Temperature=c("degC"), logfO2=c("logfO2"), check.names=F)
    for(species in species_names){
      if(!(species %in% c("O2(g)", "H2O", "water"))){
        input_template[[species]] <- c("Molality")
      }
    }
    write.csv(input_template, template_name, row.names=F)
  }
}