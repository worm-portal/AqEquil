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

# print messages if 'verbose' setting >= vlevel of message.
vmessage <- function(m, vlevel, verbose){
  if(verbose >= vlevel){
    print(m)
  }
}

# function for inserting a row at index r while shifting other rows down
insertRow <- function(existingDF, newrow, r) {
  if(r < nrow(existingDF)){
    existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
    return(existingDF)
  }else{
    return(rbind(existingDF, newrow))
  }
}

# checks if sp in a dissrxn all match strict basis names
all_in_strict_basis <- function(s, strict_basis_names, fixed_species){
  s <- s[c(FALSE, TRUE)]
  s <- s[2:length(s)]
  s <- s[!s %in% fixed_species]
  return(all(s %in% strict_basis_names))
}

# function to vectorize all_in_strict_basis()
split_and_check_dissrxn <- function(dissrxn, strict_basis_names, fixed_species){
  s <- str_split(dissrxn, " ")
  s <- unlist(lapply(s, all_in_strict_basis, strict_basis_names, fixed_species))
  return(s)
}

# order thermo_df to put strict basis species in the front, aux basis species in
# the middle (and sort them to prevent EQ3 errors), and everything else at the
# end.
order_thermo_df <- function(thermo_df, fixed_species, verbose){
    
  basis_entries <- filter(thermo_df, tag=="basis")
  aux_entries   <- filter(thermo_df, tag=="aux")
  other_entries <- filter(thermo_df, tag!="basis" & tag!="aux")
    
  aux_entries_with_only_strict <- aux_entries[split_and_check_dissrxn(aux_entries[, "dissrxn"], basis_entries[, "name"], fixed_species), ]
  aux_entries_with_other_aux   <- aux_entries[!split_and_check_dissrxn(aux_entries[, "dissrxn"], basis_entries[, "name"], fixed_species), ]
  aux_entries_with_strict_aux <- aux_entries_with_other_aux[split_and_check_dissrxn(aux_entries_with_other_aux[, "dissrxn"], c(aux_entries_with_only_strict[, "name"], basis_entries[, "name"]), fixed_species), ]
    
  if(!identical(aux_entries_with_other_aux[, "dissrxn"], character(0))){
    aux_entries_with_nonstrict_aux <- aux_entries_with_other_aux[!split_and_check_dissrxn(aux_entries_with_other_aux[, "dissrxn"], c(aux_entries_with_only_strict[, "name"], basis_entries[, "name"]), fixed_species), ]
  }else{
    aux_entries_with_nonstrict_aux <- data.frame()
  }

  if(nrow(aux_entries_with_nonstrict_aux) > 0){
    # sort aux species with nonstrict aux so that their order of 
    # declaration won't cause EQPT errors.
    counter <- 0
    max_counter_steps <- 10
    while(counter < max_counter_steps){
      counter <- counter + 1
      previous_state <- aux_entries_with_nonstrict_aux[,"name"]
      for(i in 1:nrow(aux_entries_with_nonstrict_aux)){
        # for each auxiliary species with other aux species in its dissrxn...
        # get s, a vector of containing species in its dissxn
        s <- str_split(aux_entries_with_nonstrict_aux[i, "dissrxn"], " ")[[1]]
        s <- s[2:length(s)]
        s <- s[!s %in% fixed_species]
        s <- s[!s %in% c(aux_entries_with_nonstrict_aux[i, "name"])] # exclude itself
          
        for(sp in s){
          # for each species in the dissrxn...
          for(ii in 1:nrow(aux_entries_with_nonstrict_aux)){
            # for each row in the nonstrict aux df...
              
            if(sp == aux_entries_with_nonstrict_aux[ii,"name"]){
              # if the species is in the nonstrict aux df,
              # copy the original aux species directly after the aux species in
              # its dissrxn.
              if(i < ii){
                  # if the original aux species is in a row before the aux species
                  # in its dissrxn, reorder it.
                  
#                   print("Original aux species:")
#                   print(aux_entries_with_nonstrict_aux[i, "name"])
#                   print("... the aux species in its dissrxn:")
#                   print(aux_entries_with_nonstrict_aux[ii, "name"])
#                   print("The i's")
#                   print(i)
#                   print(ii)

                  aux_entries_with_nonstrict_aux <- insertRow(aux_entries_with_nonstrict_aux, aux_entries_with_nonstrict_aux[i, ], ii+1)

                  # and delete the original row
                  aux_entries_with_nonstrict_aux <- aux_entries_with_nonstrict_aux[-c(i), ]
              }
            }
          }
        }
      }
        
      #print(aux_entries_with_nonstrict_aux[, "name"])
        
      if(all(aux_entries_with_nonstrict_aux[,"name"] == previous_state)){
        break
      }
      previous_state <- aux_entries_with_nonstrict_aux[,"name"]
      
    }
    if(counter >= max_counter_steps){
      m <- paste("Warning: Could not automatically re-order the auxiliary basis species",
                 paste(aux_entries_with_nonstrict_aux[,"name"], collapse=", "), ". This might happen if there are",
                 "many interdependencies of auxiliary basis species in their dissociation reactions.")
      vmessage(m, 1, verbose)
    }
  
  }
    
  thermo_df <- rbind(basis_entries, aux_entries_with_only_strict,
                     aux_entries_with_strict_aux, aux_entries_with_nonstrict_aux,
                     other_entries)

  return(thermo_df)
}


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
    print(m)
  }
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
    
  if(length(sp_name) > 1 & length(unique(sp_name)) == 1){
    sp_name = unique(sp_name)
  }
    
  if(length(sp_name) > 0){
    # get a vector of elements that make up the (non-basis) species
    basis_elem <- (function (x) setdiff(names(unlist(makeup(info(info(x), check.it=F)$formula))), c("H", "O", "Z"))) (thermo_df[, "name"])
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
                   
#   # check that every basis element has a basis pref
#   if(!setequal(basis_elem, names(basis_pref))){
#     missing_basis <- basis_elem[!(basis_elem %in% names(basis_pref))]
#     stop(paste("Error: the element(s)", paste(missing_basis, collapse=", "), "require strict basis species in the database."))
#   }

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
                            
  # check whether basis species in the basis_list are actually given as basis
  # species in thermo_df.
  elements_without_basis <- c()
  for(elem in names(basis_list)){
    if(length(basis_list[[elem]]) == 0){
      elements_without_basis <- c(elements_without_basis, elem)
    }
  }
  if(length(elements_without_basis) > 0){
    msg <- paste("Error during database processing. The following elements",
                  "appear in chemical species but do not have representative",
                  "basis species:", paste(elements_without_basis, collapse=", "),
                  ". Are necessary basis species being excluded?")
    stop(msg)
  }
                            
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
  comp <- makeup(thermo_df[thermo_df["name"]==sp, "formula_modded"][1])
    
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
      
    sp_formula_ox = thermo_df[thermo_df["name"]==sp, "formula_ox_modded"][1]
      
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
          # Species has an element in two different oxidation states
          # e.g., hastingsite has both Fe+2 and Fe+3

          # When there is more than one instance of this element in the non-basis species,
          # find the basis species with the closest oxidation states and assign.
          chosen_basis_species_vec <- c()
          basis_ave_ox_states <- unlist(basis_df["ave_ox_state_of_elem"]) 
        
          names(basis_ave_ox_states) <- unlist(basis_df["formula_ox_modded"])
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
          # Species has an element with one oxidation state
          # e.g., aegerine has Fe+3 (not Fe+3 AND Fe+2)
            
          # Find the basis species with the closest oxidation state to the
          # element inside the non-basis species and assign.
          if(!("Z" %in% names(sp_ox_elem))){
            sp_ox_elem["Z"] <- 0
          }
          basis_ave_ox_states <- unlist(basis_df["ave_ox_state_of_elem"])
          names(basis_ave_ox_states) <-  unlist(basis_df["formula_ox_modded"])
            
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

#   print(sp)
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



                            
suppress_redox_and_generate_dissrxns <- function(thermo_df,
                              db,
                              water_model,
                              template,
                              exceed_Ttr,
                              data0_formula_ox_name,
                              suppress_redox,
                              infer_formula_ox,
                              exclude_category,
                              fixed_species=c("H2O", "H+", "O2(g)", "water", "Cl-", "e-"),
                              verbose=1){
    
  
    
  # specify molecules to balance H, O, and charge (Z)
  HOZ_balancers <- c("H+", "O2(g)", "H2O") # might be dataset-specific (e.g., "O2(g)")
  
#   # load thermodynamic data
#   thermo_df <- read.csv(filename, stringsAsFactors=F)
    
  # remove leading and trailing whitespace from dissrxn and formula_ox columns
  ws_cols <- c("dissrxn", "formula_ox")
  thermo_df[ws_cols] <- lapply(thermo_df[ws_cols], trimws)
    
  # exclude rows based on categories
  if(length(exclude_category) > 0){
    for(cat in names(exclude_category)){
      thermo_df <- thermo_df %>%
        filter(!(!!as.name(cat)) %in% exclude_category[[cat]])
    }
  }
    
  # Ensure that Cl- and other EQ3NR-required species are in thermo_df by marking
  # them as "required" in category 1 of the thermodynamic db CSV.
    
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

        thermo_df[thermo_df[, "name"] == entry, "formula_ox"] <- inferred_formula_ox[[entry]]
      } else if (entry %in% known_species){
        thermo_df[thermo_df[, "name"] == entry, "formula_ox"] <- entry
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
            
          sp_formula_ox <- thermo_df[thermo_df[, "name"] == species_name, "formula_ox"]

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

          thermo_df[thermo_df[, "name"] == species_name, "formula_ox_modded"] <- this_formula_ox_modded
            
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
                                   
      thermo_df[thermo_df[, "name"] == species_name, "formula_modded"] <- formula
    }
  }
                                   
  # reorder thermo_df so that strict basis are first, then aux that depend only
  # on strict basis, then aux that depend on other aux, then everything else.
  # Prevents EQPT errors.
  thermo_df <- order_thermo_df(thermo_df, fixed_species, verbose)
                                   
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
    thermo(OBIGT=thermo()$OBIGT[unique(info(fixed_species)), ]) # replaces the default OBIGT database with user-supplied database
    mod.OBIGT(to_mod_OBIGT, replace=TRUE) # produces a message
  })


  # begin handling basis preferences
  basis_df <- thermo_df %>%
    filter(tag=="basis")
  aux_df <- thermo_df %>%
    filter(tag=="aux")

  basis_pref <- basis_df[, "name"]
                                   
  basis_pref_elements <- lapply(lapply(lapply(basis_df[,"formula_modded"], makeup), names), setdiff, c("Z", "O", "H"))
                                   
  # handle the possibility of a basis species with more than one non-OH element,
  # e.g., glycine as a basis species has C and N.
  for(i in 1:length(basis_pref_elements)){
      
    if(length(basis_pref_elements[[i]]) > 1){
        
        msg <- paste0("Error during database processing. The strict basis ",
                     "species '", basis_pref[i], "' has more than one ",
                     "element besides O and H. This can be solved by making ",
                     "this species an auxiliary basis species with a dissociation ",
                     "reaction into strict basis species representing the elements ",
                     paste(basis_pref_elements[[i]], collapse=" "))
        stop(msg)

    }
  }
                                   
  names(basis_pref) <- basis_pref_elements
                                   
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
  if(!identical(aux_pref, character(0))){
  keep <- c()
  for(i in 1:length(aux_pref)){
    if(makeup(thermo_df[thermo_df[, "name"]==aux_pref[[i]], "formula_modded"])[aux_pref_names[[i]]] == 1){
      keep <- c(keep, i)
    }
  }
  aux_pref <- aux_pref[keep]
  aux_pref_names <- aux_pref_names[keep]
    
  names(aux_pref) <- aux_pref_names
  }else{
    aux_pref <- list()
    aux_pref_names <- c()
  }
                           
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
        
        thermo_df[thermo_df[, "name"] == valid_aux, "tag"] <- "basis"

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
  thermo_df[, "regenerate_dissrxn"] <- FALSE
  all_basis_species <- c()
  for(species in thermo_df[, "name"]){
    
    tag <- thermo_df[thermo_df[, "name"] == species, "tag"]
    if(tag != "basis"){
      dissrxn <- thermo_df[thermo_df[, "name"] == species, "dissrxn"]
      dissrxn <- strsplit(dissrxn, " ")[[1]] # split the rxn into coeffs and species
      dissrxn_names <- dissrxn[c(FALSE, TRUE)] # get names of reactants and products
      all_basis_species <- c(all_basis_species, dissrxn_names[2:length(dissrxn_names)])
      dissrxn_ispecies <- suppressMessages(info(dissrxn_names))
      dissrxn_coefs <- dissrxn[c(TRUE, FALSE)] # get coeffs of reactants and products
      dissrxn_coefs <- as.numeric(dissrxn_coefs) # convert coeffs from str to numeric
      tryCatch({
        subcrt_bal(dissrxn_ispecies, dissrxn_coefs)
      }, error=function(e){
        thermo_df[thermo_df[, "name"] == species, "regenerate_dissrxn"] <<- TRUE # assign global variable with <<- because this is within the error function
      })
    }
  }
  all_basis_species <- unique(all_basis_species)
                                   
  # check whether there are species in dissociation reactions that are missing or
  # not tagged in the database as 'basis' or 'aux'.
  missing_basis_species <- c()
  species_to_blame <- c()
  species_in_dissrxns_that_are_not_basis <- c()
  for(basis in all_basis_species){
    if(!(basis %in% thermo_df$name) && !(basis %in% fixed_species)){
       missing_basis_species <- c(missing_basis_species, basis)
    }else if(!(basis %in% basis_df$name) && !(basis %in% aux_df$name) && !(basis %in% fixed_species)){
      species_in_dissrxns_that_are_not_basis <- c(species_in_dissrxns_that_are_not_basis, basis)
    }
  }
                                   
  if(length(missing_basis_species) == 1){
    if(is.na(missing_basis_species)){
      missing_basis_species = c()
    }
  }
                                   
  if(length(missing_basis_species) > 0){
    for(species in thermo_df[, "name"]){
      tag <- thermo_df[thermo_df[, "name"] == species, "tag"]
      if(tag != "basis"){
        dissrxn <- thermo_df[thermo_df[, "name"] == species, "dissrxn"]
        dissrxn <- strsplit(dissrxn, " ")[[1]] # split the rxn into coeffs and species
        dissrxn_names <- dissrxn[c(FALSE, TRUE)] # get names of reactants and products
        for(n in dissrxn_names){
          if(n %in% missing_basis_species){
            species_to_blame <- c(species_to_blame, species)
          }
        }
      }
    }
    species_to_blame <- unique(species_to_blame)
      
    msg_missing_basis <- paste("The following species appear in dissociation",
      "reactions but are not included in the database: [", paste(missing_basis_species, collapse=", "),
      "]. These missing species are found in the dissociation reactions of [",
      paste(species_to_blame, collapse=", "), "].")
  }else{
    msg_missing_basis <- ""
  }
  
  species_to_blame_have_bad_basis <- c()
  if(length(species_in_dissrxns_that_are_not_basis) > 0){
    for(species in thermo_df[, "name"]){
      tag <- thermo_df[thermo_df[, "name"] == species, "tag"]
      if(tag != "basis"){
        dissrxn <- thermo_df[thermo_df[, "name"] == species, "dissrxn"]
        dissrxn <- strsplit(dissrxn, " ")[[1]] # split the rxn into coeffs and species
        dissrxn_names <- dissrxn[c(FALSE, TRUE)] # get names of reactants and products
        dissrxn_names <- dissrxn_names[2:length(dissrxn_names)]
        for(n in dissrxn_names){
          if(n %in% species_in_dissrxns_that_are_not_basis){
            species_to_blame_have_bad_basis <- c(species_to_blame_have_bad_basis, species)
          }
        }
      }
    }
    species_to_blame_have_bad_basis <- unique(species_to_blame_have_bad_basis)
    msg_not_basis <- paste("The following species appear in dissociation",
      "reactions but are not tagged as strict or auxiliary basis species: [",
      paste(species_in_dissrxns_that_are_not_basis, collapse=", "),
      "]. These non-basis species are found in the dissociation reactions of [",
      paste(species_to_blame_have_bad_basis, collapse=", "), "].")
  }else{
    msg_not_basis <- ""
  }
  if(length(missing_basis_species) > 0 | length(species_in_dissrxns_that_are_not_basis) > 0){
    stop(paste("One or more errors were encountered during database processing.",
               msg_missing_basis, msg_not_basis))
  }
                                   
  df_needs_dissrxns <- thermo_df %>%
    filter(tag != "basis") %>%
    filter(regenerate_dissrxn == TRUE)
      
  if(nrow(df_needs_dissrxns) > 0){
    if(length(suppress_redox) == 0){
      needs_dissrxns_message <- paste("Balanced dissociation reactions are missing for species:",
                                      paste(unique(unlist(df_needs_dissrxns["name"])), collapse=", "))
    }else{
      needs_dissrxns_message <- paste("Balanced and/or redox-suppressed dissociation reactions are missing for species:",
                                      paste(unique(unlist(df_needs_dissrxns["name"])), collapse=", "))
    }
      
    vmessage(needs_dissrxns_message, 1, verbose)
    vmessage("Generating dissociation reactions for these species using strict and auxiliary basis species containing a maximum of one atom of one element besides O and H...", 1, verbose)
  }
                    
  # generate dissociation reactions
  dissrxns <- get_dissrxn(sp_name=unlist(df_needs_dissrxns["name"]),
                                           basis_pref=basis_pref,
                                           aux_pref=aux_pref,
                                           HOZ_balancers=HOZ_balancers,
                                           thermo_df=thermo_df,
                                           verbose=verbose,
                                           redox_elem_states=redox_elem_states)
                                   
  # Produce a warning message about which dissrxns were (re)generated and what they are.
  if(nrow(df_needs_dissrxns) > 0){
    names <- df_needs_dissrxns[["name"]]
    generated_dissrxns <- dissrxns[names]
    nonbasis_idx <- unlist(lapply(lapply(lapply(generated_dissrxns, strsplit, " "), `[[`, 1), FUN=function(x) length(x)!=4))
    basis_idx <- !nonbasis_idx
    nonbasis_names <- names(nonbasis_idx)[nonbasis_idx]
    nonbasis_names <- unique(nonbasis_names)
    basis_names <- names(basis_idx)[basis_idx]
    vmessage(paste(nonbasis_names, ":", generated_dissrxns[nonbasis_names], "\n"), 1, verbose)
    if(!identical(basis_names, character(0))){
      vmessage(paste("Species that have been converted into strict basis:", paste(unique(basis_names), collapse=", ")), 1, verbose)
    }

    # replace dissrxns with any regenerated dissrxns
    for(name in names){
        
#       print(name)
#       print(dissrxns[[name]])
#       print(thermo_df[thermo_df["name"]==name, "dissrxn"])
        
      thermo_df[thermo_df["name"]==name, "dissrxn"] <- dissrxns[[name]]
    }
  }
#   print(tail(thermo_df))
                                  
  thermo_df[is.na(thermo_df)]=''
  
  out_list = list("OBIGT_df"=thermo_df, "dissrxns"=dissrxns, "basis_pref"=basis_pref)
                                  
  return(out_list)
}