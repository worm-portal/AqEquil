# hidden functions in Jeff Dick's CHNOSZ package

info.character <- function(species, state=NULL, check.protein=TRUE) {
  # returns the rownumbers of thermo()$OBIGT having an exact match of 'species' to
  # thermo()$OBIGT$[species|abbrv|formula] or NA otherwise
  # a match to thermo()$OBIGT$state is also required if 'state' is not NULL

  # (first occurence of a match to species is returned otherwise)
  thermo <- get("thermo", CHNOSZ)
  # find matches for species name, abbreviation or formula
  matches.species <- thermo$OBIGT$name==species | thermo$OBIGT$abbrv==species | thermo$OBIGT$formula==species
  # since thermo()$OBIGT$abbrv contains NAs, convert NA results to FALSE
  matches.species[is.na(matches.species)] <- FALSE
  # turn it in to no match if it's a protein in the wrong state
  ip <- pinfo(species)
  if(any(matches.species) & !is.na(ip) & !is.null(state)) {
    matches.state <- matches.species & grepl(state, thermo$OBIGT$state)
    if(!any(matches.state)) matches.species <- FALSE
  }
  # no match, not available
  if(!any(matches.species)) {
    # unless it's a protein
    if(check.protein) {
      # did we find a protein? add its properties to OBIGT
      if(!is.na(ip)) {
        # here we use a default state from thermo()$opt$state
        if(is.null(state)) state <- thermo$opt$state
        # add up protein properties
        eos <- protein.OBIGT(ip, state=state)
        # the real assignment work 
        nrows <- suppressMessages(mod.OBIGT(eos))
        thermo <- get("thermo", CHNOSZ)
        matches.species <- rep(FALSE, nrows)
        matches.species[nrows] <- TRUE
      } else return(NA)
    } else return(NA)
  }
  # do we demand a particular state
  if(!is.null(state)) {
    # special treatment for H2O: aq retrieves the liq
    if(species %in% c("H2O", "water") & state=="aq") state <- "liq"
    # the matches for both species and state
    matches.state <- matches.species & state == thermo$OBIGT$state
    if(!any(matches.state)) {
      # the requested state is not available for this species
      available.states <- thermo$OBIGT$state[matches.species]
      if(length(available.states)==1) a.s.verb <- "is" else a.s.verb <- "are"
      a.s.text <- paste("'", available.states, "'", sep="", collapse=" ")
      message("info.character: requested state '", state, "' for ", species, 
        " but only ", a.s.text, " ", a.s.verb, " available")
      # warn about looking for aqueous methane (changed to CH4) 20200707
      if(identical(species, "methane") & identical(state, "aq")) {
        warning("'methane' is not an aqueous species; use 'CH4' instead\nTo revert to the old behavior, run mod.OBIGT(info('CH4'), name = 'methane')")
      }
      return(NA)
    }
    matches.species <- matches.state
  }
  # all of the species that match
  ispecies.out <- ispecies <- which(matches.species)
  # processing for more than one match
  if(length(ispecies) > 1) {
    # if a single name matches, use that one (useful for distinguishing pseudo-H4SiO4 and H4SiO4) 20171020
    matches.name <- matches.species & thermo$OBIGT$name==species
    if(sum(matches.name)==1) ispecies.out <- which(matches.name)
    else ispecies.out <- ispecies[1]  # otherwise, return only the first species that matches
    # let user know if there is more than one state for this species
    mystate <- thermo$OBIGT$state[ispecies.out]
    ispecies.other <- ispecies[!ispecies %in% ispecies.out]
    otherstates <- thermo$OBIGT$state[ispecies.other]
    # for minerals (cr), use the word "phase"; otherwise, use "state" 20190209
    word <- "state"
    # substitute the mineral name for "cr" 20190121
    if(mystate == "cr" | sum(otherstates=="cr") > 1) {
      word <- "phase"
      otherstates[otherstates=="cr"] <- thermo$OBIGT$name[ispecies.other[otherstates=="cr"]]
    }
    transtext <- othertext <- ""
    # we count, but don't show the states for phase transitions (cr2, cr3, etc)
    istrans <- otherstates %in% c("cr2", "cr3", "cr4", "cr5", "cr6", "cr7", "cr8", "cr9")
    if(mystate=="cr") {
      # if we are "cr" we show the number of phase transitions
      ntrans <- sum(istrans)
      if(ntrans == 1) transtext <- paste(" with", ntrans, "phase transition")
      else if(ntrans > 1) transtext <- paste(" with", ntrans, "phase transitions")
      # if it's not already in the species name, substitute the mineral name for "cr" 20190121
      if(species != thermo$OBIGT$name[ispecies.out]) mystate <- thermo$OBIGT$name[ispecies.out]
    }
    otherstates <- otherstates[!istrans]
    if(length(otherstates) == 1) othertext <- paste0("; other available ", word, " is ", otherstates)
    if(length(otherstates) > 1) othertext <- paste0("; other available ", word, "s are ", paste(otherstates, collapse=", "))
    if(transtext != "" | othertext != "") {
      starttext <- paste0("info.character: found ", species, "(", mystate, ")")
      message(starttext, transtext, othertext)
    }
  }
  return(ispecies.out)
}

can.be.numeric <- function(x) {
  # return FALSE if length of argument is zero
  if(length(x) == 0) FALSE else
  if(length(x) > 1) as.logical(sapply(x, can.be.numeric)) else {
    if(is.numeric(x)) TRUE else
    if(!is.na(as.numeric.nowarn(x))) TRUE else
    if(x %in% c('.','+','-')) TRUE else FALSE
  }
}
        
# something like R's as.numeric(), but without the "NAs introduced by coercion" warnings
# (needed because testthat somehow detects the warnings suppressed by suppressWarnings) 20170427
as.numeric.nowarn <- function(x) {
  if(length(x) == 0) numeric() else
  if(length(x) > 1) sapply(x, as.numeric.nowarn) else
  # http://stackoverflow.com/questions/12643009/regular-expression-for-floating-point-numbers
  if(grepl("^[+-]?([0-9]*[.])?[0-9]+$", x)) as.numeric(x) else NA_real_
}

outvert <- function(value,units) {
  # converts the given value from the given units to
  # those specified in thermo()$opt
  units <- tolower(units)
  opt <- get("thermo", CHNOSZ)$opt
  if(units %in% c('c','k')) {
    if(units=='c' & opt$T.units=='K') return(convert(value,'k'))
    if(units=='k' & opt$T.units=='C') return(convert(value,'c'))
  }
  if(units %in% c('j','cal')) {
    if(units=='j' & opt$E.units=='Cal') return(convert(value,'cal'))
    if(units=='cal' & opt$E.units=='J') return(convert(value,'j'))
  }
  if(units %in% c('bar','mpa')) {
    if(units=='mpa' & opt$P.units=='bar') return(convert(value,'bar'))
    if(units=='bar' & opt$P.units=='MPa') return(convert(value,'mpa'))
  }
  return(value)
}

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

subcrt_bal <- function(species, coeff, state, P){
    
  thermo <- get("thermo", CHNOSZ)
  
    
#   if(is.numeric(species[1])) {
    ispecies <- species
    species <- as.character(thermo$OBIGT$name[ispecies])
    state <- as.character(thermo$OBIGT$state[ispecies])
#   } else {
#     # from names, get species indices and states and possibly
#     # keep track of phase species (cr,cr2 ...)
#     ispecies <- numeric()
#     newstate <- character()
#     for(i in 1:length(species)) {
#       # get the species index for a named species
#       if(!can.be.numeric(species[i])) si <- info.character(species[i], state[i])
#       else {
#         # check that a numeric argument is a rownumber of thermo()$OBIGT
#         si <- as.numeric(species[i])
#         if(!si %in% 1:nrow(thermo$OBIGT)) stop(paste(species[i], "is not a row number of thermo()$OBIGT"))
#       }
#       # that could have the side-effect of adding a protein; re-read thermo
#       thermo <- get("thermo", CHNOSZ)
#       if(is.na(si[1])) stop('no info found for ',species[i],' ',state[i])
#       if(!is.null(state[i])) is.cr <- state[i]=='cr' else is.cr <- FALSE
#       if(thermo$OBIGT$state[si[1]]=='cr' & (is.null(state[i]) | is.cr)) {
#         newstate <- c(newstate,'cr')
#         ispecies <- c(ispecies,si[1])
#       } else {
#         newstate <- c(newstate,as.character(thermo$OBIGT$state[si[1]]))
#         ispecies <- c(ispecies,si[1])
#       }
#     }
#   }
          
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
      #message("subcrt: reaction is not balanced; it is missing this composition:")
      # we have to do this awkward dance to send a formatted matrix to message
      #message(paste(capture.output(print(miss)), collapse="\n"))
      # look for basis species that have our compositoin
      tb <- thermo$basis
      if(!is.null(tb)) {
        if(all(names(miss) %in% colnames(tb)[1:nrow(tb)])) {
          # the missing composition as formula
          ft <- as.chemical.formula(miss)
          # the basis species needed to supply it
          bc <- species.basis(ft)
          # drop zeroes
          bc.new <- bc[,(bc[1,]!=0),drop=FALSE]
          # and get the states
          b.state <- as.character(thermo$basis$state)[bc[1,]!=0]
          bc <- bc.new
          # special thing for Psat
          if(identical(P[[1]], "Psat")) P <- "Psat"
          else P <- outvert(P,"bar")
          # add to logact values if present
          if(!is.null(logact)) {
            ila <- match(colnames(bc),rownames(thermo$basis))
            nla <- !(can.be.numeric(thermo$basis$logact[ila]))
            if(any(nla)) warning('subcrt: logact values of basis species',
              c2s(rownames(thermo$basis)[ila]),'are NA.')
            logact <- c(logact,thermo$basis$logact[ila])
          }
          # warn user and do it!
          ispecies.new <- tb$ispecies[match(colnames(bc),rownames(tb))]
          b.species <- thermo$OBIGT$formula[ispecies.new]
          if(identical(species,b.species) & identical(state,b.state))
            message("subcrt: balanced reaction, but it is a non-reaction; restarting...")
          #else message('subcrt: adding missing composition from basis definition and restarting...')
          newspecies <- c(species, tb$ispecies[match(colnames(bc), rownames(tb))])
          newcoeff <- c(coeff, as.numeric(bc[1, ]))
          newstate <- c(state, b.state)
          return(list('newspecies'=newspecies, 'newcoeff'=newcoeff, 'newstate'=newstate))
        } else warnings <- c(warnings, paste('reaction among', paste(species, collapse = ","), 'was unbalanced, missing', as.chemical.formula(miss)))
      } else warnings <- c(warnings, paste('reaction among', paste(species, collapse = ","), 'was unbalanced, missing', as.chemical.formula(miss)))
    }
}


# basis("CHNOS+")

# start_time <- Sys.time()
# subcrt_bal(species=c(18, 16), coeff=c(-1, 1), state=c("aq", "aq"))
# print(Sys.time()-start_time)
                    
# start_time <- Sys.time()
# subcrt(species=c(18, 16), coeff=c(-1, 1), state=c("aq", "aq"))$reaction
# print(Sys.time()-start_time)