# preprocess_for_EQ3 -----------------------------------------------------------
# Summary: import, clean, and generate a batch of EQ3 input files from a
#   user-supplied csv table of samples.
# Author: Grayson Boyer (gmboyer@asu.edu)
# Created for the ENKI project
# Last updated Jan 31, 2021


suppressMessages({
  library(CHNOSZ)
  library(dplyr)
})


vprint <- function(string, verbose=2){
  # print messages depending on desired level of 'verbose'-ness.
  if(verbose == 2){
    message(string)
  } else if(verbose == 1){
    if(grepl("warning", tolower(string)) | grepl("error", tolower(string))){
      message(string)
    }
  }
}


header_unit <- function(header, keepcase=FALSE){
  # get units from a header
  # (unit should always be at end of header after an underscore)
  split_header <- strsplit(header, "_")[[1]] # split header by underscore
  if(!keepcase){
    tolower(split_header[length(split_header)]) # get last string in split
  } else {
    split_header[length(split_header)] # get last string in split
  }
}

uc_molal <- function(value=1, chemical="H+", unit="ppm"){
  # convert unit name to lowercase
  unit <- tolower(unit)

  # convert ppb to ppm
  if(unit == "ppb"){
    value <- value/1000
  }

  # convert ppm or mg/L to molality
  if(unit == "ppb" | unit == "ppm" | unit == "mg/l"){
    value <- value*0.001/mass(info(info("water"))$formula)
    return(value)
  } else if(unit == "molality" | unit == "molal"){
    return(value)
  } else {
    stop("Error in uc_molal(): unit must be either ppb, ppm, or mg/L")
  }
}


header_species <- function(header){
    # get chemical species from header
    # (should always be at start of header before an underscore)
  split_header <- strsplit(header, "_")[[1]][1] # split header by underscore
}

preprocess <- function(input_filename,
                       exclude,
                       redox_flag,
                       default_logfO2,
                       charge_balance_on,
                       suppress_missing,
                       suppress,
                       verbose=2){
    # Start the timer:
    ptm <- proc.time()
    
    # function to convert ppb, ppm, or mg/L to molality
    acceptable_units <- c("ppb", "ppm", "mg/l", "molality", "molal")
    EQ3_jflags <- c("Suppressed", "Molality", "Molarity", "mg/L",
                    "mg/kg.sol", "Alk., eq/kg.H2O", "Alk., eq/L",
                    "Alk., eq/kg.sol", "Alk., mg/L CaCO3", "Alk., mg/L HCO3-",
                    "Log activity", "Log act combo", "Log mean act", "pX",
                    "pH", "pHCl", "pmH", "pmX", "Hetero. equil.",
                    "Homo. equil.", "Make non-basis")

    # Read input -------------------------------------------------------------------
    vprint("Reading input file...", verbose=verbose)

    r <- read.csv(input_filename,
                  header=T,
                  check.names=F,
                  stringsAsFactors=F)

    df <- as.data.frame(r)                                     

    # delete blank columns if they exist
    # df <- df[, names(df) != ""]

    num_cols <- names(df[!(names(df) %in% exclude)])
    num_cols <- match(num_cols, colnames(df))
    num_cols <- num_cols[num_cols != 1]


    # Fix headers ------------------------------------------------------------------
    vprint("Handling headers...", verbose=verbose)

    # specify headers and subheaders
    header1 <- names(df)
    header2 <- df[1, ]

    # create concatenated header_subheader column names
    for(i in 1:length(header1)){
      if(is.na(header2[i])){
        header2[i] <- ""
      }
      if(header2[i] != ""){
        colnames(df)[i] <- paste(header1[i], header2[i], sep="_")
      }else{
        colnames(df)[i] <- header1[i]
      }
        
    }
    #colnames(df) <- paste(header1, header2, sep="_")
    df = df[-1, ] # remove subheader row
    colnames(df)[1] <- "Sample" # rename "_" as sample name header

    # Convert data to numeric ------------------------------------------------------
    vprint("Ensuring data is numeric...", verbose=verbose)

    # convert columns to numeric
    df[, num_cols] <- mutate_all(df[, num_cols], function(x) as.numeric(as.character(x)))

    # convert these special columns to numeric even if excluded from aqueous block
    special_numeric <- c("rho_", "Pressure_", "redox_", "Temperature_", "pe_", "Eh_")               

    for(col in special_numeric){
      if(TRUE %in% grepl(col, colnames(df))){
        this_match <- grep(col, colnames(df))
        df[, this_match] <- as.numeric(as.character(df[, this_match]))
      }
    }              

    # Handle exclusions -----------------------------------------------------------
    # add columns to exclusion list
    exclude <- c(exclude, "", "Sample", "rho", "Pressure", "redox", "Temperature", "pe", "Eh")                             

    # Calculate density of water (rho) ---------------------------------------------
    vprint("Calculating density of water (rho)...", verbose=verbose)

    # find temperature (in degK)
    if("Temperature_degC" %in% names(df)){
      temp_degC <- df[, "Temperature_degC"]
      temp_degK <- temp_degC + 273.15
    } else if("Temperature_degK" %in% names(df)){
      temp_degK <- df[, "Temperature_degK"]
      temp_degC <- temp_degK - 273.15
    } else {
      stop("Error: need a column for 'Temperature' with units in 'degC' or 'degK'")
    }

    # find pressure (in bar)
    if("Pressure_bar" %in% names(df)){
      pressure_bar <- df[, "Pressure_bar"]
    } else {
      vprint(paste("No column called 'Pressure' with unit 'bar' found.",
                    "Using a pressure of 1 bar for each sample."), verbose=verbose)
      psat_press <- water("Psat", T=temp_degK)[[1]]
      pressure_bar <- psat_press
    }

    # calculate rho and append
    df <- mutate(df, rho=water("rho", T=temp_degK, P=pressure_bar)[[1]])
    df$rho <- df$rho/1000 #convert rho to g/cm3


    # Handle redox block -----------------------------------------------------------
    vprint("Handling redox block options...", verbose=verbose)

    df$redox_flag <- redox_flag
    df$redox_value <- NA
    df$redox_unit <- NA

    # set to TRUE if a warning is given
    # about a missing column to prevent
    # repeated warnings for each sample
    warned_about_redox_column <- FALSE

    for (row in 1:nrow(df)){
      this_redox_flag <- redox_flag
      assigned <- FALSE

      # use gaseous O2 in aqueous species block
      if(redox_flag == -3){
        #redox_col_index <- grep("O2,g", names(df))
        redox_col_index <- grep("O2,g", names(df), fixed = T, value = T)
        redox_col_index <- redox_col_index[substring(redox_col_index, 1, nchar("O2,g")) == "O2,g"]
        if(length(redox_col_index) > 0){
          if(!is.na(df[row, redox_col_index])){
            this_redox_value <- sprintf("%.4E", df[row, redox_col_index])
            this_redox_unit <- "using O2,g in aqueous species block"
          } else {
            vprint(paste0("Warning: non-numeric 'O2,g' value in sample ",
                           df[row, "Sample"], ". Resorting to using ",
                           "Log fO2 (log bars) with a value of ", default_logfO2),
                           verbose=verbose)
            this_redox_flag <- 0
            this_redox_value <- sprintf("%.4E", default_logfO2)
            this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
            assigned <- TRUE
          }
        } else {
          if(!warned_about_redox_column){
            vprint(paste("Warning: no 'O2,g' column found. Resorting to using",
                         "Log fO2 (log bars) with a value of", default_logfO2),
                         verbose=verbose)
            warned_about_redox_column <- TRUE
          }
          this_redox_flag <- 0
          this_redox_value <- sprintf("%.4E", default_logfO2)
          this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
          assigned <- TRUE
        }
      }

      # specify pe in pe units
      if(redox_flag == -2){
        redox_col_index <- grep("^pe_", names(df))
        if(length(redox_col_index) > 0){
          if(!is.na(df[row, redox_col_index])){
            this_redox_value <- sprintf("%.4E", df[row, redox_col_index])
            this_redox_unit <- "pe, pe units"
          } else {
            vprint(paste0("Warning: non-numeric pe value in sample ",
                           df[row, "Sample"], ". Resorting to using ",
                           "Log fO2 (log bars) with a value of ", default_logfO2),
                   verbose=verbose)
            this_redox_flag <- 0
            this_redox_value <- sprintf("%.4E", default_logfO2)
            this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
            assigned <- TRUE
          }
        } else {
          if(!warned_about_redox_column){
            vprint(paste0("Warning: no 'pe' column found. Resorting to using",
                          "Log fO2 (log bars) with a value of ", default_logfO2),
                   verbose=verbose)
            warned_about_redox_column <- TRUE
          }
          this_redox_flag <- 0
          this_redox_value <- sprintf("%.4E", default_logfO2)
          this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
          assigned <- TRUE
        }
      }

      # specify Eh in volts
      if(redox_flag == -1){
        redox_col_index <- grep("Eh_volts", names(df))
        if(length(redox_col_index) > 0){
          if(!is.na(df[row, redox_col_index])){
            this_redox_value <- sprintf("%.4E", df[row, redox_col_index])
            this_redox_unit <- "Eh, volts"
          } else {
            vprint(paste0("Warning: non-numeric Eh value in sample ",
                           df[row, "Sample"], ". Resorting to using ",
                           "Log fO2 (log bars) with a value of ", default_logfO2),
                   verbose=verbose)
            this_redox_flag <- 0
            this_redox_value <- sprintf("%.4E", default_logfO2)
            this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
            assigned <- TRUE
          }
        } else {
          if(!warned_about_redox_column){
            vprint(paste0("Warning: no 'Eh' column found with 'volts' units. ",
                          "Resorting to using Log fO2 (log bars) with a value of ",
                          default_logfO2),
                   verbose=verbose)
            warned_about_redox_column <- TRUE
          }
          this_redox_flag <- 0
          this_redox_value <- sprintf("%.4E", default_logfO2)
          this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
          assigned <- TRUE
        }
      }

      if(redox_flag == 0 && assigned == FALSE){
        redox_col_index <- grep("logfO2", names(df))
        if(length(redox_col_index) > 0){
          if(!is.na(df[row, redox_col_index])){
            this_redox_value <- sprintf("%.4E", df[row, redox_col_index])
            this_redox_unit <- "Log fO2 (log bars) [from logfO2 column]"
          } else {
            vprint(paste0("Warning: non-numeric logfO2 value in sample ",
                           df[row, "Sample"], ". Resorting to using ",
                           "a logfO2 value of ", default_logfO2),
                   verbose=verbose)
            this_redox_flag <- 0
            this_redox_value <- sprintf("%.4E", default_logfO2)
            this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
            assigned <- TRUE
          }
        } else {
          if(!warned_about_redox_column){
            vprint(paste("Warning: no 'logfO2' column found. Attempting to find a",
                        "column for aqueous O2 to estimate logfO2 at sample temperature and",
                        "pressure..."),
                   verbose=verbose)
            warned_about_redox_column <- TRUE
          }
          redox_col_index <- grep("(^O2,AQ)|(^O2_)", names(df)) # TODO: make more flexible (June 25, 2020)
          if(length(redox_col_index) > 0){
            if(!is.na(df[row, redox_col_index])){
              header_name <- names(df)[redox_col_index]
              this_header_unit <- header_unit(header_name)
              if(this_header_unit %in% acceptable_units){
                O2_molal <- uc_molal(value=df[row, redox_col_index],
                                     chemical="O2", unit=this_header_unit)
                if(is.na(temp_degC[row]) | is.na(pressure_bar[row])){
                  vprint(paste0("Warning: non-numeric temperature or pressure value ",
                                 "in sample ", df[row, "Sample"], ". Resorting to ",
                                 "using Log fO2 (log bars) with a value of ", default_logfO2),
                         verbose=verbose)
                  this_redox_flag <- 0
                  this_redox_value <- sprintf("%.4E", default_logfO2)
                  this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
                  assigned <- TRUE
                } else {
                  suppressMessages({
                  logfO2 <- log10(O2_molal*10^subcrt(c("O2", "O2"),
                                                     c(-1, 1),
                                                     c("aq", "g"),
                                                     T=temp_degC[row],
                                                     P=pressure_bar[row])$out$logK)
                  })
                  this_redox_value <- sprintf("%.4E", logfO2)
                  this_redox_unit <- paste("Log fO2 (log bars) [calculated from O2(aq) = O2(g) at",
                                           "temperature and pressure of sample]")
                  assigned <- TRUE
                }
              } else {
                vprint(paste("Warning: column found for aqueous O2, but units are not ppb, ppm",
                               "mg/L, or molality. Resorting to using Log fO2",
                               "(log bars) with a value of", default_logfO2, "for all samples."),
                       verbose=verbose)
                this_redox_flag <- 0
                this_redox_value <- sprintf("%.4E", default_logfO2)
                this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
                assigned <- TRUE
              }
            } else {
              vprint(paste0("Warning: non-numeric aqueous O2 value in sample ",
                             df[row, "Sample"], ". Resorting to using ",
                             "Log fO2 (log bars) with a value of ", default_logfO2),
                             verbose=verbose)
              this_redox_flag <- 0
              this_redox_value <- sprintf("%.4E", default_logfO2)
              this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
              assigned <- TRUE
            }
          } else {
              if(!warned_about_redox_column){
                vprint(paste("Warning: a column for aqueous O2 was not found. Resorting to",
                           "using Log fO2 (log bars) with a value of ",
                           default_logfO2),
                       verbose=verbose)
                warned_about_redox_column <- TRUE
              }
              this_redox_flag <- 0
              this_redox_value <- sprintf("%.4E", default_logfO2)
              this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
              assigned <- TRUE
          }
        }
      }

      # specify aux species redox couple
      if(redox_flag == 1){
        redox_col_index <- grep(redox_aux, names(df), fixed = T, value = T)
        redox_col_index <- redox_col_index[substring(redox_col_index, 1, nchar(redox_aux)) == redox_aux]
        if(length(redox_col_index) > 0){
          if(!is.na(df[row, redox_col_index])){
            this_redox_value <- sprintf("%.4E", df[row, redox_col_index])
            this_redox_unit <- paste(redox_aux, "aux. sp.")
          } else {
            vprint(paste0("Warning: non-numeric ", redox_aux, " value in sample ",
                           df[row, "Sample"], ". Resorting to using ",
                           "Log fO2 (log bars) with a value of ", default_logfO2),
                   verbose=verbose)
            this_redox_flag <- 0
            this_redox_value <- sprintf("%.4E", default_logfO2)
            this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
            assigned <- TRUE
          }
        } else {
              if(!warned_about_redox_column){
                vprint(paste("Warning: no", redox_aux, "column found. Resorting to using",
                             "Log fO2 (log bars) with a value of ", default_logfO2),
                       verbose=verbose)
                warned_about_redox_column <- TRUE
              }
          this_redox_flag <- 0
          this_redox_value <- sprintf("%.4E", default_logfO2)
          this_redox_unit <- paste0("Log fO2 (log bars) [default = ", default_logfO2, "]")
          assigned <- TRUE
        }
      }

      # append redox values
      df[row, "redox_flag"] <- this_redox_flag
      df[row, "redox_value"] <- this_redox_value
      df[row, "redox_unit"] <- this_redox_unit

    } # end loop for redox block

    # Handle charge balance block --------------------------------------------------
    vprint("Handling charge balance options...", verbose=verbose)

    if(charge_balance_on == "none"){
      eq3.cb_block <- paste("\n|Electrical balancing option (iebal3):                                         |",
    "|  [x] ( 0) No balancing is done                                               |",
    "|  [ ] ( 1) Balance on species |None                    | (uebal)              |", sep="\n")

    } else {
      eq3.cb_block <- paste(
      "\n|Electrical balancing option (iebal3):                                         |",
      "|  [ ] ( 0) No balancing is done                                               |",
      "|  [x] ( 1) Balance on species |", format(charge_balance_on, width=24), "| (uebal)              |",
      sep="\n")
    }

    # Handle input file directory --------------------------------------------------
    vprint("Deleting previous rxn_3i folder...", verbose=verbose)

    # get user's working directory
    wd <- getwd()

    # delete the rxn_3i directory if it exists
    unlink("rxn_3i", recursive = TRUE)

    # create a fresh rxn_3i directory
    vprint("Creating fresh rxn_3i folder...", verbose=verbose)
    dir.create("rxn_3i", showWarnings = FALSE)

    # get path name of rxn_3i folder
    input_path <- paste(wd, "/rxn_3i", sep="")

    # set working directory to rxn_3i folder
    setwd(input_path)

    # create unique df row names that will become .3i filenames
    row_order <- make.names(df$Sample, unique=T)
    rownames(df) <- row_order


    # Write .3i input file ---------------------------------------------------------
    vprint("Creating 3i files...", verbose=verbose)
    for (row in 1:nrow(df)){
      vprint(paste0("Working on file ", row, " out of ", nrow(df)), verbose=verbose)

      eq3.filename <- paste0(c(rownames(df)[row], ".3i"), collapse = "") #create the name of the file

      eq3.header1 <- paste("|------------------------------------------------------------------------------|",
    "| Title                  | (utitl(n))                                          |",
    "|------------------------------------------------------------------------------|",
    "|                                                                              |",
    "|", sep="\n")

      eq3.samplename <- paste(c("Sample:", format(df[row, "Sample"], width=70)), collapse = " ")

      eq3.header2 <- paste("|",
    "|                                                                              |",
    "|------------------------------------------------------------------------------|",
    "|Special Basis Switches (for model definition only)       | (nsbswt)           |",
    "|------------------------------------------------------------------------------|",
    "|Replace |None                                            | (usbsw(1,n))       |",
    "|   with |None                                            | (usbsw(2,n))       |",
    "|------------------------------------------------------------------------------|",
    "|Temperature (C)         | ", sep="\n")

      eq3.temperature <- sprintf("%.5E", temp_degC[row])

      eq3.header3 <-              paste("| (tempc)                                |",
    "|------------------------------------------------------------------------------|",
    "|Pressure option (jpres3):                                                     |",
    "|  [ ] ( 0) Data file reference curve value                                    |",
    "|  [ ] ( 1) 1.013-bar/steam-saturation curve value                             |",
    paste0("|  [x] ( 2) Value (bars) | ", sprintf("%.5E", pressure_bar[row]), "| (press)                                |"),
    "|------------------------------------------------------------------------------|",
    "|Density (g/cm3)         | ", sep="\n")

      eq3.density <- sprintf("%.5E", df[row, "rho"])

      eq3.header4 <-               paste("| (rho)                                  |",
    "|------------------------------------------------------------------------------|",
    "|Total dissolved solutes option (itdsf3):                                      |",
    "|  [x] ( 0) Value (mg/kg.sol) | 0.00000E+00| (tdspkg)                          |",
    "|  [ ] ( 1) Value (mg/L)      | 0.00000E+00| (tdspl)                           |",
    "|------------------------------------------------------------------------------|",
                                         sep="\n")

                     # cb_block will be pasted here (one for all samples)

      eq3.header5 <- paste("\n|------------------------------------------------------------------------------|",
        "|Default redox constraint (irdxc3):                                            |", sep="\n")

      default_redox_minus3 <- "\n|  [ ] (-3) Use O2(g) line in the aqueous basis species block                  |"
      default_redox_minus2 <-   "|  [ ] (-2) pe (pe units)      | 0.00000E+00| (pei)                            |"
      default_redox_minus1 <-   "|  [ ] (-1) Eh (volts)         | 0.00000E+00| (ehi)                            |"
      default_redox_0      <-   "|  [ ] ( 0) Log fO2 (log bars) | 0.00000E+00| (fo2lgi)                         |"
      default_redox_1      <-   "|  [ ] ( 1) Couple (aux. sp.)  |None                    | (uredox)             |"


      if(df[row, "redox_flag"] == -3){
        default_redox_minus3 <- "\n|  [x] (-3) Use O2(g) line in the aqueous basis species block                  |"
      } else if(df[row, "redox_flag"] == -2){
        default_redox_minus2 <- paste0("|  [x] (-2) pe (pe units)      |",
                              format(df[row, "redox_value"], width=12, justify="right"),
                              "| (pei)                            |")
      } else if(df[row, "redox_flag"] == -1){
        default_redox_minus1 <- paste0("|  [x] (-1) Eh (volts)         |",
                              format(df[row, "redox_value"], width=12, justify="right"),
                              "| (ehi)                            |")
      } else if(df[row, "redox_flag"] == 0){
        default_redox_0 <- paste0("|  [x] ( 0) Log fO2 (log bars) |",
                              format(df[row, "redox_value"], width=12, justify="right"),
                              "| (fo2lgi)                         |")
      } else if(df[row, "redox_flag"] == 1){
        default_redox_1 <- paste0("|  [x] ( 1) Couple (aux. sp.)  | ",
                              format(redox_aux, width=23, justify="left"),
                              "| (uredox)             |")
      } else {
        stop(paste0("Error when writing .3i file for sample ", df[row, "Sample"],
                   ". Redox flag was not recognized! Choose -3, -2, -1, 0, or 1."))
      }
      redox_block <- paste(default_redox_minus3, default_redox_minus2,
                           default_redox_minus1, default_redox_0,
                           default_redox_1, sep = "\n")

      eq3.header6 <- paste("\n|------------------------------------------------------------------------------|",
                             "|Aqueous Basis Species/Constraint Species        |Conc., etc. |Units/Constraint|",
                             "| (uspeci(n)/ucospi(n))                          | (covali(n))|(ujf3(jflgi(n)))|",
                             "|------------------------------------------------------------------------------|",
                             sep = "\n")

      # handle aqueous block
      aqueous_lines <- c()
      for(column in names(df)){
        if(suppress_missing & is.na(df[row, column])){
          df[row, column] <- 0
        }
        if(!is.na(df[row, column])){
          species_name  <- header_species(column)
          if(!(species_name %in% exclude)){
           species_value <- df[row, column]
           species_unit  <- header_unit(column, keepcase=T)
           if(!(species_unit %in% EQ3_jflags)){
             if(tolower(species_unit) == "ppb"){
               species_value <- species_value/1000
               species_unit <- "mg/L"
             } else if(tolower(species_unit) == "ppm"){
               species_unit <- "mg/L"
             } else {
               vprint(paste("Error creating .3i file:", species_unit,
               "is not a recognized aqueous block jflag. Try checking",
               "capitalization and spelling to match one of the following:",
               EQ3_jflags), verbose=verbose)
             }
           }
           species_value <- format(sprintf("%.5E", as.numeric(species_value)), width=12, justify="right")
           this_aq_line <- paste0("\n|",
                                  format(species_name, width=48), "|",
                                  species_value,  "|",
                                  format(species_unit, width=16),   "|")
           aqueous_lines <- c(aqueous_lines, this_aq_line)
          }
        }
      }

      # suppressing species
      for(species in suppress){
        this_aq_line <- paste0("\n|",
                               format(species, width=48), "| ",
                               sprintf("%.5E", 0),  "|",
                               format("Suppressed", width=16),   "|")
        aqueous_lines <- c(aqueous_lines, this_aq_line)
      }

      aqueous_block <- paste(aqueous_lines, collapse="")

      eq3.ender <- paste("\n|------------------------------------------------------------------------------|",
    "* Valid jflag strings (ujf3(jflgi(n))) are:                                    *",
    "*    Suppressed          Molality            Molarity                          *",
    "*    mg/L                mg/kg.sol           Alk., eq/kg.H2O                   *",
    "*    Alk., eq/L          Alk., eq/kg.sol     Alk., mg/L CaCO3                  *",
    "*    Alk., mg/L HCO3-    Log activity        Log act combo                     *",
    "*    Log mean act        pX                  pH                                *",
    "*    pHCl                pmH                 pmX                               *",
    "*    Hetero. equil.      Homo. equil.        Make non-basis                    *",
    "*------------------------------------------------------------------------------*",
    "|Create Ion Exchangers  | (net)                                                |",
    "|------------------------------------------------------------------------------|",
    "|Advisory: no exchanger creation blocks follow on this file.                   |",
    "|Option: on further processing (writing a PICKUP file or running XCON3 on the  |",
    "|present file), force the inclusion of at least one such block (qgexsh):       |",
    "|  [ ] (.true.)                                                                |",
    "|------------------------------------------------------------------------------|",
    "|Ion Exchanger Compositions      | (neti)                                      |",
    "|------------------------------------------------------------------------------|",
    "|Exchanger phase |None                    | (ugexpi(n))                        |",
    "|------------------------------------------------------------------------------|",
    "|->|Moles/kg.H2O    |  0.0000    | (cgexpi(n))                                 |",
    "|------------------------------------------------------------------------------|",
    "|->|Exchange site   |None    | (ugexji(j,n))                                   |",
    "|------------------------------------------------------------------------------|",
    "|--->|Exchange species        |Eq. frac.   | (this is a table header)          |",
    "|------------------------------------------------------------------------------|",
    "|--->|None                    | 0.00000E+00| (ugexsi(i,j,n), egexsi(i,j,n))    |",
    "|------------------------------------------------------------------------------|",
    "|Solid Solution Compositions     | (nxti)                                      |",
    "|------------------------------------------------------------------------------|",
    "|Solid Solution          |None                    | (usoli(n))                 |",
    "|------------------------------------------------------------------------------|",
    "|->|Component               |Mole frac.  | (this is a table header)            |",
    "|------------------------------------------------------------------------------|",
    "|->|None                    | 0.00000E+00| (umemi(i,n), xbari(i,n))            |",
    "|------------------------------------------------------------------------------|",
    "|Alter/Suppress Options  | (nxmod)                                             |",
    "|------------------------------------------------------------------------------|",
    "|Species                                         |Option          |Alter value |",
    "| (uxmod(n))                                     |(ukxm(kxmod(n)))| (xlkmod(n))|",
    "|------------------------------------------------------------------------------|",
    "|None                                            |None            | 0.00000E+00|",
    "|------------------------------------------------------------------------------|",
    "* Valid alter/suppress strings (ukxm(kxmod(n))) are:                           *",
    "*    Suppress            Replace             AugmentLogK                       *",
    "*    AugmentG                                                                  *",
    "*------------------------------------------------------------------------------*",
    "|Iopt Model Option Switches (\"( 0)\" marks default choices)                     |",
    "|------------------------------------------------------------------------------|",
    "|iopt(4) - Solid Solutions:                                                    |",
    "|  [x] ( 0) Ignore                                                             |",
    "|  [ ] ( 1) Permit                                                             |",
    "|------------------------------------------------------------------------------|",
    "|iopt(11) - Auto Basis Switching in pre-N-R Optimization:                      |",
    "|  [x] ( 0) Turn off                                                           |",
    "|  [ ] ( 1) Turn on                                                            |",
    "|------------------------------------------------------------------------------|",
    "|iopt(17) - PICKUP File Options:                                               |",
    "|  [ ] (-1) Don't write a PICKUP file                                          |",
    "|  [x] ( 0) Write a PICKUP file                                                |",
    "|------------------------------------------------------------------------------|",
    "|iopt(19) - Advanced EQ3NR PICKUP File Options:                                |",
    "|  [x] ( 0) Write a normal EQ3NR PICKUP file                                   |",
    "|  [ ] ( 1) Write an EQ6 INPUT file with Quartz dissolving, relative rate law  |",
    "|  [ ] ( 2) Write an EQ6 INPUT file with Albite dissolving, TST rate law       |",
    "|  [ ] ( 3) Write an EQ6 INPUT file with Fluid 1 set up for fluid mixing       |",
    "|------------------------------------------------------------------------------|",
    "|Iopg Activity Coefficient Option Switches (\"( 0)\" marks default choices)      |",
    "|------------------------------------------------------------------------------|",
    "|iopg(1) - Aqueous Species Activity Coefficient Model:                         |",
    "|  [ ] (-1) The Davies equation                                                |",
    "|  [x] ( 0) The B-dot equation                                                 |",
    "|  [ ] ( 1) Pitzer's equations                                                 |",
    "|  [ ] ( 2) HC + DH equations                                                  |",
    "|------------------------------------------------------------------------------|",
    "|iopg(2) - Choice of pH Scale (Rescales Activity Coefficients):                |",
    "|  [ ] (-1) \"Internal\" pH scale (no rescaling)                                 |",
    "|  [x] ( 0) NBS pH scale (uses the Bates-Guggenheim equation)                  |",
    "|  [ ] ( 1) Mesmer pH scale (numerically, pH = -log m(H+))                     |",
    "|------------------------------------------------------------------------------|",
    "|Iopr Print Option Switches (\"( 0)\" marks default choices)                     |",
    "|------------------------------------------------------------------------------|",
    "|iopr(1) - Print All Species Read from the Data File:                          |",
    "|  [x] ( 0) Don't print                                                        |",
    "|  [ ] ( 1) Print                                                              |",
    "|------------------------------------------------------------------------------|",
    "|iopr(2) - Print All Reactions:                                                |",
    "|  [x] ( 0) Don't print                                                        |",
    "|  [ ] ( 1) Print the reactions                                                |",
    "|  [ ] ( 2) Print the reactions and log K values                               |",
    "|  [ ] ( 3) Print the reactions, log K values, and associated data             |",
    "|------------------------------------------------------------------------------|",
    "|iopr(3) - Print the Aqueous Species Hard Core Diameters:                      |",
    "|  [x] ( 0) Don't print                                                        |",
    "|  [ ] ( 1) Print                                                              |",
    "|------------------------------------------------------------------------------|",
    "|iopr(4) - Print a Table of Aqueous Species Concentrations, Activities, etc.:  |",
    "|  [ ] (-3) Omit species with molalities < 1.e-8                               |",
    "|  [ ] (-2) Omit species with molalities < 1.e-12                              |",
    "|  [ ] (-1) Omit species with molalities < 1.e-20                              |",
    "|  [x] ( 0) Omit species with molalities < 1.e-100                             |",
    "|  [ ] ( 1) Include all species                                                |",
    "|------------------------------------------------------------------------------|",
    "|iopr(5) - Print a Table of Aqueous Species/H+ Activity Ratios:                |",
    "|  [x] ( 0) Don't print                                                        |",
    "|  [ ] ( 1) Print cation/H+ activity ratios only                               |",
    "|  [ ] ( 2) Print cation/H+ and anion/H+ activity ratios                       |",
    "|  [ ] ( 3) Print ion/H+ activity ratios and neutral species activities        |",
    "|------------------------------------------------------------------------------|",
    "|iopr(6) - Print a Table of Aqueous Mass Balance Percentages:                  |",
    "|  [ ] (-1) Don't print                                                        |",
    "|  [x] ( 0) Print those species comprising at least 99% of each mass balance   |",
    "|  [ ] ( 1) Print all contributing species                                     |",
    "|------------------------------------------------------------------------------|",
    "|iopr(7) - Print Tables of Saturation Indices and Affinities:                  |",
    "|  [ ] (-1) Don't print                                                        |",
    "|  [x] ( 0) Print, omitting those phases undersaturated by more than 10 kcal   |",
    "|  [ ] ( 1) Print for all phases                                               |",
    "|------------------------------------------------------------------------------|",
    "|iopr(8) - Print a Table of Fugacities:                                        |",
    "|  [ ] (-1) Don't print                                                        |",
    "|  [x] ( 0) Print                                                              |",
    "|------------------------------------------------------------------------------|",
    "|iopr(9) - Print a Table of Mean Molal Activity Coefficients:                  |",
    "|  [x] ( 0) Don't print                                                        |",
    "|  [ ] ( 1) Print                                                              |",
    "|------------------------------------------------------------------------------|",
    "|iopr(10) - Print a Tabulation of the Pitzer Interaction Coefficients:         |",
    "|  [x] ( 0) Don't print                                                        |",
    "|  [ ] ( 1) Print a summary tabulation                                         |",
    "|  [ ] ( 2) Print a more detailed tabulation                                   |",
    "|------------------------------------------------------------------------------|",
    "|iopr(17) - PICKUP file format (\"W\" or \"D\"):                                   |",
    "|  [x] ( 0) Use the format of the INPUT file                                   |",
    "|  [ ] ( 1) Use \"W\" format                                                     |",
    "|  [ ] ( 2) Use \"D\" format                                                     |",
    "|------------------------------------------------------------------------------|",
    "|Iodb Debugging Print Option Switches (\"( 0)\" marks default choices)           |",
    "|------------------------------------------------------------------------------|",
    "|iodb(1) - Print General Diagnostic Messages:                                  |",
    "|  [x] ( 0) Don't print                                                        |",
    "|  [ ] ( 1) Print Level 1 diagnostic messages                                  |",
    "|  [ ] ( 2) Print Level 1 and Level 2 diagnostic messages                      |",
    "|------------------------------------------------------------------------------|",
    "|iodb(3) - Print Pre-Newton-Raphson Optimization Information:                  |",
    "|  [x] ( 0) Don't print                                                        |",
    "|  [ ] ( 1) Print summary information                                          |",
    "|  [ ] ( 2) Print detailed information (including the beta and del vectors)    |",
    "|  [ ] ( 3) Print more detailed information (including matrix equations)       |",
    "|  [ ] ( 4) Print most detailed information (including activity coefficients)  |",
    "|------------------------------------------------------------------------------|",
    "|iodb(4) - Print Newton-Raphson Iteration Information:                         |",
    "|  [x] ( 0) Don't print                                                        |",
    "|  [ ] ( 1) Print summary information                                          |",
    "|  [ ] ( 2) Print detailed information (including the beta and del vectors)    |",
    "|  [ ] ( 3) Print more detailed information (including the Jacobian)           |",
    "|  [ ] ( 4) Print most detailed information (including activity coefficients)  |",
    "|------------------------------------------------------------------------------|",
    "|iodb(6) - Print Details of Hypothetical Affinity Calculations:                |",
    "|  [x] ( 0) Don't print                                                        |",
    "|  [ ] ( 1) Print summary information                                          |",
    "|  [ ] ( 2) Print detailed information                                         |",
    "|------------------------------------------------------------------------------|",
    "|Numerical Parameters                                                          |",
    "|------------------------------------------------------------------------------|",
    "| Beta convergence tolerance      | 0.00000E+00| (tolbt)                       |",
    "| Del convergence tolerance       | 0.00000E+00| (toldl)                       |",
    "| Max. Number of N-R Iterations   |   0        | (itermx)                      |",
    "|------------------------------------------------------------------------------|",
    "|Ordinary Basis Switches (for numerical purposes only)    | (nobswt)           |",
    "|------------------------------------------------------------------------------|",
    "|Replace |None                                            | (uobsw(1,n))       |",
    "|   with |None                                            | (uobsw(2,n))       |",
    "|------------------------------------------------------------------------------|",
    "|Sat. flag tolerance     | 0.00000E+00| (tolspf)                               |",
    "|------------------------------------------------------------------------------|",
    "|Aq. Phase Scale Factor  | 1.00000E+00| (scamas)                               |",
    "|------------------------------------------------------------------------------|",
    "|End of problem                                                                |",
    "|------------------------------------------------------------------------------|", sep="\n")


      this_file <- paste0(eq3.header1, eq3.samplename, eq3.header2,
                          eq3.temperature, eq3.header3, eq3.density,
                          eq3.header4, eq3.cb_block, eq3.header5,
                          redox_block, eq3.header6, aqueous_block,
                          eq3.ender, collapse = "")


      write(this_file, eq3.filename, append=FALSE)

    } # end loop for .3i file writing

    # restore original working directory
    setwd(wd)
    vprint("Finished creating EQ3 input files in rxn_3i folder.", verbose=verbose)

    this_time <- proc.time() - ptm                             

    message(paste("Preprocessing done! Took", round(this_time["elapsed"], 1), "seconds."))
                                 
    return(df)
}