library(CHNOSZ)

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

fill_data0_head <- function(data0_template, db, grid_temps, grid_press, water_model){
    
  data0_template <- paste(data0_template, collapse="\n")
    
  # make sure the TP grid is in the correct form, esp. for single TP values
  grid_temps <- unlist(grid_temps)
  grid_press <- unlist(grid_press)
  grid_temps_original <- grid_temps
  if(length(grid_temps) == 1){
    grid_temps <- grid_temps + 0:7 # only the first T value is valid, but this is needed for EQ3
    grid_press <- rep(grid_press, 8)
  }
    
  # calculate debye huckel a and b parameters for the grid
  if(length(grid_temps_original) == 8){
    A_DH_grid <- unlist(water("A_DH", T=273.15+grid_temps, P=grid_press))
    B_DH_grid <- unlist(water("B_DH", T=273.15+grid_temps, P=grid_press)*10^-8)
  }else if(length(grid_temps_original) == 1){
    A_DH_grid <- unlist(water("A_DH", T=273.15+grid_temps[1], P=grid_press[1]))
    A_DH_grid <- c(A_DH_grid, rep(0, 7))
    B_DH_grid <- unlist(water("B_DH", T=273.15+grid_temps[1], P=grid_press[1])*10^-8)
    B_DH_grid <- c(B_DH_grid, rep(0, 7))
  }
      
  # format grid values
  grid_temps_f <- as.character(format(round(grid_temps, 4), nsmall = 4, scientific=F))
  grid_press_f <- as.character(format(round(grid_press, 4), nsmall = 4, scientific=F))
  A_DH_grid_f <- as.character(format(round(A_DH_grid, 4), nsmall = 4, scientific=F))
  B_DH_grid_f <- as.character(format(round(B_DH_grid, 4), nsmall = 4, scientific=F))

  # calculate bdot parameter
  if(length(grid_temps_original) == 1){
    bdot_grid_f <- c(as.character(format(round(calc_bdot(grid_temps[1]), 4), nsmall = 4, scientific=F)), rep("0.0000", 7))
  }else{
    bdot_grid_f <- as.character(format(round(calc_bdot(grid_temps), 4), nsmall = 4, scientific=F))
  }
    
  # cco2 (coefficients for the drummond (1981) polynomial)
  # GB note: might not change with T or P?
  # Examination of various data0 files seems to indicate that DBcreate does not change these values.

  # Calculate the "log k for eh reaction" grid.
  # From eq. 9 in EQPT (version 7.0) user manual, part 2, by Wolery:
  if(length(grid_temps_original) == 8){
    logK_Eh_vals <- subcrt(c("H2O", "O2", "e-", "H+"),
                           c(-2, 1, 4, 4),
                           c("liq", "gas", "aq", "aq"),
                           T=grid_temps,
                           P=round(grid_press, 9),
                           exceed.rhomin=TRUE,
                           exceed.Ttr=TRUE)$out$logK
  }else if(length(grid_temps_original) == 1){
    logK_Eh_vals <- subcrt(c("H2O", "O2", "e-", "H+"),
                           c(-2, 1, 4, 4),
                           c("liq", "gas", "aq", "aq"),
                           T=grid_temps[1],
                           P=round(grid_press[1], 9),
                           exceed.rhomin=TRUE,
                           exceed.Ttr=TRUE)$out$logK
    logK_Eh_vals <- c(logK_Eh_vals, rep(0, 7))
  }
      
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

  # insert pressure grid values into data0 template
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
  desc <- "data0.%s\nWater model: %s\nTP points: %s"
  min_desc <- "data0.min\nminimal working data0 file"
  data0_template <- sub(min_desc, sprintf(desc, db, water_model, length(grid_temps_original)), data0_template)

  return(data0_template) # returns a single string
    
}