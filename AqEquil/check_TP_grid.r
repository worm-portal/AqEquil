# load packages
suppressMessages({
  library(CHNOSZ)
  library(dplyr)
  library(comprehenr)
  library(stringr)
})

# print messages if 'verbose' setting >= vlevel of message.
vmessage <- function(m, vlevel, verbose){
  if(verbose >= vlevel){
    print(m)
  }
}

# function to round up the last digit of a four-digit number
# e.g., 52.6112 becomes 52.6113.
# This function is useful for reporting rounded PSAT pressures that keeps water in the liquid phase
roundup <- function(x, n){
  return(round(x, n) + 10^-n)
}

check_TP_grid <- function(grid_temps, grid_press, P1, water_model="SUPCRT92", check_for_errors=TRUE, verbose=1){

  grid_temps <- as.numeric(grid_temps)
    
  # set water model
  suppressMessages(water(water_model))
  
  # round grid temperatures to four decimal places
  grid_temps <- round(grid_temps, 4)

  # calculate PSAT pressure if specified by user or if pressure grid
  # has a number of values that does not equal temperature grid length.
  if(tolower(grid_press) == "psat"){
      vmessage("Calculating pressure grid along liquid-vapor saturation curve...", 2, verbose)
      grid_press <- water("Psat", T=grid_temps+273.15, P1=P1)[[1]]
  }else{
    grid_press <- as.numeric(grid_press)
  }

  # check TP polynomial
  if(length(grid_temps) == 8){
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
  }else{
    poly_coeffs_1 <- 'None'
    poly_coeffs_2 <- 'None'
  }

  if(check_for_errors){
    # check that water is a liquid at each T-P point
    if(length(grid_press) > 1){
      TP_grid_errors <- c()
      for(i in 1:length(grid_temps)){
                 
        tryCatch({
          psat_press <<- suppressMessages(water("Psat", T=grid_temps[i]+273.15, P1=P1)[[1]])
        }, error=function(e){
          psat_press <<- NA
        })
        if(is.na(psat_press)){
          tryCatch({
              rho_val <<- suppressMessages(water("rho", T=grid_temps[i]+273.15, P=grid_press[i], P1=P1)[[1]])
          }, error=function(e){
              rho_val <<- NA
          })
          if(is.na(rho_val)){
            TP_grid_errors <- c(TP_grid_errors,
                                paste("\nWater density could not be calculated at",
                                      grid_temps[i], "degrees C and",
                                      grid_press[i], "bar with the water model", water_model))
          }else if(rho_val <= 310){
            TP_grid_errors <- c(TP_grid_errors,
                                paste("\nWater density is", roundup(rho_val, 3), "kg m^3 at",
                                      grid_temps[i], "degrees C and",
                                      grid_press[i], "bar. This is too low (< 310 kg m^3)",
                                      "for this calculation. Increase pressure or decrease temperature."))
          }
        }else if(grid_press[i] < psat_press){
            
          TP_grid_errors <- c(TP_grid_errors,
                              paste("\n", grid_press[i], "bar is below liquid-vapor",
                                    "saturation pressure", roundup(psat_press, 4),
                                    "bar at", grid_temps[i], "degrees C."))
        }
      }
      if(length(TP_grid_errors) > 0){
        stop(paste(paste(TP_grid_errors, collapse="\n"),
                   "\n\nIncrease the pressure at these temperature points in 'grid_press'",
                   "to keep water in a liquid state."))
      }
      
  
        
    }
  }

  return(list("grid_temps"=grid_temps, "grid_press"=grid_press,
              "poly_coeffs_1"=poly_coeffs_1, "poly_coeffs_2"=poly_coeffs_2))
}