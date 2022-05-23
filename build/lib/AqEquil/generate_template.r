generate_template <- function(thermo_df, template_name, template_type, fixed_species){
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
    if(!(species %in% fixed_species)){
      input_template[[species]] <- c("Molality")
    }
  }
    
  write.csv(input_template, template_name, row.names=F)
}