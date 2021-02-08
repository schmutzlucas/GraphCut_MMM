library(ncdf4)
library(fields)
library(stringr)

cdo_extract <- function(
  var,
  model,
  model_version,
  grid
){

  exercise_output <- "hist_rcp"

  path_input_hist <- paste("/bdd/CMIP5/output", model, model_version, "historical/mon/atmos/Amon/r1i1p1/latest", var, sep="/")
  path_input_rcp <- paste("/bdd/CMIP5/output", model, model_version, "rcp85/mon/atmos/Amon/r1i1p1/latest", var, sep="/")
  name_file <- paste(var, freq, model_version, exercise_output, run, sep="_")
  name_file <- paste0(name_file, ".nc")
  list_name_file_hist <- list.files(path_input_hist, full.names = TRUE)
  list_name_file_rcp <- list.files(path_input_rcp, full.names = TRUE)
  print(path_input_rcp)
  print(path_input_hist)
  print(list_name_file_hist)
  print(list_name_file_rcp)

  cmd_cdo <- ""
  for(file in c(list_name_file_hist, list_name_file_rcp)){
    cmd_cdo <- paste0(cmd_cdo,paste0("-remapbil,",grid," -yearmean ", file," "))
  }
  cmd_cdo <- paste0(
      "cdo mergetime ",
      cmd_cdo,
      name_file
  )

  print(cmd_cdo)
  system(cmd_cdo)

}
