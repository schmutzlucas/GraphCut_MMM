source("function_extract_hist_rcp.R")

models_list <- readRDS("models_list.rds")
head(models_list)


### Inputs
run <- "r1i1p1"
freq <- "Amon"
grid <- "my_grid.txt"

variables <- c("tas", "pr")

for(var in variables){
  for(m in 1:nrow(models_list)) {
    cdo_extract(var, models_list[m,  1], models_list[m, 2], grid)
  }
}

