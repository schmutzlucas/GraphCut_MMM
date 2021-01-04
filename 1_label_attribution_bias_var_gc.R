source("functions/gc_function.R")
source("functions/bias_var_function.R")
source("functions/mmm_function.R")
source("functions/mmmo_function.R")
source("functions/min_bias_function.R")

# Variable utilisée : "tas" ou "pr"
var <- "tas" 


load("input/models_list.rdata")
load(paste0("input/load_",var,"_cv.rdata")) # useless?

# Names of models used in the ensemble
model_names <- TAS_MODELS[c(1,3,4,5,7,8,11,12,13,16,17,20,22,24,26,30,31,33,34,37),2] 


ref <- get(paste0("tas_",model_names[[1]]))
var_present <- var_future <- array(0,c(nrow = nrow(ref),ncol = ncol(ref),length(model_names)))
for(i in 1:length(model_names)){
  var_present[,,i] <- get(paste0(var,"_",model_names[[i]]))
  var_future[,,i] <- get(paste0(var,"_",model_names[[i]],"_2100"))
}


# List of label attribution matrices obtained with GraphCut
mmcombi <- list()
for(m in 1:length(model_names)){
  mmcombi[[m]] <- list()
  #### GraphCut ####
  gc_result <- graphcut(
    ref.datacost = var_present[,, m],
    models.datacost = var_present[,, -m],
    models.smoothcost = var_present[,, -m],
    weight.data = 1,
    weight.smooth = 1
  )
  gc_bias_var <- bias_var(
    ref_future = var_future[,, m],
    var_future = var_future[,, -m],
    labeling = gc_result$`Label attribution`
  )
  mmcombi[[m]][["gc_present"]] <- c(gc_bias_var, gc_result) 

  gc_result <- graphcut(
    ref.datacost = var_present[,, m],
    models.datacost = var_present[,, -m],
    models.smoothcost = var_future[,, -m],
    weight.data = 1,
    weight.smooth = 1
  )
  gc_bias_var <- bias_var(
    ref_future = var_future[,, m],
    var_future = var_future[,, -m],
    labeling = gc_result$`Label attribution`
  )
  mmcombi[[m]][["gc_hybrid"]] <- c(gc_bias_var, gc_result)

  gc_result <- graphcut(
    ref.datacost = var_future[,, m],
    models.datacost = var_future[,, -m],
    models.smoothcost = var_future[,, -m],
    weight.data = 1,
    weight.smooth = 1
  )
  gc_bias_var <- bias_var(
    ref_future = var_future[,, m],
    var_future = var_future[,, -m],
    labeling = gc_result$`Label attribution`
  )
  mmcombi[[m]][["gc_future"]] <- c(gc_bias_var, gc_result)

  ##### Min bias #####
  minbias_result <- min_bias(
    ref = var_present[,, m],
    var = var_present[,, -m]
  )
  minbias_bias_var <- bias_var(
    ref_future = var_future[,, m],
    var_future = var_future[,, -m],
    labeling = minbias_result
    )
  mmcombi[[m]][["min_bias"]] <- c(minbias_bias_var, list("Label attribution" = minbias_result))

  ##### MMM #####
  mmcombi[[m]][["mmm"]] <- mmm(
    ref = var_future[,, m],
    var = var_future[,, -m]
  )

  ##### MMM optimized #####
  mmcombi[[m]][["mmmo"]] <- mmmo(
    ref_present = var_present[,, m],
    ref_future = var_future[,, m],
    var_present = var_present[,, -m],
    var_future = var_future[,, -m]
  )

}

names(mmcombi) <- model_names
saveRDS(
  mmcombi,
  file = paste0("output/", var, "_mmcombi_perfectmodel.rds")
)
saveRDS(
  var_present,
  file = paste0("output/", var, "_var_present_cmip5.rds")
)
saveRDS(
  var_future,
  file = paste0("output/", var, "_var_futur_cmip5.rds")
)

