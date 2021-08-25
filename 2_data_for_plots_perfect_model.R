source("functions/gradients.R")


# A refaire et à connecter avec le script précédent
args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
var = args[1]
print(var)


var_future <- readRDS(
  file = paste0("output/", var, "_var_futur_perfectmodel.rds")
)

mmcombi <- readRDS(file = paste0("output/", var, "_mmcombi_perfectmodel.rds"))
model_names <- names(mmcombi)


bias_map <- lapply(mmcombi, function(m){
  lapply(m, function(mmc){
    mmc$Bias
  })
})

mae_bias_map <- lapply(mmcombi, function(m){
  lapply(m, function(mmc){
    abs(mmc$Bias)
  })
})

mae_bias <- sapply(mmcombi, function(m){
  sapply(m, function(mmc){
    mean(abs(mmc$Bias))
  })
})



mae_gradient_map <- lapply(seq_along(mmcombi), function(m){
  lapply(mmcombi[[m]], function(mmc){
    gradient_mae(var_future[,, m], mmc$Var)
  })
})

mae_gradient <- sapply(mae_gradient_map, function(m){
  sapply(m, function(mmc){
    mean(mmc)
  })
})



#### Summary

mae_gradient <- cbind(mae_gradient, "summary" = apply(mae_gradient, 1, mean))


mae_gradient_map[["summary"]] <- lapply(seq_along(mae_gradient_map[[1]]), function(i){
  apply(
    array(
      unlist(lapply(mae_gradient_map,  function(map){
        map[[i]]
      })),
      dim = dim(var_future)
    ),
    1:2, mean
  )
})

mae_bias <- cbind(mae_bias, "summary"  = apply(mae_bias, 1, mean))

mae_bias_map[["summary"]] <- lapply(seq_along(mmcombi[[1]]), function(i){
  apply(
    array(
      unlist(lapply(mmcombi,  function(map){
        abs(map[[i]]$Bias)
      })),
      dim = dim(var_future)
    ),
    1:2, mean
  )
})

saveRDS(mae_bias,  paste0("output/", var, "_mae_bias_perfectmodel.rds"))
saveRDS(mae_bias_map,  paste0("output/", var, "_mae_bias_map_perfectmodel.rds"))
saveRDS(bias_map,  paste0("output/", var, "_bias_map_perfectmodel.rds"))
saveRDS(mae_gradient,  paste0("output/", var, "_mae_gradient_perfectmodel.rds"))
saveRDS(mae_gradient_map,  paste0("output/", var, "_mae_gradient_map_perfectmodel.rds"))
