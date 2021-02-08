#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

var="pr"

# Rscript  --vanilla  1_label_attribution_bias_var_gc.R $var
Rscript  --vanilla  2_data_for_plots_perfect_model.R $var
Rscript  --vanilla  3_plot_graphs_maps_perfert_model.R $var
# Rscript  --vanilla  4_label_attribution_era5.R $var
Rscript  --vanilla  5_data_for_plots_era5.R $var
Rscript  --vanilla  6_plot_graphs_maps_era5.R $var
Rscript  --vanilla  7_projection_era5.R $var
