python download_era5.py

cdo -b 32 chname,t2m,tas -remapbil,netcdf/my_grid.txt -yearmean tas_era5.nc netcdf/tas_era5.nc
cdo -b 32 chname,tp,pr -chunit,m,mm/s -mulc,0.01157407407 -remapbil,netcdf/my_grid.txt -yearmean pr_era5.nc netcdf/pr_era5.nc

