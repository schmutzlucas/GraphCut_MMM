#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

pattern="../netcdf/pr*.nc"

echo $pattern

for file in $pattern
do
    if [ "$file" != "../netcdf/pr_era5.nc" ]; then
    echo "processing $file"
    cdo -b 32  -chunit,"kg m-2 s-1","mm/d" -mulc,86400 $file tmp.nc
    mv tmp.nc $file
    fi
done

cdo -b 32  -chunit,"mm/s","mm/d" -mulc,86400 ../netcdf/pr_era5.nc tmp.nc
mv tmp.nc pr_era5.nc
