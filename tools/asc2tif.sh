#!/bin/bash

# Convert .asc to .fig using gdal_translate

dir_in="/project/mugi/nas/PAPER2/CCLM-DCEP-Tree/ahf/asc"
dir_out="/project/mugi/nas/PAPER2/CCLM-DCEP-Tree/ahf/tif"

rm $dir_out/*

for f in $dir_in/*; do
	gdalwarp -tr 0.0025 0.0025 -r near "$dir_in/${f##*/}" "$dir_out/${f##*/}.tif"
done

