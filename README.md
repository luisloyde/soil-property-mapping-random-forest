# Soil Property Mapping using Random Forest

This project provides an R pipeline for spatial modeling and mapping of soil properties using Random Forest machine learning, based on standardized field sampling, DEM-derived topographic variables, soil geospatial data from Mexico, and the Harmonized World Soil Database (FAO). The workflow includes model calibration, raster prediction, clipping to study area boundaries, grid visualization, and exporting zonal statistics by soil unit.

## Features

- Reads field data with coordinates and soil properties.
- Imports and processes a DEM to extract topographic variables (elevation, slope).
- Fits Random Forest models for each soil property and depth stratum.
- Predicts continuous raster maps for each combination.
- Clips output rasters to a study area boundary.
- Generates grid visualizations with custom labels.
- Calculates polygonal medians per edaphological unit and exports as Shapefile/GeoPackage.

## Requirements

- R (>= 4.0)
- R packages:
    - `terra`
    - `sf`
    - `dplyr`
    - `randomForest`
    - `ggplot2`
    - `caret`
    - `viridis`
    - `here`
    - `gridExtra`
    - `exactextractr`
    - `DBI`, `odbc` (for linking to Access/SQL databases)

Install dependencies in R:
```r
install.packages(c("terra", "sf", "dplyr", "randomForest", "ggplot2", "caret", "viridis", "here", "gridExtra", "exactextractr", "DBI", "odbc"))
```

## Usage

1. Place your standardized sampling file (`muestreo_estandarizado.csv`), DEM (`dem_wgs84.tif`), study area boundary (`limite_wgs84.shp`), and soil units shapefile (`suelo_wgs84.shp`) in the project folder.
2. Edit the script to match your data column names and CRSs if necessary.
3. Run the R script in RStudio or from the command line.
4. Output rasters, maps, and statistics will be generated in the working directory.

## Data format notes

- The CSV file must have columns for coordinates (`X`, `Y`), a column `Estrato_cm` for depth (e.g., `0-30 cm`), and columns for each soil property (see script for variable names).
- All spatial files must be in the same coordinate reference system (see `crs` settings in the script).

## License

This project is licensed under the GPL-3.0 License. See [LICENSE](LICENSE) for details.

## Disclaimer

This code is provided as-is, without warranty. Please check the compatibility of all R packages with your system.

## Outputs example

- Raster maps for each property and depth (GeoTIFF and PNG)
<img width="1067" height="770" alt="image" src="https://github.com/user-attachments/assets/6412a9b1-f4c6-4aeb-aeef-4cc288b9c5c0" />


- Clipped rasters to the study area boundary
- Shapefile of classified soil units with SWAT linkage
- Lookup table for SWAT soil codes
