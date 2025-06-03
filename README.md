# Soil Property Mapping using Random Forest

This project provides an R pipeline for spatial modeling and mapping of soil properties using Random Forest machine learning, based on standardized field sampling and DEM-derived topographic variables. The workflow includes model calibration, raster prediction, clipping to study area boundaries, grid visualization, and exporting zonal statistics by soil unit.

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
    - `this.path`
    - `gridExtra`
    - `exactextractr`

Install dependencies in R:
```r
install.packages(c("terra", "sf", "dplyr", "randomForest", "ggplot2", "caret", "viridis", "this.path", "gridExtra", "exactextractr"))
```

## Usage

1. Place your standardized sampling file (`muestreo_estandarizado.csv`), DEM (`dem.tif`), study area boundary (`limite.shp`), and soil units shapefile (`suelo.shp`) in the project folder.
2. Edit the script to match your data column names and CRS if necessary.
3. Run the R script in RStudio or from the command line.
4. Output rasters and statistics will be generated in the working directory.

## Data format notes

- The CSV file must have columns for coordinates (`X`, `Y`), a column `Estrato_cm` for depth (e.g., `0-30 cm`), and columns for each soil property (see script for variable names).
- All spatial files must be in the same CRS (see `crs` settings in the script).

## License

This project is licensed under the GPL-3.0 License. See [LICENSE](LICENSE) for details.

## Disclaimer

This code is provided as-is, without warranty. Please check the compatibility of all R packages with your system.

## Outputs example

![mapas_propiedades_grid_etiquetas](https://github.com/user-attachments/assets/bd3907c6-26a0-4178-8ba0-589eb52c492e)

