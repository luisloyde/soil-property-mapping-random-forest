# Digital Soil Property Mapping from Field Sampling and DEM

## 1. Libraries and Configuration

rm(list = ls())

library(terra)
library(sf)
library(dplyr)
library(randomForest)
library(ggplot2)
library(caret)
library(viridis)
library(gridExtra)
library(exactextractr)

# Set working directory to project root
setwd(here::here())

## 2. Load and Process Field Sampling Data

data <- read.csv("muestreo_estandarizado.csv")
# Convert to sf object with appropriate CRS (preferably UTM)
data_sf <- st_as_sf(data, coords = c("X", "Y"), crs = 32614)

## 3. Load and Process DEM

dem <- rast("dem_wgs84.tif")
crs(dem) <- "EPSG:4326" # Preferably reproject to WGS84

## 4. Derive Topographic Variables

# Calculate slope from DEM
slope <- terrain(dem, v = "slope", unit = "degrees")
topo_vars <- c(dem, slope)
names(topo_vars) <- c("elevation", "slope")

# Extract elevation and slope values for each sampling point
extracted_vals <- extract(topo_vars, vect(data_sf))
data <- cbind(data, extracted_vals)

## 5. Define Variables and Labels

variables <- c("SOL_BD","SOL_AWC","SOL_K","SOL_CBN","SOL_CLAY","SOL_SILT","SOL_SAND","USLE_K")
labels <- c(
  expression("Bulk Density (g/cm"^3*")"),
  expression("Available Water Capacity (mm"[H[2]*O] * "/mm"[soil] * ")"),
  "Hydraulic Conductivity K (mm/h)",
  "Organic Carbon Content (%)",
  "Clay (%)",
  "Silt (%)",
  "Sand (%)",
  "USLE K Factor"
)
names(labels) <- variables

depths <- unique(data$Estrato_cm)

## 6. Modeling and Raster Map Generation

raster_files <- list()
map_names <- list()
map_labels <- list()
boundary <- st_read("limite_wgs84.shp") # Should be in CRS 4326 (WGS84)
boundary_sf <- vect(boundary)

for (var in variables) {
  cat("Processing variable:", var, "\n")
  for (depth in depths) {
    cat("Processing depth:", depth, "for variable:", var, "\n")
    data_depth <- data %>% filter(Estrato_cm == depth)
    data_depth <- data_depth %>%
      filter(!is.na(.data[[var]]), !is.na(elevation), !is.na(slope))
    if (nrow(data_depth) < 10) {
      cat("Not enough data for depth", depth, "and variable", var, "\n")
      next
    }
    set.seed(123)
    train_index <- createDataPartition(data_depth[[var]], p = 0.5, list = FALSE)
    train_data <- data_depth[train_index, ]
    test_data <- data_depth[-train_index, ]
    rf_model <- randomForest(as.formula(paste(var, "~ slope + elevation")), 
                             data = train_data, importance = TRUE)
    test_pred <- predict(rf_model, newdata = test_data)
    mse <- mean((test_data[[var]] - test_pred)^2)
    cat("MSE:", mse,"\n")
    cat("RMSE:", sqrt(mse), "\n")
    
    # Prediction and direct clipping
    var_map <- predict(topo_vars, rf_model)
    var_map_clip <- mask(crop(var_map, boundary_sf), boundary_sf)
    
    output_clip_filename <- paste0("map_", var, "_", gsub("-", "_", depth), "_clip.tif")
    writeRaster(var_map_clip, output_clip_filename, overwrite = TRUE)
    
    # Generate PNG/JPG map using ggplot2
    df_raster <- as.data.frame(var_map_clip, xy = TRUE)
    colnames(df_raster)[3] <- "Value"
    
    # Mathematical label as title
    label_title <- as.expression(substitute(e * "\n" * d, list(e = labels[[var]], d = depth)))
    
    p <- ggplot(df_raster, aes(x = x, y = y, fill = Value)) +
      geom_raster() +
      scale_fill_viridis(option = "D", na.value = "white") +
      coord_fixed() +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
      ) +
      ggtitle(label_title)
    
    # Save image file
    ggsave(sprintf("map_%s_%s.png", var, gsub("-", "_", depth)), p, width = 5, height = 3.75, dpi = 300)
    
    raster_files[[length(raster_files) + 1]] <- output_clip_filename
    map_names[[length(map_names) + 1]] <- paste(var, depth, sep = " - ")
  }
}

## 8. Edaphological Classification and Vectorization

# Stack clipped rasters (each band is a variable/depth)
stack <- rast(unlist(raster_files))
names(stack) <- variables

# Read soil units shapefile and simplify
soil_units <- st_read("suelo_wgs84.shp")
soil_units_simple <- soil_units %>%
  group_by(Grupo1) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  filter(!(Grupo1 %in% c("NA")))

# Stack alignment to DEM
stack_aligned <- resample(stack, dem, method = "bilinear")
stack_full <- c(stack_aligned, dem, slope)

# Convert to data frame, keep NAs
pixel_data <- as.data.frame(stack_full, na.rm=FALSE)

# Use complete pixels for clustering
valid_rows <- complete.cases(pixel_data)
set.seed(123)
k = nrow(soil_units_simple)

kmeans_res <- kmeans(pixel_data[valid_rows, ], centers = k)
# Assign clusters back to raster
cluster_map <- rep(NA, nrow(pixel_data))
cluster_map[valid_rows] <- kmeans_res$cluster
r_class <- rast(stack_full, nlyr=1)
values(r_class) <- cluster_map

# Smoothing
r_class_smooth <- focal(r_class, w=9, fun=mean, na.policy="omit", na.rm=TRUE) # adjust w and fun as needed

# Calculate median for each polygon
units_class <- as.polygons(r_class_smooth, dissolve = TRUE)
medians_shp <- extract(stack_full, units_class, fun = median, na.rm = TRUE)
units_class <- cbind(units_class, medians_shp[,-1])

# Link each unit to SWAT database
library(here)
library(DBI)
library(odbc)

file_path <- here::here("SWAT2012.mdb")

con <- dbConnect(
  odbc(),
  Driver = "Microsoft Access Driver (*.mdb, *.accdb)",
  DBQ = file_path
)

swat_df <- dbReadTable(con, "usersoil")

swat_medians <- swat_df %>%
  rowwise() %>%
  mutate(
    SOL_BD   = mean(c(SOL_BD1,   SOL_BD2),   na.rm = TRUE),
    SOL_AWC  = mean(c(SOL_AWC1,  SOL_AWC2),  na.rm = TRUE),
    SOL_K    = mean(c(SOL_K1,    SOL_K2),    na.rm = TRUE),
    SOL_CBN  = mean(c(SOL_CBN1,  SOL_CBN2),  na.rm = TRUE),
    SOL_CLAY = mean(c(CLAY1,     CLAY2),     na.rm = TRUE),
    SOL_SILT = mean(c(SILT1,     SILT2),     na.rm = TRUE),
    SOL_SAND = mean(c(SAND1,     SAND2),     na.rm = TRUE),
    USLE_K   = mean(c(USLE_K1,   USLE_K2),   na.rm = TRUE)
  ) %>%
  select(SEQN, SNAM, SOL_BD, SOL_AWC, SOL_K, SOL_CBN, SOL_CLAY, SOL_SILT, SOL_SAND, USLE_K)

# Calculate Euclidean distances to assign most similar soil type

library(purrr)

# Extract attributes as data.frame
df_units <- as.data.frame(units_class)

# Select only numeric columns needed
df_units_num <- df_units[variables]

# Ensure numeric
df_units_num <- data.frame(lapply(df_units_num, function(x) as.numeric(as.character(x))))

mat_units <- as.matrix(df_units_num)
mat_swat <- as.matrix(swat_medians[variables])

closest_indices <- apply(mat_units, 1, function(x) {
  dists <- apply(mat_swat, 1, function(y) sqrt(sum((x - y)^2)))
  which.min(dists)
})

# Assign most similar soils
df_units$SEQN_match <- swat_medians$SEQN[closest_indices]
df_units$SNAM_match <- swat_medians$SNAM[closest_indices]

# Add fields to SpatVector
units_class$SEQN <- df_units$SEQN_match
units_class$SNAM <- df_units$SNAM_match
writeVector(units_class, "edafologia_SWAT.shp", overwrite=TRUE)

# Convert SpatVector to raster
ref_raster <- dem
r_seqn <- rasterize(units_class, ref_raster, field = "SEQN")
writeRaster(r_seqn, "edafologia_SWAT_SEQN.tif", overwrite = TRUE)

# Generate .txt lookup table
swat_table <- unique(data.frame(
  Value = units_class$SEQN,
  Name  = units_class$SNAM
))
swat_table <- swat_table[order(swat_table$Value), ]
write.table(swat_table, "SWAT_soil_table.txt", sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)