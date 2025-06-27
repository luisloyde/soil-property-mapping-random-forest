
# Modelado y Mapeo Digital de Propiedades de Suelo a partir de Muestreo de Campo y DEM

## **1. Librerías y Configuración**

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

# Definir el directorio de trabajo como el del proyecto
setwd(here::here())


## **2. Cargar y Procesar Datos de Muestreo**

datos <- read.csv("muestreo_estandarizado.csv")
# Convertir a objeto sf con CRS adecuado (ajusta según corresponda)
datos_sf <- st_as_sf(datos, coords = c("X", "Y"), crs = 32614)

## **3. Cargar y Procesar DEM**

dem <- rast("dem.tif")
crs(dem) <- "EPSG:32614"

## **4. Derivar Variables Topográficas**

# Calcular pendiente a partir del DEM
pendiente <- terrain(dem, v = "slope", unit = "degrees")
variables_topograficas <- c(dem, pendiente)
names(variables_topograficas) <- c("altitud", "pendiente")

# Extraer valores de altitud y pendiente para cada punto de muestreo
valores_extraidos <- extract(variables_topograficas, vect(datos_sf))
datos <- cbind(datos, valores_extraidos)


## **5. Definir Variables y Etiquetas**


variables <- c("Profundidad", "DAP", "CC", "PMP", "Arena", "Limo", "Arcilla", "Arenas_finas", "CO", "MO", "Cond_Hidr")
etiquetas <- c(
  "Profundidad (cm)",
  "Densidad Aparente Promedio (g/cm³)",
  "Capacidad de Campo (%HG)",
  "Punto de Marchitez Permanente (%HG)",
  "Arena (%)",
  "Limo (%)",
  "Arcilla (%)",
  "Arenas finas (g)",
  "Carbono orgánico (%)",
  "Materia orgánica (%)",
  "Conductividad hidráulica K (mm/h)"
)

variables <- c("SOL_BD","SOL_AWC","SOL_K","SOL_CBN","SOL_CLAY","SOL_SILT","SOL_SAND","USLE_K")
etiquetas <- c(
  expression("Densidad Aparente Promedio (g/cm"^3*")"),
  expression("Capacidad de agua disponible (mm[H[2]*O]/mm[suelo]"),
  expression("Conductividad hidráulica K (mm/h)"),
  "Contenido de carbono orgánico (%)",
  "Arcilla (%)",
  "Limo (%)",
  "Arena (%)",
  "Factor de erosibilidad del suelo"
)

names(etiquetas) <- variables

estratos <- unique(datos$Estrato_cm)


## **6. Modelado y Generación de Mapas Raster**


archivos_raster <- list()
nombres_mapa <- list()
etiquetas_mapa <- list()
limite <- st_read("limite_v2.shp")
limite_sf <- vect(limite)

for (var in variables) {
  cat("Procesando la variable:", var, "\n")
  for (estrato in estratos) {
    cat("Procesando el estrato:", estrato, "para la variable:", var, "\n")
    datos_estrato <- datos %>% filter(Estrato_cm == estrato)
    datos_estrato <- datos_estrato %>%
      filter(!is.na(.data[[var]]), !is.na(altitud), !is.na(pendiente))
    if (nrow(datos_estrato) < 10) {
      cat("No hay suficientes datos para el estrato", estrato, "y variable", var, "\n")
      next
    }
    set.seed(123)
    train_index <- createDataPartition(datos_estrato[[var]], p = 0.5, list = FALSE)
    datos_train <- datos_estrato[train_index, ]
    datos_test <- datos_estrato[-train_index, ]
    modelo_rf <- randomForest(as.formula(paste(var, "~ altitud + pendiente")), 
                              data = datos_train, importance = TRUE)
    predicciones_test <- predict(modelo_rf, newdata = datos_test)
    mse <- mean((datos_test[[var]] - predicciones_test)^2)
    cat("MSE:", mse,"\n")
    cat("RMSE:", sqrt(mse), "\n")
    
    # Predicción y recorte directo
    mapa_var <- predict(variables_topograficas, modelo_rf)
    mapa_var_clip <- mask(crop(mapa_var, limite_sf), limite_sf)
    
    output_clip_filename <- paste0("mapa_", tolower(var), "_", gsub("-", "_", estrato), "_clip.tif")
    writeRaster(mapa_var_clip, output_clip_filename, overwrite = TRUE)
    
    archivos_raster[[length(archivos_raster) + 1]] <- output_clip_filename
    nombres_mapa[[length(nombres_mapa) + 1]] <- paste(var, estrato, sep = " - ")
    etiquetas_mapa[[length(etiquetas_mapa) + 1]] <- paste0(etiquetas[var], "\n", estrato)
  }
}


## **8. Clasificación Edafológica y Digitalización Vectorial**

# Apilar rásteres recortados (cada banda representa una variable/estrato)
stack <- rast(unlist(archivos_raster))
names(stack) <- variables

# Leer shapefile de unidades edafológicas y simplificarlas
unidades <- st_read("suelo.shp")
unidades_simple <- unidades %>%
  group_by(GRUPO1) %>%  #PUEDE CAMBIARSE A GRUPO1
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  filter(!(GRUPO1 %in% c("NA")))

# Pila todos los raster, asegurando mismo extent/resolución
stack_full <- c(stack, dem)
# Convierte a tabla, pero NO elimina NAs
pixel_data <- as.data.frame(stack_full, na.rm=FALSE)
# Usar píxeles completos para clustering
valid_rows <- complete.cases(pixel_data)
set.seed(123)
k=nrow(unidades_simple)
kmeans_res <- kmeans(pixel_data[valid_rows, ], centers = k)
# Asignar los clusters de vuelta como arriba
cluster_map <- rep(NA, nrow(pixel_data))
cluster_map[valid_rows] <- kmeans_res$cluster
r_class <- rast(stack_full, nlyr=1)
values(r_class) <- cluster_map
# 5. Suavizado y exportación
r_class_smooth <- focal(r_class, w=9, fun=mean, na.policy="omit", na.rm=TRUE) #modificar w y fun
plot(r_class_smooth)
plot(unidades_simple, add = TRUE, border = "red", col = NA)

unidades_edafo <- as.polygons(r_class_smooth)                  # Convierte a polígonos
plot(unidades_edafo)


writeRaster(r_class_smooth, "kmeans_clusters_suavizado.tif", overwrite=TRUE)












stack_full <- c(stack, dem_clip)
names(stack_full)[nlyr(stack_full)] <- "DEM"

# 3. Extrae valores de todos los píxeles, elimina NAs
pixel_data <- as.data.frame(stack_full, na.rm = TRUE)

# 4. Realiza k-means
set.seed(123)
k <- 15  # El número de clusters que prefieras
kmeans_res <- kmeans(pixel_data, centers = k)

# 5. Crear raster de clusters
r_class <- rast(stack_full, nlyr = 1)
values(r_class) <- NA
valid_idx <- which(complete.cases(as.data.frame(stack_full)))
values(r_class)[valid_idx] <- kmeans_res$cluster

# 6. (Opcional) Suavizado usando filtro de mayoría
r_class_smooth <- focal(r_class, w = 3, fun = modal, na.policy = "omit", na.rm=TRUE)

# 7. Guarda el raster final
writeRaster(r_class_smooth, "kmeans_clusters_suavizado.tif", overwrite=TRUE)
writeRaster(r_class, "kmeans_clusters_bruto.tif", overwrite=TRUE)

## **10. Exportar imágenes de mapas**

profundidades <- unique(estratos)
propiedades <- variables

plots_list <- list()
for (i in seq_along(archivos_raster)) {
  raster_file <- archivos_raster[[i]]
  etiqueta_mapa <- etiquetas_mapa[[i]]
  raster_data <- rast(raster_file)
  df <- as.data.frame(raster_data, xy = TRUE)
  colnames(df)[3] <- "valor"
  p <- ggplot(df, aes(x = x, y = y, fill = valor)) +
    geom_raster() +
    scale_fill_viridis(option = "D", na.value = "white") +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
    ) +
    ggtitle(etiqueta_mapa)
  plots_list[[i]] <- p
}

n_filas <- length(propiedades)
n_columnas <- length(profundidades)
grid_plots <- vector("list", n_filas * n_columnas)
for (i in seq_along(nombres_mapa)) {
  name_parts <- strsplit(nombres_mapa[[i]], " - ")[[1]]
  prop_idx <- match(name_parts[1], propiedades)
  prof_idx <- match(name_parts[2], profundidades)
  grid_plots[[(prop_idx - 1) * n_columnas + prof_idx]] <- plots_list[[i]]
}

ancho_px <- 1000 * n_columnas
alto_px  <- 750 * n_filas

jpeg("mapas_propiedades_grid_etiquetas.jpeg", width = ancho_px, height = alto_px, res = 300)
do.call("grid.arrange", c(grid_plots, nrow = n_filas, ncol = n_columnas))
dev.off()

png("mapas_propiedades_grid_etiquetas.png", width = ancho_px, height = alto_px, res = 300)
do.call("grid.arrange", c(grid_plots, nrow = n_filas, ncol = n_columnas))
dev.off()

