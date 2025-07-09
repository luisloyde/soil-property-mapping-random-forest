
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
limite <- st_read("limite.shp")
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
    modelo_rf <- randomForest(as.formula(paste(var, "~ pendiente + altitud")), 
                              data = datos_train, importance = TRUE)
    predicciones_test <- predict(modelo_rf, newdata = datos_test)
    mse <- mean((datos_test[[var]] - predicciones_test)^2)
    cat("MSE:", mse,"\n")
    cat("RMSE:", sqrt(mse), "\n")
    
    # Predicción y recorte directo
    mapa_var <- predict(variables_topograficas, modelo_rf)
    mapa_var_clip <- mask(crop(mapa_var, limite_sf), limite_sf)
    
    output_clip_filename <- paste0("mapa_", var, "_", gsub("-", "_", estrato), "_clip.tif")
    writeRaster(mapa_var_clip, output_clip_filename, overwrite = TRUE)
    
    archivos_raster[[length(archivos_raster) + 1]] <- output_clip_filename
    nombres_mapa[[length(nombres_mapa) + 1]] <- paste(var, estrato, sep = " - ")
  }
}


## **8. Clasificación Edafológica y Digitalización Vectorial**

# Apilar rásteres recortados (cada banda representa una variable/estrato)
stack <- rast(unlist(archivos_raster))
names(stack) <- variables

# Leer shapefile de unidades edafológicas y simplificarlas
unidades <- st_read("suelo.shp")
unidades_simple <- unidades %>%
  group_by(GRUPO1) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  filter(!(GRUPO1 %in% c("NA")))

# Pila todos los raster, asegurando mismo extent/resolución
stack_full <- c(stack, dem, pendiente)
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
# SuavizadO
r_class_smooth <- focal(r_class, w=9, fun=mean, na.policy="omit", na.rm=TRUE) #modificar w y fun

# Calcular mediana de cada polígono
unidades_class <- as.polygons(r_class_smooth, dissolve = TRUE)
medianas_shp <- extract(stack_full, unidades_class, fun = median, na.rm = TRUE)
unidades_class <- cbind(unidades_class, medianas_shp[,-1])

# Relacionar cada unidad a base de datos de SWAT
library(here)
library(DBI)
library(odbc)

ruta_archivo <- here::here("SWAT2012.mdb")

con <- dbConnect(
  odbc(),
  Driver = "Microsoft Access Driver (*.mdb, *.accdb)",
  DBQ = ruta_archivo
)

swat_df <- dbReadTable(con, "usersoil")

swat_medianas <- swat_df %>%
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

# Calcular distancias euclidianas para obtener tipo de suelo más parecido

library(purrr)

# Extrae los atributos como data.frame
df_unidades <- as.data.frame(unidades_class)

# Ahora selecciona solo las columnas numéricas que necesitas
df_unidades_num <- df_unidades[variables]

# Asegúrate que son numéricas
df_unidades_num <- data.frame(lapply(df_unidades_num, function(x) as.numeric(as.character(x))))

# Ahora ya puedes usar as.matrix
mat_unidades <- as.matrix(df_unidades_num)

mat_swat <- as.matrix(swat_medianas[variables])

closest_indices <- apply(mat_unidades, 1, function(x) {
  dists <- apply(mat_swat, 1, function(y) sqrt(sum((x - y)^2)))
  which.min(dists)
})

# Asignar los suelos más parecidos
df_unidades$SEQN_match <- swat_medianas$SEQN[closest_indices]
df_unidades$SNAM_match <- swat_medianas$SNAM[closest_indices]

# Añadir los campos al SpatVector
unidades_class$SEQN <- df_unidades$SEQN_match
unidades_class$SNAM <- df_unidades$SNAM_match
writeVector(unidades_class, "edafologia_SWAT.shp", overwrite=TRUE)

# Pasar spatVector a raster
library(terra)

ref_raster <- dem
r_seqn <- rasterize(unidades_class, ref_raster, field = "SEQN")
writeRaster(r_seqn, "edafologia_SWAT_SEQN.tif", overwrite = TRUE)









# 1. Define etiquetas como lista de expressions
etiquetas <- list(
  SOL_BD   = expression("Densidad Aparente Promedio (g/cm"^3*")"),
  SOL_AWC  = expression("Capacidad de agua disponible (mm[H[2]*O]/mm[suelo]"),
  SOL_K    = expression("Conductividad hidráulica K (mm/h)"),
  SOL_CBN  = expression("Contenido de carbono orgánico (%)"),
  SOL_CLAY = expression("Arcilla (%)"),
  SOL_SILT = expression("Limo (%)"),
  SOL_SAND = expression("Arena (%)"),
  USLE_K   = expression("Factor de erosibilidad del suelo")
)

# 2. Construye etiquetas_mapa como lista de expressions (no strings)
etiquetas_mapa <- list()
for (var in variables) {
  for (estrato in estratos) {
    etiquetas_mapa[[length(etiquetas_mapa) + 1]] <-
      as.expression(substitute(e * "\n" * p, list(e = etiquetas[[var]], p = estrato)))
  }
}

# 3. En el loop de gráficos, úsala directamente en ggtitle
plots_list <- list()
for (i in seq_along(archivos_raster)) {
  raster_file <- archivos_raster[[i]]
  etiqueta_mapa <- etiquetas_mapa[[i]]
  raster_data <- rast(raster_file)
  df <- as.data.frame(raster_data, xy = TRUE)
  colnames(df)[3] <- "Valor"
  p <- ggplot(df, aes(x = x, y = y, fill = Valor)) +
    geom_raster() +
    scale_fill_viridis(option = "D", na.value = "white") +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
    ) +
    ggtitle(etiqueta_mapa)  # <- ¡Aquí va como expression!
  plots_list[[i]] <- p
  ggsave(sprintf("mapa_%02d.png", i), p, width=5, height=3.75, dpi=300)
}










## **10. Exportar imágenes de mapas**

profundidades <- unique(estratos)
propiedades <- variables

plots_list <- list()
for (i in seq_along(archivos_raster)) {
  raster_file <- archivos_raster[[i]]
  etiqueta_mapa <- etiquetas_mapa[[i]]
  raster_data <- rast(raster_file)
  df <- as.data.frame(raster_data, xy = TRUE)
  colnames(df)[3] <- "Valor"
  p <- ggplot(df, aes(x = x, y = y, fill = Valor)) +
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

ancho_px <- 500 * n_columnas
alto_px  <- 375 * n_filas

for (i in seq_along(grid_plots)) {
  ggsave(sprintf("mapa_%02d.png", i), grid_plots[[i]], width=5, height=3.75, dpi=300)
}


