library(terra)

# Load the raster
raster_file <- "F:\\ArcPro\\Stinknet_RF_Updtetd\\ONPI RF\\ENV_Rasts\\DEM_repro.tif"
r <- rast(raster_file)

# Load the shapefile
shapefile <- "F:\\ArcPro\\Stinknet_RF_Updtetd\\ONPI RF\\AZBound.shp"
shp <- vect(shapefile)

# Clip the raster using the shapefile
clipped_raster <- crop(r, shp)
clipped_raster <- mask(clipped_raster, shp)

# Save the clipped raster
writeRaster(clipped_raster, "F:\\ArcPro\\Stinknet_RF_Updtetd\\ONPI RF\\ENV_Cliped\\DEM.tif",  overwrite = TRUE)
