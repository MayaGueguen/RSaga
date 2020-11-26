
rm(list=ls())
library(raster)
library(rgdal)

path.to.SAGA = "C:/Program Files (x86)/SAGA-GIS/"
path.to.data = "Z:/macroeco/GIS_DATA/CHELSA_DOWNSCALING/LAND_SURFACE_TEMPERATURE/"
setwd(path.to.SAGA)

proj.value = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "


###################################################################################################

list.fi = list.files(path = paste0(path.to.data, "Alps_Mercator_resolution25/")
                     , pattern = ".sgrd$"
                     , full.names = TRUE)

for(proj.res in c("100", "250"))
{
  new.dir = paste0(path.to.data, "Alps_Mercator_resolution", proj.res, "/")
  if (!dir.exists(new.dir)) dir.create(new.dir)
  
  for(fi in list.fi)
  {
    param.input = fi
    param.output = gsub("Alps_Mercator_resolution25", paste0("Alps_Mercator_resolution", proj.res), fi)
    
    system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
                            , paste0("\"", proj.value, "\"")
                            , " -SOURCE="
                            , paste0("\"", param.input, "\"")
                            , " -GRIDS="
                            , paste0("\"", param.output, "\"")
                            , " -TARGET_USER_SIZE="
                            , proj.res
                            , " -RESAMPLING=3") ## B-spline interpolation
    system(system.command)
  }
}