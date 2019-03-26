
### TO BE DONE ONLY ONCE !!
### INFORMATION : ERA INTERIM levels selected are going to 9000 m elevation
### (L60 model, going from 60 to 33 levels)

rm(list=ls())
library(raster)
library(rgdal)
library(ncdf4)

# path.to.data = "C:/Users/gueguen/Documents/CLIMATE_DOWNSCALING/"
# path.to.SAGA = "C:/Program Files (x86)/SAGA-GIS/"
path.to.data = "/media/gueguen/equipes/macroeco/GIS_DATA/CHELSA_DOWNSCALING/"
path.to.SAGA = path.to.data

zone_name.clouds = "World"
zone_name.tempCHELSA = "World"
zone_name.tempERA = "World"
zone_name.GMTED = "FID30"
proj.res.clouds = 6000
proj.res.tempCHELSA = 4000
proj.res.tempERA = "100000"
proj.res.GMTED = "310"
proj.name = "Mercator"
proj.value = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

setwd(path.to.SAGA)


###################################################################
### REPROJECT INPUT data (must be a conserving angle projection !!)
###################################################################


new.folder.name = paste0("../", zone_name.clouds, "_", proj.name, "_resolution", proj.res.clouds, "/")
if (!dir.exists(paste0(path.to.data, "CLOUDS/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "CLOUDS/RAW/", new.folder.name))
}

### Monthly cloud coverage
for (mm in 1:12)
{
  cat("\n ==> Reproject EarthEnv clouds coverage for month ", mm, "\n")
  
  input.name = paste0("CLOUDS/RAW/MODCF_monthlymean_", mm, ".tif")
  new.file.name = paste0("CLOUDS_", zone_name.clouds, "_", proj.name, "_resolution", proj.res.clouds, "_", mm, ".sgrd")
  output.name = sub(
    basename(input.name),
    paste0(new.folder.name, new.file.name),
    input.name
  )
  
  if (!file.exists(paste0(path.to.data, output.name)))
  {
    system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
                            , paste0("\"", proj.value, "\"")
                            , " -SOURCE="
                            , paste0("\"", path.to.data, input.name, "\"")
                            , " -GRIDS="
                            , paste0("\"", path.to.data, output.name, "\"")
                            , " -RESAMPLING=3") ## B-spline interpolation
    
    system(system.command)
  }
}


###################################################################

new.folder.name = paste0("../", zone_name.tempCHELSA, "_", proj.name, "_resolution", proj.res.tempCHELSA, "/")
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name))
}

### CHELSA Temperature : mean, max, min : CURRENT
for (mm in 1:12)
{
  for (i in 1:3)
  {
    cat("\n ==> Reproject CHELSA ", c("MEAN","MAX","MIN")[i], "temperature for month ", mm, "\n")
    
    input.name = paste0("TEMPERATURE/RAW/CHELSA_", c("temp","tmax","tmin")[i], "10_", mm, "_1979-2013_V1.2_land.tif")
    new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone_name.tempCHELSA, "_", proj.name, "_resolution", proj.res.tempCHELSA, "_", mm, ".sgrd")
    output.name = sub(
      basename(input.name),
      paste0(new.folder.name, new.file.name),
      input.name
    )
    
    if (!file.exists(paste0(path.to.data, output.name)))
    {
      system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
                              , paste0("\"", proj.value, "\"")
                              , " -SOURCE="
                              , paste0("\"", path.to.data, input.name, "\"")
                              , " -GRIDS="
                              , paste0("\"", path.to.data, output.name, "\"")
                              , " -RESAMPLING=3") ## B-spline interpolation
      
      system(system.command)
      
      system.command = paste0("saga_cmd grid_calculus 1 -GRIDS="
                              , paste0("\"", path.to.data, output.name, "\"")
                              , " -XGRIDS=NULL -RESAMPLING=3 -RESULT="
                              , paste0("\"", path.to.data, output.name, "\"")
                              , " -FORMULA=\"g1 / 10\""
                              , " -NAME="
                              , paste0("\"", sub(extension(output.name), "", basename(output.name)), "\"")
                              , " -TYPE=7")
      system(system.command) 
    }
  }
}


###################################################################

new.folder.name.1 = paste0("../", zone_name.tempERA, "_longlat_resolution0.75/")
if (!dir.exists(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name.1)))
{
  dir.create(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name.1))
}
new.folder.name.2 = paste0("../", zone_name.tempERA, "_", proj.name,"_resolution", proj.res.tempERA, "/")
if (!dir.exists(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name.2)))
{
  dir.create(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name.2))
}

### ERA-interim temperature
for (mm in 1:12)
{
  input.name.nc = paste0("LAPSE_RATE/RAW/ERAinterim_modelLevels_Temperature_2017_", mm, ".nc")
  
  for (lev in 1:28)
  {
    cat("\n ==> Extract and average ERA-interim temperature for month ", mm, " level ", lev, "\n")
    
    setwd(path.to.data)
    input.ras = brick(input.name.nc, level=lev)
    input.ras = stack(input.ras)
    input.ras = mean(input.ras, na.rm = T)
    new.file.name = paste0("ERAinterim_TEMP_MEAN_", zone_name.tempERA, "_longlat_resolution", unique(res(input.ras)),"_", mm,"_LEVEL", lev, ".img")
    output.name = sub(
      basename(input.name.nc),
      paste0(new.folder.name.1, new.file.name),
      input.name.nc
    )
    if (!file.exists(paste0(path.to.data, output.name)))
    {
      writeRaster(input.ras, filename = output.name)
    }
    
    ### NOT WORKING : gives a 2 column band
    # setwd(path.to.SAGA)
    # 
    # cat("\n ==> Reproject ERA-interim temperature for month ", mm, " level ", lev, "\n")
    # 
    # input.name = output.name
    # new.file.name = paste0("ERAinterim_TEMP_MEAN_", zone_name.tempERA, "_", proj.name, "_resolution", proj.res.tempERA,"_", mm, "_LEVEL", lev, ".sgrd")
    # output.name = sub(
    #   basename(input.name.nc),
    #   paste0(new.folder.name.2, new.file.name),
    #   input.name.nc
    # )
    # 
    # if (!file.exists(paste0(path.to.data, output.name)))
    # {
    #   system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
    #                           , paste0("\"", proj.value, "\"")
    #                           , " -SOURCE="
    #                           , paste0("\"", path.to.data, input.name, "\"")
    #                           , " -GRIDS="
    #                           , paste0("\"", path.to.data, output.name, "\"")
    #                           , " -RESAMPLING=3") ## B-spline interpolation
    #   
    #   system(system.command)
    # }
  }
}

###################################################################

new.folder.name = paste0("../", zone_name.GMTED, "_", proj.name, "_resolution", proj.res.GMTED, "/")
if (!dir.exists(paste0(path.to.data, "DEM/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "DEM/RAW/", new.folder.name))
}

### GMTED2010 DEM
cat("\n ==> Reproject GMTED2010 DEM \n")

input.name = paste0("DEM/RAW/30N000E_20101117_gmted_mea075.tif")
new.file.name = paste0("DEM_REF_", zone_name.GMTED, "_", proj.name, "_resolution", proj.res.GMTED, ".sgrd")
output.name = sub(
  basename(input.name),
  paste0(new.folder.name, new.file.name),
  input.name
)

if (!file.exists(paste0(path.to.data, output.name)))
{
  system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
                          , paste0("\"", proj.value, "\"")
                          , " -SOURCE="
                          , paste0("\"", path.to.data, input.name, "\"")
                          , " -GRIDS="
                          , paste0("\"", path.to.data, output.name, "\"")
                          , " -RESAMPLING=3") ## B-spline interpolation
  
  system(system.command)
}


###################################################################
### LAPSE RATE
###################################################################

new.folder.name = paste0("../", zone_name.tempERA, "_", proj.name,"_resolution", proj.res.tempERA, "/")

### Monthly temperature with model levels
for (mm in 1:12)
{
  cat("\n ==> Compute lapse-rate (temperature ~ elevation level) for month ", mm, "\n")
  
  input.name = paste0("ERAinterim_TEMP_MEAN_", zone_name.tempERA, "_", proj.name, "_resolution", proj.res.tempERA, "_", mm)
  input.name = paste0(input.name, "_LEVEL", 1:28)
  input.name = paste0(input.name, ".img")
  input.name = paste0("LAPSE_RATE/RAW/", new.folder.name, input.name)
  
  new.file.name = paste0("LAPSE_RATE_", zone_name.tempERA, "_", proj.name, "_resolution", proj.res.tempERA, "_", mm, ".sgrd")
  output.name = paste0("LAPSE_RATE/RAW/", new.folder.name, new.file.name)
  output.name.1 = sub(extension(output.name), "_coeff1.sgrd", output.name)
  output.name.2 = sub(extension(output.name), "_coeff2.sgrd", output.name)
  
  L60_model_levels = "ERA-interim_L60_model_levels.txt"
  
  if (!file.exists(paste0(path.to.data, output.name.1)))
  {
    system.command = paste0("saga_cmd statistics_regression 9 -Y_GRIDS="
                            , paste0("\"", paste0(path.to.data, input.name, collapse = ";"), "\"")
                            , " -COEFF="
                            , paste0("\"", paste0(path.to.data, c(output.name.1, output.name.2), collapse = ";"), "\"")
                            , " -ORDER=1 -XSOURCE=1 -X_TABLE="
                            , paste0("\"", path.to.data, L60_model_levels, "\""))
    
    system(system.command)
  }
}
