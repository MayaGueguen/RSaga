
### TO BE DONE ONLY ONCE !!
### INFORMATION : ERA INTERIM levels selected are going to 9000 m elevation
### (L60 model, going from 60 to 33 levels)

rm(list=ls())
library(raster)
library(rgdal)
library(ncdf4)

# path.to.data = "C:/Users/gueguen/Documents/CLIMATE_DOWNSCALING/"
# path.to.SAGA = "C:/Program Files (x86)/SAGA-GIS/"
# path.to.data = "/run/user/1001/gvfs/smb-share:server=129.88.191.70,share=equipes/macroeco/GIS_DATA/CHELSA_DOWNSCALING/"
path.to.data = "/home/gueguema/Documents/CHELSA_DOWNSCALING/"
path.to.SAGA = path.to.data

zone_name.tempERA = "World"
proj.res.clouds = 1000 ## 6000
proj.res.tempCHELSA = 1000 ## 4000
proj.res.tempERA = 30000
proj.res.GMTED = 30 ## 310
proj.name = "Mercator"
proj.value = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
proj.longlat = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

list.mm = c(paste0("0", 1:9), 10:12)
past.years = 1979:2013

fut.scenarios = c("CMCC-CM", "ACCESS1-0", "MIROC5", "CESM1-BGC")
fut.rcp = c("45", "85")
fut.years = c("2041-2060", "2061-2080")

fut.ts.scenarios = c("CMCC-CM", "ACCESS1-3", "MIROC5", "CESM1-BGC")
fut.ts.rcp = c("45", "85")
fut.ts.years = seq(2010, 2100, 10)

setwd(path.to.SAGA)

###################################################################
### FUNCTION : CLIP & REPROJECT

clipReproject = function(param.input, param.output, param.extent, param.proj, param.res)
{
  system.command = paste0("saga_cmd grid_tools 31 -GRIDS="
                          , paste0("\"", path.to.data, param.input, "\"")
                          , " -CLIPPED="
                          , paste0("\"", path.to.data, param.output, "\"")
                          , " -EXTENT=0 -XMIN="
                          , param.extent[1]
                          , " -XMAX="
                          , param.extent[2]
                          , " -YMIN="
                          , param.extent[3]
                          , " -YMAX="
                          , param.extent[4])
  system(system.command)
  
  system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
                          , paste0("\"", param.proj, "\"")
                          , " -SOURCE="
                          , paste0("\"", path.to.data, param.output, "\"")
                          , " -GRIDS="
                          , paste0("\"", path.to.data, param.output, "\"")
                          , " -TARGET_USER_SIZE="
                          , unique(param.res)
                          , " -RESAMPLING=3") ## B-spline interpolation
  system(system.command)
}

###################################################################

zone_name = "Greenland"
DEM_name = "DEM/RAW/DEM_Greenland_ETRS89_resolution90.img"

DEM_ras = raster(DEM_name)
proj.res = unique(res(DEM_ras))[1]

# tmp = projectExtent(DEM_ras, crs = proj.longlat)
# proj.extent = extent(tmp)

zone.file.name = paste0(zone_name, "_", proj.name, "_resolution", proj.res)
zone.folder.name = paste0(zone.file.name, "/")

## Reproject DEM file
new.folder.name = paste0("../", zone.folder.name)
if (!dir.exists(paste0(path.to.data, "DEM/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "DEM/RAW/", new.folder.name))
}

input.name = DEM_name
output.name = sub(basename(input.name)
                  , paste0(new.folder.name, "DEM_", zone.file.name, ".sgrd")
                  , input.name)
DEM_name = output.name

if (!file.exists(paste0(path.to.data, output.name)))
{
  cat("\n ==> Reproject DEM into ", proj.name, " projection and .sgrd file \n")
  
  system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
                          , paste0("\"", proj.value, "\"")
                          , " -SOURCE="
                          , paste0("\"", path.to.data, input.name, "\"")
                          , " -GRIDS="
                          , paste0("\"", path.to.data, output.name, "\"")
                          , " -TARGET_USER_SIZE="
                          , unique(proj.res)
                          , " -RESAMPLING=3") ## B-spline interpolation
  system(system.command)
  
  cat("\n ==> Check if extent is too large \n")
  
  ras = raster(sub(".sgrd$", ".sdat", output.name))
  ind.na = rowColFromCell(ras, which(is.na(ras[])))
  tab.row = table(ind.na[, 1])
  tab.col = table(ind.na[, 2])
  ind.row = which(tab.row == ncol(ras))
  ind.col = which(tab.col == nrow(ras))
  if (length(ind.row) > 0 || length(ind.col) > 0)
  {
    cat("\n ==> Create new DEM raster with adequate extent \n")
    ind.row.toKeep = 1:nrow(ras)
    if (length(ind.row) > 0) ind.row.toKeep = ind.row.toKeep[-ind.row]
    ind.col.toKeep = 1:ncol(ras)
    if (length(ind.col) > 0) ind.col.toKeep = ind.col.toKeep[-ind.col]
    new_ras = crop(ras, extent(ras
                               , min(ind.row.toKeep), max(ind.row.toKeep)
                               , min(ind.col.toKeep), max(ind.col.toKeep)))
    new_extent = extent(new_ras)
    rm(ras)
    rm(new_ras)
    
    system.command = paste0("saga_cmd grid_tools 31 -GRIDS="
                            , paste0("\"", path.to.data, output.name, "\"")
                            , " -CLIPPED="
                            , paste0("\"", path.to.data, output.name, "\"")
                            , " -EXTENT=0 -XMIN="
                            , new_extent[1]
                            , " -XMAX="
                            , new_extent[2]
                            , " -YMIN="
                            , new_extent[3]
                            , " -YMAX="
                            , new_extent[4])
    system(system.command)
  }
}

DEM_ras = raster(sub(".sgrd$", ".sdat", DEM_name))
tmp = projectExtent(DEM_ras, crs = proj.longlat)
proj.extent = extent(tmp)


###############################################################################
## LEAF AREA INDEX for LST calculation
input.name.lai = sub(basename(DEM_name), sub("DEM_", "LAI_0.01_", basename(DEM_name)), DEM_name)
if (!file.exists(paste0(path.to.data, input.name.lai)))
{
  system.command = paste0("saga_cmd grid_calculus 1 -GRIDS="
                          , paste0("\"", path.to.data, DEM_name, "\"")
                          , " -XGRIDS=NULL -RESAMPLING=3 -RESULT="
                          , paste0("\"", path.to.data, input.name.lai, "\"")
                          , " -FORMULA=0.01 -NAME="
                          , sub(".sgrd", "", sub("DEM_", "LAI_0.01_", basename(DEM_name)))
                          , " -TYPE=7")
  system(system.command) 
}

## FLAT DEM (1)
input.name.DEM.flat = sub(basename(DEM_name), sub("DEM_", "DEM_FLAT_", basename(DEM_name)), DEM_name)
if (!file.exists(paste0(path.to.data, input.name.DEM.flat)))
{
  system.command = paste0("saga_cmd grid_calculus 1 -GRIDS="
                          , paste0("\"", path.to.data, DEM_name, "\"")
                          , " -XGRIDS=NULL -RESAMPLING=3 -RESULT="
                          , paste0("\"", path.to.data, input.name.DEM.flat, "\"")
                          , " -FORMULA=1 -NAME="
                          , sub(".sgrd", "", sub("DEM_", "DEM_FLAT_", basename(DEM_name)))
                          , " -TYPE=7")
  system(system.command) 
}

###############################################################################
### CLOUDS (EarthEnv)
## CLIP & REPROJECT INPUT data (must be a conserving angle projection !!)
## DOWNSCALE (geographically weighted regression)
###############################################################################

new.folder.name1 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res.clouds, "/")
if (!dir.exists(paste0(path.to.data, "CLOUDS/RAW/", new.folder.name1)))
{
  dir.create(paste0(path.to.data, "CLOUDS/RAW/", new.folder.name1))
}
new.folder.name2 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res, "/")
if (!dir.exists(paste0(path.to.data, "CLOUDS/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "CLOUDS/RAW/", new.folder.name2))
}

### Monthly cloud coverage
for (mm in list.mm)
{
  cat("\n ==> Reproject EarthEnv clouds coverage for month ", mm, "\n")
  
  input.name = paste0("CLOUDS/RAW/MODCF_monthlymean_", mm, ".tif")
  new.file.name = paste0("CLOUDS_", zone_name, "_", proj.name, "_resolution", proj.res.clouds, "_", mm, ".sgrd")
  output.name = sub(
    basename(input.name),
    paste0(new.folder.name1, new.file.name),
    input.name
  )
  
  if (!file.exists(paste0(path.to.data, output.name)))
  {
    clipReproject(param.input = input.name
                  , param.output = output.name
                  , param.extent = proj.extent
                  , param.proj = proj.value
                  , param.res = proj.res.clouds)
  }
  
  cat("\n ==> GWR of EarthEnv clouds coverage in function of DEM for month ", mm, "\n")
  
  input.name = output.name
  output.name = sub(proj.res.clouds, proj.res, basename(input.name))
  output.name.1 = paste0("CLOUDS/RAW/", new.folder.name2, sub(extension(output.name), "_regression.sgrd", output.name))
  output.name.2 = paste0("CLOUDS/RAW/", new.folder.name2, sub(extension(output.name), "_regression_rescorr.sgrd", output.name))
  
  if (!file.exists(paste0(path.to.data, output.name.1)))
  {
    system.command = paste0("saga_cmd statistics_regression 14 -PREDICTORS="
                            , paste0("\"", path.to.data, DEM_name, "\"")
                            , " -REGRESSION="
                            , paste0("\"", path.to.data, output.name.1, "\"")
                            , " -REG_RESCORR="
                            , paste0("\"", path.to.data, output.name.2, "\"")
                            , " -DEPENDENT="
                            , paste0("\"", path.to.data, input.name, "\""))
    
    system(system.command) 
  }
}


###############################################################################
### TEMPERATURE (CHELSA)
## CLIP & REPROJECT INPUT data (must be a conserving angle projection !!)
## DOWNSCALE (geographically weighted regression)
###############################################################################

### CHELSA Temperature : mean, max, min : CURRENT -----------------------------
new.folder.name1 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res.tempCHELSA, "/")
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name1)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name1))
}
new.folder.name2 = paste0("../", zone.folder.name)
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name2))
}

for (mm in list.mm)
{
  for (i in 2:3)
  {
    cat("\n ==> Reproject CHELSA ", c("MEAN","MAX","MIN")[i], "temperature for month ", mm, "\n")
    
    input.name = paste0("TEMPERATURE/RAW/CHELSA_", c("temp","tmax","tmin")[i], "10_", mm, "_1979-2013_V1.2_land.tif")
    new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone_name, "_", proj.name, "_resolution", proj.res.tempCHELSA, "_", mm, ".sgrd")
    output.name = sub(
      basename(input.name),
      paste0(new.folder.name1, new.file.name),
      input.name
    )
    
    if (!file.exists(paste0(path.to.data, output.name)))
    {
      clipReproject(param.input = input.name
                    , param.output = output.name
                    , param.extent = proj.extent
                    , param.proj = proj.value
                    , param.res = proj.res.tempCHELSA)
      
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
    
    cat("\n ==> Downscale CHELSA ", c("MEAN","MAX","MIN")[i], "temperature for month ", mm, "\n")
    
    input.name = output.name
    new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone.file.name, "_", mm, ".sgrd")
    output.name = paste0("TEMPERATURE/RAW/", new.folder.name2, new.file.name)
    
    if (!file.exists(paste0(path.to.data, output.name)))
    {
      system.command = paste0("saga_cmd grid_tools 0 -INPUT="
                              , paste0("\"", path.to.data, input.name, "\"")
                              , " -OUTPUT="
                              , paste0("\"", path.to.data, output.name, "\"")
                              , " -SCALE_DOWN=3"
                              , " -TARGET_DEFINITION=1"
                              , " -TARGET_TEMPLATE="
                              , paste0("\"", path.to.data, DEM_name, "\""))
      
      system(system.command)
    }
  }
}

### CHELSA Temperature : mean, max, min : PAST TIMESERIES ---------------------
new.folder.name1 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res.tempCHELSA, "_TS_PAST/")
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name1)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name1))
}
new.folder.name2 = paste0("../", sub("/$", "_TS_PAST/", zone.folder.name))
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name2))
}

for(ye in past.years)
{
  for (mm in list.mm)
  {
    for (i in 2:3)
    {
      cat("\n ==> Reproject CHELSA ", c("MEAN","MAX","MIN")[i], "temperature for year ", ye, " and month ", mm, "\n")

      input.name = paste0("TEMPERATURE/RAW_TS_PAST/CHELSA_", c("temp","tmax","tmin")[i], "_", ye, "_", mm, "_V1.2.1.tif")
      new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone_name, "_"
                             , proj.name, "_resolution", proj.res.tempCHELSA, "_", ye, "_", mm, ".sgrd")
      output.name = sub(
        basename(input.name),
        paste0(new.folder.name1, new.file.name),
        input.name
      )

      if (!file.exists(paste0(path.to.data, output.name)))
      {
        clipReproject(param.input = input.name
                      , param.output = output.name
                      , param.extent = proj.extent
                      , param.proj = proj.value
                      , param.res = proj.res.tempCHELSA)

        system.command = paste0("saga_cmd grid_calculus 1 -GRIDS="
                                , paste0("\"", path.to.data, output.name, "\"")
                                , " -XGRIDS=NULL -RESAMPLING=3 -RESULT="
                                , paste0("\"", path.to.data, output.name, "\"")
                                , " -FORMULA=\"g1 / 10 - 273.15\""
                                , " -NAME="
                                , paste0("\"", sub(extension(output.name), "", basename(output.name)), "\"")
                                , " -TYPE=7")
        system(system.command)
      }
      
      cat("\n ==> Downscale CHELSA ", c("MEAN","MAX","MIN")[i], "temperature for year ", ye, " and month ", mm, "\n")
      
      input.name = output.name
      new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone.file.name, "_", ye, "_", mm, ".sgrd")
      output.name = paste0("TEMPERATURE/RAW/", new.folder.name2, new.file.name)
      
      if (!file.exists(paste0(path.to.data, output.name)))
      {
        system.command = paste0("saga_cmd grid_tools 0 -INPUT="
                                , paste0("\"", path.to.data, input.name, "\"")
                                , " -OUTPUT="
                                , paste0("\"", path.to.data, output.name, "\"")
                                , " -SCALE_DOWN=3"
                                , " -TARGET_DEFINITION=1"
                                , " -TARGET_TEMPLATE="
                                , paste0("\"", path.to.data, DEM_name, "\""))
        
        system(system.command)
      }
    }
  }
}

### CHELSA Temperature : mean, max, min : FUTURE ------------------------------
new.folder.name1 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res.tempCHELSA, "_FUTURE/")
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name1)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name1))
}
new.folder.name2 = paste0("../", sub("/$", "_FUTURE/", zone.folder.name))
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name2))
}

for (sce in fut.scenarios)
{
  for (rcp in fut.rcp)
  {
    for (ye in fut.years)
    {
      for (mm in 1:12)
      {
        for (i in 2:3)
        {
          cat("\n ==> Reproject CHELSA ", c("MEAN","MAX","MIN")[i], "temperature for ", sce, rcp, ye, " and month ", mm, "\n")
          
          input.name = paste0("TEMPERATURE/RAW_FUTURE/CHELSA_", c("tas","tasmax","tasmin")[i]
                              , "_mon_", sce, "_rcp", rcp, "_r1i1p1_g025.nc_", mm, "_", ye, "_V1.2.tif")
          new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone_name, "_", proj.name
                                 , "_resolution", proj.res.tempCHELSA, "_", sce, "_rcp", rcp, "_", mm, "_", ye, ".sgrd")
          output.name = sub(
            basename(input.name),
            paste0(new.folder.name1, new.file.name),
            input.name
          )
          
          if (!file.exists(paste0(path.to.data, output.name)))
          {
            clipReproject(param.input = input.name
                          , param.output = output.name
                          , param.extent = proj.extent
                          , param.proj = proj.value
                          , param.res = proj.res.tempCHELSA)
            
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
          
          cat("\n ==> Downscale CHELSA ", c("MEAN","MAX","MIN")[i], "temperature for ", sce, rcp, ye, " and month ", mm, "\n")
          
          input.name = output.name
          new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone.file.name, "_", sce, "_rcp", rcp, "_", mm, "_", ye, ".sgrd")
          output.name = paste0("TEMPERATURE/RAW/", new.folder.name2, new.file.name)
          
          if (!file.exists(paste0(path.to.data, output.name)))
          {
            system.command = paste0("saga_cmd grid_tools 0 -INPUT="
                                    , paste0("\"", path.to.data, input.name, "\"")
                                    , " -OUTPUT="
                                    , paste0("\"", path.to.data, output.name, "\"")
                                    , " -SCALE_DOWN=3"
                                    , " -TARGET_DEFINITION=1"
                                    , " -TARGET_TEMPLATE="
                                    , paste0("\"", path.to.data, DEM_name, "\""))
            
            system(system.command)
          }
        }
      }
    }
  }
}

### CHELSA Temperature : mean, max, min : FUTURE TIMESERIES -------------------
new.folder.name1 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res.tempCHELSA, "_TS_FUTURE/")
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name1)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name1))
}
new.folder.name2 = paste0("../", sub("/$", "_TS_FUTURE/", zone.folder.name))
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name2))
}

for (sce in fut.ts.scenarios)
{
  for (rcp in fut.ts.rcp)
  {
    for (ye in fut.ts.years)
    {
      for (mm in list.mm)
      {
        for (i in 2:3)
        {
          cat("\n ==> Reproject CHELSA ", c("MEAN","MAX","MIN")[i], "temperature for ", sce, rcp, ye, " and month ", mm, "\n")
          
          input.name = paste0("TEMPERATURE/RAW_TS_FUTURE/CHELSA_", sce, "_rcp", rcp, "_", ye, "_", mm, "_", c("tas","tmax","tmin")[i], ".tif")
          new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone_name, "_", proj.name
                                 , "_resolution", proj.res.tempCHELSA, "_", sce, "_rcp", rcp, "_", mm, "_", ye, ".sgrd")
          output.name = sub(
            basename(input.name),
            paste0(new.folder.name1, new.file.name),
            input.name
          )
          
          if (!file.exists(paste0(path.to.data, output.name)))
          {
            clipReproject(param.input = input.name
                          , param.output = output.name
                          , param.extent = proj.extent
                          , param.proj = proj.value
                          , param.res = proj.res.tempCHELSA)
            
            system.command = paste0("saga_cmd grid_calculus 1 -GRIDS="
                                    , paste0("\"", path.to.data, output.name, "\"")
                                    , " -XGRIDS=NULL -RESAMPLING=3 -RESULT="
                                    , paste0("\"", path.to.data, output.name, "\"")
                                    , " -FORMULA=\"g1 / 10 - 273.15\""
                                    , " -NAME="
                                    , paste0("\"", sub(extension(output.name), "", basename(output.name)), "\"")
                                    , " -TYPE=7")
            system(system.command)
          }
          
          cat("\n ==> Downscale CHELSA ", c("MEAN","MAX","MIN")[i], "temperature for ", sce, rcp, ye, " and month ", mm, "\n")
          
          input.name = output.name
          new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone.file.name, "_", sce, "_rcp", rcp, "_", mm, "_", ye, ".sgrd")
          output.name = paste0("TEMPERATURE/RAW/", new.folder.name2, new.file.name)
          
          if (!file.exists(paste0(path.to.data, output.name)))
          {
            system.command = paste0("saga_cmd grid_tools 0 -INPUT="
                                    , paste0("\"", path.to.data, input.name, "\"")
                                    , " -OUTPUT="
                                    , paste0("\"", path.to.data, output.name, "\"")
                                    , " -SCALE_DOWN=3"
                                    , " -TARGET_DEFINITION=1"
                                    , " -TARGET_TEMPLATE="
                                    , paste0("\"", path.to.data, DEM_name, "\""))
            
            system(system.command)
          }
        }
      }
    }
  }
}


###############################################################################
### GMTED2010
## CLIP & REPROJECT INPUT data (must be a conserving angle projection !!)
## DOWNSCALE (geographically weighted regression)
###############################################################################

new.folder.name1 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res.GMTED, "/")
if (!dir.exists(paste0(path.to.data, "DEM/RAW/", new.folder.name1)))
{
  dir.create(paste0(path.to.data, "DEM/RAW/", new.folder.name1))
}
new.folder.name2 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res, "/")
if (!dir.exists(paste0(path.to.data, "DEM/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "DEM/RAW/", new.folder.name2))
}


cat("\n ==> Reproject GMTED2010 DEM \n")

input.name = "DEM/RAW/mn30_grd.img"
new.file.name = paste0("DEM_REF_", zone_name, "_", proj.name, "_resolution", proj.res.GMTED, ".sgrd")
output.name = sub(
  basename(input.name),
  paste0(new.folder.name1, new.file.name),
  input.name
)

if (!file.exists(paste0(path.to.data, output.name)))
{
  clipReproject(param.input = input.name
                , param.output = output.name
                , param.extent = proj.extent
                , param.proj = proj.value
                , param.res = proj.res.tempCHELSA)
}

cat("\n ==> Downscale GMTED2010 DEM \n")

input.name = output.name
new.file.name = paste0("DEM_REF_", zone_name, "_", proj.name, "_resolution", proj.res, ".sgrd")
output.name = paste0("DEM/RAW/", new.folder.name2, new.file.name)

if (!file.exists(paste0(path.to.data, output.name)))
{
  system.command = paste0("saga_cmd grid_tools 0 -INPUT="
                          , paste0("\"", path.to.data, input.name, "\"")
                          , " -OUTPUT="
                          , paste0("\"", path.to.data, output.name, "\"")
                          , " -SCALE_DOWN=3"
                          , " -TARGET_DEFINITION=1"
                          , " -TARGET_TEMPLATE="
                          , paste0("\"", path.to.data, DEM_name, "\""))
  
  system(system.command)
}


###############################################################################
### LAPSE-RATE (ERA5)
## CLIP & REPROJECT INPUT data (must be a conserving angle projection !!)
## DOWNSCALE (geographically weighted regression)
###############################################################################

new.folder.name1 = paste0("../", zone_name.tempERA, "_longlat_resolution0.25/")
if (!dir.exists(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name1)))
{
  dir.create(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name1))
}
new.folder.name2 = paste0("../", zone_name, "_", proj.name,"_resolution", proj.res.tempERA, "/")
if (!dir.exists(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name2))
}


### ERA5 temperature
ERA5.years = 1979:2019
levels.pressure = rev(c(seq(250, 750, 50), seq(775, 1000, 25)))
for (lev in 1:length(levels.pressure))
{
  input.name.nc = paste0("LAPSE_RATE/RAW/ERA5_modelLevels_Temperature_1979-2019_LEVEL", lev, "_", levels.pressure[lev], "hPa.nc")
  
  for (ye in ERA5.years)
  {
    cat("\n ==> Extract and average ERA-interim temperature for level ", lev, " and year ", ye, "\n")
    
    input.ras = brick(input.name.nc)
    input.ras = stack(input.ras)
    input.ras = input.ras[[grep(paste0("^X", ye), names(input.ras))]]
    origin(input.ras) = c(0,0)
    
    new.mm = sapply(names(input.ras), function(x) strsplit(x, "[.]")[[1]][2])
    new.file.name = paste0("ERA5_TEMP_MEAN_", zone_name.tempERA, "_longlat_resolution", unique(res(input.ras))
                           , "_", as.numeric(new.mm) , "_", ye, "_LEVEL", lev, ".img")
    output.name = sapply(new.file.name
                         , function(x) sub(
                           basename(input.name.nc),
                           paste0(new.folder.name1, x),
                           input.name.nc
                         ))
    
    if (!file.exists(paste0(path.to.data, output.name[1])))
    {
      input.pts = as.data.frame(rasterToPoints(input.ras))
      input.pts$x[which(input.pts$x > 180)] =  (input.pts$x[which(input.pts$x > 180)] - 360)
      new_ras = rasterFromXYZ(input.pts, res = unique(res(input.ras)), crs= CRS(proj.longlat))
      
      writeRaster(new_ras, filename = output.name, bylayer = TRUE)
    }
    
    cat("\n ==> Clip and reproject ERA-interim temperature for level ", lev, " and year ", ye, "\n")
    
    input.name = output.name
    new.file.name = paste0("ERA5_TEMP_MEAN_", zone_name, "_", proj.name, "_resolution", proj.res.tempERA
                           , "_", as.numeric(new.mm) , "_", ye, "_LEVEL", lev, ".sgrd")
    output.name = sapply(new.file.name
                         , function(x) sub(
                           basename(input.name.nc),
                           paste0(new.folder.name2, x),
                           input.name.nc
                         ))

    for (i in 1:length(input.name))
    {
      if (!file.exists(paste0(path.to.data, output.name[i])))
      {
        clipReproject(param.input = input.name[i]
                      , param.output = output.name[i]
                      , param.extent = proj.extent
                      , param.proj = proj.value
                      , param.res = proj.res.tempERA)
      }
    }
  }
}

### Monthly temperature with model levels
new.folder.name1 = paste0("../", zone_name, "_", proj.name,"_resolution", proj.res.tempERA, "/")
new.folder.name2 = paste0("../", zone_name, "_", proj.name,"_resolution", proj.res, "/")
if (!dir.exists(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name2))
}

for (ye in ERA5.years)
{
  for (mm in 1:12)
  {
    cat("\n ==> Compute lapse-rate (temperature ~ elevation level) for year ", ye, " and month ", mm, "\n")
    
    input.name = paste0("ERA5_TEMP_MEAN_", zone_name, "_", proj.name, "_resolution", proj.res.tempERA, "_", mm, "_", ye)
    input.name = paste0(input.name, "_LEVEL", 1:length(levels.pressure))
    input.name = paste0(input.name, ".sgrd")
    input.name = paste0("LAPSE_RATE/RAW/", new.folder.name1, input.name)
    
    new.file.name = paste0("LAPSE_RATE_", zone_name, "_", proj.name, "_resolution", proj.res.tempERA, "_", mm, "_", ye, ".sgrd")
    output.name = paste0("LAPSE_RATE/RAW/", new.folder.name1, new.file.name)
    output.name.1 = sub(extension(output.name), "_coeff1.sgrd", output.name)
    output.name.2 = sub(extension(output.name), "_coeff2.sgrd", output.name)
    
    L60_model_levels = "ERA5_L60_model_levels.txt"
    
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
    
    input.name = output.name.2
    new.file.name = paste0("LAPSE_RATE_", zone_name, "_", proj.name, "_resolution", proj.res, "_", mm, "_", ye, ".sgrd")
    new.file.name = sub(extension(new.file.name), "_coeff2.sgrd", new.file.name)
    output.name = paste0("LAPSE_RATE/RAW/", new.folder.name2, new.file.name)

    if (!file.exists(paste0(path.to.data, output.name)))
    {
      system.command = paste0("saga_cmd grid_tools 0 -INPUT="
                              , paste0("\"", path.to.data, input.name, "\"")
                              , " -OUTPUT="
                              , paste0("\"", path.to.data, output.name, "\"")
                              , " -SCALE_DOWN=3"
                              , " -TARGET_DEFINITION=1"
                              , " -TARGET_TEMPLATE="
                              , paste0("\"", path.to.data, DEM_name, "\""))

      system(system.command)
    }
  }
}


###################################################################
### SKY VIEW FACTOR
###################################################################

### DEM or DEM FLAT
for (VAR in c(DEM_name, input.name.DEM.flat))
{
  input.name = VAR
  output.name.vis = sub(extension(input.name), "_VISIBLE.sgrd", input.name)
  output.name.svf = sub(extension(input.name), "_SVF.sgrd", input.name)
  
  if (!file.exists(paste0(path.to.data, output.name.svf)))
  {
    cat("\n ==> Calculating sky view factor \n")
    
    system.command = paste0("saga_cmd ta_lighting 3 -DEM="
                            , paste0("\"", path.to.data, input.name, "\"")
                            , " -VISIBLE="
                            , paste0("\"", path.to.data, output.name.vis, "\"")
                            , " -SVF="
                            , paste0("\"", path.to.data, output.name.svf, "\""))
    
    system(system.command)
  }
}


###################################################################
### SOLAR RADIATION
###################################################################

new.folder.name = paste0("SOLAR_RADIATION/", zone_name, "_", proj.name,"_resolution", proj.res, "/")
if (!dir.exists(paste0(path.to.data, new.folder.name)))
{
  dir.create(paste0(path.to.data, new.folder.name))
}

# input.name.DEM.flat = sub(basename(DEM_name), sub("DEM_", "DEM_FLAT_", basename(DEM_name)), DEM_name)
for (VAR in c(DEM_name, input.name.DEM.flat))
{
  input.name.dem = VAR
  input.name.svf = sub(extension(input.name.dem), "_SVF.sgrd", input.name.dem)
  
  # for (ye in c(1979:2020, seq(2030, 2100, 10)))
  ye = 2018
  {
    for (mm in 1:12)
    {
      dd = 1
      {
        cat("\n ==> Calculate solar radiation for month ", mm, " from day ", dd, "\n")
        mm = as.numeric(mm)
        mm.end = mm + 1
        if (mm == 12) mm.end = 1
        yy.end = ifelse(mm == 12, ye + 1, ye)
        dd.end = 1
        nb.days = nrow(as.data.frame(seq.POSIXt(from = ISOdate(ye, mm, dd),
                                                to = ISOdate(yy.end, mm.end, dd.end),
                                                by = "day"))) - 1
        
        cat("STARTING POINT : ", ISOdate(ye, mm, dd), "\n")
        cat("ENDING POINT : ", ISOdate(yy.end, mm.end, dd.end), "\n")
        cat("NB DAYS : ", nb.days, "\n")
        
        output.name.direct = paste0(new.folder.name, "DirectRad_", zone_name, "_", proj.name,"_resolution", proj.res, "_", mm, ".sgrd")
        output.name.diffus = paste0(new.folder.name, "DiffuseRad_", zone_name, "_", proj.name,"_resolution", proj.res, "_", mm, ".sgrd")
        output.name.total = paste0(new.folder.name, "TotalRad_", zone_name, "_", proj.name,"_resolution", proj.res, "_", mm, ".sgrd")
        if (length(grep("FLAT", VAR)) > 0){
          output.name.direct = sub(extension(output.name.direct), "_DEM_FLAT.sgrd", output.name.direct)
          output.name.diffus = sub(extension(output.name.diffus), "_DEM_FLAT.sgrd", output.name.diffus)
          output.name.total = sub(extension(output.name.total), "_DEM_FLAT.sgrd", output.name.total)
        }
        
        if (!file.exists(paste0(path.to.data, output.name.total)))
        {
          if (mm <= 9 && nchar(mm) == 1) mm = paste0("0", mm)
          system.command = paste0("saga_cmd ta_lighting 2 -GRD_DEM="
                                  , paste0("\"", path.to.data, input.name.dem, "\"")
                                  , " -GRD_SVF="
                                  , paste0("\"", path.to.data, input.name.svf, "\"")
                                  , " -GRD_DIRECT="
                                  , paste0("\"", path.to.data, output.name.direct, "\"")
                                  , " -GRD_DIFFUS="
                                  , paste0("\"", path.to.data, output.name.diffus, "\"")
                                  , " -GRD_TOTAL="
                                  , paste0("\"", path.to.data, output.name.total, "\"")
                                  , " -LOCATION=1 -PERIOD=2 -DAY=", ye, "-", mm, "-1 -DAY_STOP=", ye, "-", mm, "-", nb.days
                                  , " -DAYS_STEP=1")
          
          
          system(system.command)
        }
      }
    }
  }
}