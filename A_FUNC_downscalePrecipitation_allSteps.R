
rm(list=ls())
library(raster)
library(rgdal)

path.to.data = "/home/gueguema/Documents/CHELSA_DOWNSCALING/"
path.to.SAGA = path.to.data
setwd(path.to.data)

# -------------------------------------------------------------------- #
# SPECIFY STUDIED AREA

# zone_name = "Bauges"
# DEM_name = "DEM/RAW/DEM_Bauges.img"
# zone_name = "Lautaret"
# DEM_name = "DEM/RAW/DEM_Lautaret.img"
# zone_name = "Alps"
# DEM_name = "DEM/RAW/DEM_Alps_ETRS89_resolution25.img"
zone_name = "Greenland"
DEM_name = "DEM/RAW/DEM_Greenland_ETRS89_resolution90.img"

# -------------------------------------------------------------------- #
# SPECIFY FIXED PARAMETERS

zone_name.precip = "World" ## DO NOT CHANGE !
proj.res.precipCHELSA = 1000  ## DO NOT CHANGE !
proj.name = "ETRS89"
proj.value = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "
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

###############################################################################
### DEM
## REPROJECT INPUT data (must be an EQUAL AREA projection !!)
###############################################################################

DEM_ras = raster(DEM_name)
proj.res = unique(res(DEM_ras))[1]

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
}

DEM_ras = raster(sub(".sgrd$", ".sdat", DEM_name))
tmp = projectExtent(DEM_ras, crs = proj.longlat)
proj.extent = extent(tmp)


###############################################################################
### PRECIPITATION (CHELSA)
## CLIP & REPROJECT INPUT data (must be an EQUAL AREA projection !!)
## DOWNSCALE (geographically weighted regression)
###############################################################################

### CHELSA Precipitation : CURRENT --------------------------------------------
new.folder.name1 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res.precipCHELSA, "/")
if (!dir.exists(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name1)))
{
  dir.create(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name1))
}
new.folder.name2 = paste0("../", zone.folder.name)
if (!dir.exists(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name2))
}

for (mm in list.mm)
{
  cat("\n ==> Reproject CHELSA precipitation for month ", mm, "\n")
  
  input.name = paste0("PRECIPITATION/RAW/CHELSA_prec_", mm, "_V1.2_land.tif")
  new.file.name = paste0("PREC_", zone_name, "_", proj.name, "_resolution", proj.res.precipCHELSA, "_", mm, ".sgrd")
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
                  , param.res = proj.res.precipCHELSA)
  }
  
  cat("\n ==> GWR of CHELSA precipitation in function of DEM for month ", mm, "\n")
  
  input.name = output.name
  new.file.name = paste0("PREC_", zone.file.name, "_", mm, ".sgrd")
  output.name = paste0("PRECIPITATION/RAW/", new.folder.name2, new.file.name)
  output.name.1 = sub(extension(output.name), "_regression.sgrd", output.name)
  output.name.2 = sub(extension(output.name), "_regression_rescorr.sgrd", output.name)
  
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

### CHELSA Precipitation : PAST TIMESERIES ------------------------------------
new.folder.name1 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res.precipCHELSA, "_TS_PAST/")
if (!dir.exists(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name1)))
{
  dir.create(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name1))
}
new.folder.name2 = paste0("../", sub("/$", "_TS_PAST/", zone.folder.name))
if (!dir.exists(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name2))
}

for(ye in past.years)
{
  for (mm in list.mm)
  {
    cat("\n ==> Reproject CHELSA precipitation for year ", ye, " and month ", mm, "\n")
    
    input.name = paste0("PRECIPITATION/RAW_TS_PAST/CHELSA_prec_", ye, "_", mm, "_V1.2.1.tif")
    new.file.name = paste0("PREC_", zone_name, "_", proj.name, "_resolution"
                           , proj.res.precipCHELSA, "_", ye, "_", mm, ".sgrd")
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
                    , param.res = proj.res.precipCHELSA)
    }
    
    cat("\n ==> GWR of CHELSA precipitation in function of DEM for year ", ye, " and month ", mm, "\n")
    
    input.name = output.name
    new.file.name = paste0("PREC_", zone.file.name, "_", ye, "_", mm, ".sgrd")
    output.name = paste0("PRECIPITATION/RAW/", new.folder.name2, new.file.name)
    output.name.1 = sub(extension(output.name), "_regression.sgrd", output.name)
    output.name.2 = sub(extension(output.name), "_regression_rescorr.sgrd", output.name)
    
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
}

### CHELSA Precipitation : FUTURE ---------------------------------------------
new.folder.name1 = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res.precipCHELSA, "_FUTURE/")
if (!dir.exists(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name1)))
{
  dir.create(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name1))
}
new.folder.name2 = paste0("../", sub("/$", "_FUTURE/", zone.folder.name))
if (!dir.exists(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name2)))
{
  dir.create(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name2))
}

for (sce in fut.scenarios)
{
  for (rcp in fut.rcp)
  {
    for (ye in fut.years)
    {
      for (mm in 1:12)
      {
        cat("\n ==> Reproject CHELSA precipitation for ", sce, rcp, ye, " and month ", mm, "\n")
        
        input.name = paste0("PRECIPITATION/RAW_FUTURE/CHELSA_pr_mon_"
                            , sce, "_rcp", rcp, "_r1i1p1_g025.nc_", mm, "_", ye, ".tif")
        new.file.name = paste0("PREC_", zone_name, "_", proj.name, "_resolution", proj.res.precipCHELSA
                               , "_", sce, "_rcp", rcp, "_", mm, "_", ye, ".sgrd")
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
                        , param.res = proj.res.precipCHELSA)
        }
        
        cat("\n ==> GWR of CHELSA precipitation in function of DEM for for ", sce, rcp, ye, " and month ", mm, "\n")
        
        input.name = output.name
        new.file.name = paste0("PREC_", zone.file.name, "_", sce, "_rcp", rcp, "_", mm, "_", ye, ".sgrd")
        output.name = paste0("PRECIPITATION/RAW/", new.folder.name2, new.file.name)
        output.name.1 = sub(extension(output.name), "_regression.sgrd", output.name)
        output.name.2 = sub(extension(output.name), "_regression_rescorr.sgrd", output.name)
        
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
    }
  }
}

### Not good enough ? Maybe try :
### Clouds in function of DEM
### Monthly precipitations in function of DEM and Clouds
