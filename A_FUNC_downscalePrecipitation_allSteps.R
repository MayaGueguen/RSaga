
rm(list=ls())
# .libPaths("/bettik/emabio/R_PKG_NIX/")
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
zone_name = "Alps"
DEM_name = "DEM/RAW/DEM_Alps_ETRS89_resolution25.img"
DEM_ras = raster(DEM_name)
proj.res = unique(res(DEM_ras))


# -------------------------------------------------------------------- #
# SPECIFY FIXED PARAMETERS

zone_name.precip = "World" ## DO NOT CHANGE !
proj.res.precip = 850  ## DO NOT CHANGE !

proj.name = "ETRS89"
proj.value = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

list.mm = c(paste0("0", 1:9), 10:12)
past.years = 1979:2013
fut.scenarios = c("CMCC-CM", "ACCESS1-0", "MIROC5", "CESM1-BGC")
fut.rcp = c("45", "85")
fut.years = c("2041-2060", "2061-2080")
fut.ts.scenarios = c("CMCC-CM", "ACCESS1-3", "MIROC5", "CESM1-BGC")
fut.ts.rcp = c("45", "85")
fut.ts.years = seq(2010, 2100, 10)

setwd(path.to.SAGA)

### !!!!
## PROBLEM : Code seems fine, but reprojected maps are cut
## They have be done separately with ARCGIS
### !!!!

###################################################################
### REPROJECT INPUT data (must be EQUAL AREA projection !!)
### CURRENT
###################################################################


new.folder.name = paste0("../", zone_name.precip, "_", proj.name, "_resolution", proj.res.precip, "/")
if (!dir.exists(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name))
}


### Monthly precipitations
for (mm in 1:12)
{
  cat("\n ==> Reproject CHELSA precipitations for month ", mm, "\n")
  
  input.name = paste0("PRECIPITATION/RAW/CHELSA_prec_", mm, "_V1.2_land.tif")
  new.file.name = paste0("PRECIP_", zone_name.precip, "_", proj.name, "_resolution", proj.res.precip, "_", mm, ".sgrd")
  output.name = sub(
    basename(input.name),
    paste0(new.folder.name, new.file.name),
    input.name
  )
  
  if (!file.exists(paste0(path.to.data, output.name)))
  {
    # system.command = paste0("saga_cmd grid_tools 12 -INPUT="
    #                         , paste0("\"", path.to.data, input.name, "\"")
    #                         , " -OUTPUT="
    #                         , paste0("\"", path.to.data, output.name, "\"")
    #                         , " -METHOD=0 -IDENTITY=\"new_param.txt\"")
    # 
    # system(system.command)
    # 
    # system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
    #                         , paste0("\"", proj.value, "\"")
    #                         , " -SOURCE="
    #                         , paste0("\"", path.to.data, output.name, "\"")
    #                         , " -GRIDS="
    #                         , paste0("\"", path.to.data, output.name, "\"")
    #                         , " -RESAMPLING=3") ## B-spline interpolation
    # 
    # system(system.command)
  }
}

###################################################################
### REPROJECT INPUT data (must be EQUAL AREA projection !!)
### FUTURE
###################################################################


new.folder.name = paste0("../", zone_name.precip, "_", proj.name, "_resolution", proj.res.precip, "_FUTURE/")
if (!dir.exists(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name))
}

# for (sce in fut.scenarios)
# {
#   for (rcp in fut.rcp)
#   {
#     for (ye in fut.years)
#     {
#       ### Monthly precipitations
#       for (mm in 1:12)
#       {
#         cat("\n ==> Reproject CHELSA precipitations for month ", mm, "\n")
#         
#         input.name = paste0("PRECIPITATION/RAW/FUTURE/CHELSA_pr_mon_", sce, "_rcp", rcp, "_r1i1p1_g025.nc_", mm, "_", ye, ".tif")
#         new.file.name = paste0("PRECIP_", zone_name.precip, "_", proj.name, "_resolution", proj.res.precip, "_"
#                                , sce, "_rcp", rcp, "_", mm, "_", ye, ".sgrd")
#         output.name = sub(
#           basename(input.name),
#           paste0(new.folder.name, new.file.name),
#           input.name
#         )
#         
#         if (!file.exists(paste0(path.to.data, output.name)))
#         {
#           system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
#                                   , paste0("\"", proj.value, "\"")
#                                   , " -SOURCE="
#                                   , paste0("\"", path.to.data, input.name, "\"")
#                                   , " -GRIDS="
#                                   , paste0("\"", path.to.data, output.name, "\"")
#                                   , " -RESAMPLING=3") ## B-spline interpolation
#           
#           system(system.command)
#         }
#       }
#     }
#   }
# }

###################################################################
### REPROJECT INPUT data (must be EQUAL AREA projection !!)
###################################################################

### DEM
new.folder.name = paste0("../", zone_name, "_", proj.name, "_resolution", proj.res, "/")
if (!dir.exists(paste0(path.to.data, "DEM/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "DEM/RAW/", new.folder.name))
}

input.name = DEM_name
output.name = sub(basename(input.name)
                  , paste0(new.folder.name, "DEM_", zone_name, "_", proj.name, "_resolution", proj.res,".sgrd")
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


###################################################################
### CLIP INPUT data
###################################################################

DEM_ras = raster(paste0(path.to.data, sub(extension(DEM_name), ".sdat", DEM_name)))

precip.folder.name = paste0("PRECIPITATION/", zone_name.precip, "_", proj.name, "_resolution", proj.res.precip, "/")
if (!dir.exists(paste0(path.to.data, sub(zone_name.precip, zone_name, precip.folder.name))))
{
  dir.create(paste0(path.to.data, sub(zone_name.precip, zone_name, precip.folder.name)))
}

### Monthly precipitations : CURRENT
for (mm in 1:12)
{
  cat("\n ==> Clip CHELSA precipitations for month ", mm, "\n")
  
  precip.file.name = paste0("PRECIP_", zone_name.precip, "_", proj.name, "_resolution", proj.res.precip, "_", mm, ".tif")
  input.name = paste0(precip.folder.name, precip.file.name)
  output.name = paste0(sub(zone_name.precip, zone_name, precip.folder.name),
                       sub(zone_name.precip, zone_name, precip.file.name))
  
  if (!file.exists(paste0(path.to.data, output.name)))
  {
    system.command = paste0("saga_cmd grid_tools 31 -GRIDS="
                            , paste0("\"", path.to.data, input.name, "\"")
                            , " -CLIPPED="
                            , paste0("\"", path.to.data, output.name, "\"")
                            , " -EXTENT=0 -XMIN="
                            , extent(DEM_ras)[1]
                            , " -XMAX="
                            , extent(DEM_ras)[2]
                            , " -YMIN="
                            , extent(DEM_ras)[3]
                            , " -YMAX="
                            , extent(DEM_ras)[4]
                            , " -BUFFER=0.000000")
    system(system.command)
  }
}


precip.folder.name = paste0("PRECIPITATION/", zone_name.precip, "_", proj.name, "_resolution", proj.res.precip, "_FUTURE/")
if (!dir.exists(paste0(path.to.data, sub(zone_name.precip, zone_name, precip.folder.name))))
{
  dir.create(paste0(path.to.data, sub(zone_name.precip, zone_name, precip.folder.name)))
}

### Monthly precipitations : FUTURE
for (sce in fut.scenarios)
{
  for (rcp in fut.rcp)
  {
    for (ye in fut.years)
    {
      for (mm in 1:12)
      {
        cat("\n ==> Clip CHELSA precipitations for ", sce, rcp, ye, " and month ", mm, "\n")
        
        precip.file.name = paste0("PRECIP_", zone_name.precip, "_", proj.name, "_resolution", proj.res.precip, "_"
                               , sce, "_rcp", rcp, "_", mm, "_", ye, ".tif")
        input.name = paste0(precip.folder.name, precip.file.name)
        output.name = paste0(sub(zone_name.precip, zone_name, precip.folder.name),
                             sub(zone_name.precip, zone_name, precip.file.name))
        
        if (!file.exists(paste0(path.to.data, output.name)))
        {
          system.command = paste0("saga_cmd grid_tools 31 -GRIDS="
                                  , paste0("\"", path.to.data, input.name, "\"")
                                  , " -CLIPPED="
                                  , paste0("\"", path.to.data, output.name, "\"")
                                  , " -EXTENT=0 -XMIN="
                                  , extent(DEM_ras)[1]
                                  , " -XMAX="
                                  , extent(DEM_ras)[2]
                                  , " -YMIN="
                                  , extent(DEM_ras)[3]
                                  , " -YMAX="
                                  , extent(DEM_ras)[4]
                                  , " -BUFFER=0.000000")
          system(system.command)
        }
      }
    }
  }
}

###################################################################
### GEOGRAPHICALLY weighted regression
###################################################################

precip.folder.name = paste0("PRECIPITATION/", zone_name, "_", proj.name, "_resolution", proj.res.precip, "/")
new.folder.name = paste0("PRECIPITATION/", zone_name, "_", proj.name, "_resolution", proj.res, "/")
if (!dir.exists(paste0(path.to.data, new.folder.name)))
{
  dir.create(paste0(path.to.data, new.folder.name))
}

### Monthly precipitations in function of DEM : CURRENT
for (mm in 1:12)
{
  cat("\n ==> GWR of CHELSA precipitations in function of DEM for month ", mm, "\n")
  
  predic.name = DEM_name
  
  precip.file.name = paste0("PRECIP_", zone_name, "_", proj.name, "_resolution", proj.res.precip, "_", mm, ".sgrd")
  input.name = paste0(precip.folder.name, precip.file.name)
  output.name = sub(proj.res.precip, proj.res, precip.file.name)
  output.name.1 = paste0(new.folder.name, sub(extension(output.name), "_regression.sgrd", output.name))
  output.name.2 = paste0(new.folder.name, sub(extension(output.name), "_regression_rescorr.sgrd", output.name))
  
  if (!file.exists(paste0(path.to.data, output.name.1)))
  {
    system.command = paste0("saga_cmd statistics_regression 14 -PREDICTORS="
                            , paste0("\"", path.to.data, predic.name, "\"")
                            , " -REGRESSION="
                            , paste0("\"", path.to.data, output.name.1, "\"")
                            , " -REG_RESCORR="
                            , paste0("\"", path.to.data, output.name.2, "\"")
                            , " -DEPENDENT="
                            , paste0("\"", path.to.data, input.name, "\""))
    
    system(system.command) 
  }
}

### Not good enough ? Maybe try :
### Clouds in function of DEM
### Monthly precipitations in function of DEM and Clouds

precip.folder.name = paste0("PRECIPITATION/", zone_name, "_", proj.name, "_resolution", proj.res.precip, "_FUTURE/")
new.folder.name = paste0("PRECIPITATION/", zone_name, "_", proj.name, "_resolution", proj.res, "_FUTURE/")
if (!dir.exists(paste0(path.to.data, new.folder.name)))
{
  dir.create(paste0(path.to.data, new.folder.name))
}

### Monthly precipitations in function of DEM : FUTURE
for (sce in fut.scenarios)
{
  for (rcp in fut.rcp)
  {
    for (ye in fut.years)
    {
      for (mm in 1:12)
      {
        cat("\n ==> GWR of CHELSA precipitations in function of DEM for ", sce, rcp, ye, " and month ", mm, "\n")
        
        predic.name = DEM_name
        
        precip.file.name = paste0("PRECIP_", zone_name, "_", proj.name, "_resolution", proj.res.precip, "_"
                                  , sce, "_rcp", rcp, "_", mm, "_", ye, ".sgrd")
        input.name = paste0(precip.folder.name, precip.file.name)
        output.name = sub(proj.res.precip, proj.res, precip.file.name)
        output.name.1 = paste0(new.folder.name, sub(extension(output.name), "_regression.sgrd", output.name))
        output.name.2 = paste0(new.folder.name, sub(extension(output.name), "_regression_rescorr.sgrd", output.name))
        
        if (!file.exists(paste0(path.to.data, output.name.1)))
        {
          system.command = paste0("saga_cmd statistics_regression 14 -PREDICTORS="
                                  , paste0("\"", path.to.data, predic.name, "\"")
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
