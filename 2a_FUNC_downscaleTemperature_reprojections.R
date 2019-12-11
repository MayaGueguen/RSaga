
### TO BE DONE ONLY ONCE !!
### INFORMATION : ERA INTERIM levels selected are going to 9000 m elevation
### (L60 model, going from 60 to 33 levels)

rm(list=ls())
library(raster)
library(rgdal)
library(ncdf4)

# path.to.data = "C:/Users/gueguen/Documents/CLIMATE_DOWNSCALING/"
# path.to.SAGA = "C:/Program Files (x86)/SAGA-GIS/"
path.to.data = "/run/user/1001/gvfs/smb-share:server=129.88.191.70,share=equipes/macroeco/GIS_DATA/CHELSA_DOWNSCALING/"
path.to.data = "/home/gueguema/Documents/CHELSA_DOWNSCALING/"
path.to.SAGA = path.to.data

zone_name.clouds = "World"
zone_name.tempCHELSA = "World"
zone_name.tempERA = "World"
zone_name.GMTED = "FID30"
proj.res.clouds = 4000
proj.res.tempCHELSA = 4000
proj.res.tempERA = 30000
proj.res.GMTED = 310
proj.name = "Mercator"
proj.value = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

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
### REPROJECT INPUT data (must be a conserving angle projection !!)
###################################################################


new.folder.name = paste0("../", zone_name.clouds, "_", proj.name, "_resolution", proj.res.clouds, "/")
if (!dir.exists(paste0(path.to.data, "CLOUDS/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "CLOUDS/RAW/", new.folder.name))
}

### Monthly cloud coverage
for (mm in list.mm)
{
  cat("\n ==> Reproject EarthEnv clouds coverage for month ", mm, "\n")
  
  input.name = paste0("CLOUDS/RAW/MODCF_monthlymean_", mm, ".tif")
  new.file.name = paste0("CLOUDS_", zone_name.clouds, "_", proj.name, "_resolution", proj.res.clouds, "_", mm, ".sgrd")
  output.name = sub(
    basename(input.name),
    paste0(new.folder.name, new.file.name),
    input.name
  )
  
  ## DO NOT WORK
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


###################################################################

new.folder.name = paste0("../", zone_name.tempCHELSA, "_", proj.name, "_resolution", proj.res.tempCHELSA, "/")
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name))
}

### CHELSA Temperature : mean, max, min : CURRENT
for (mm in list.mm)
{
  for (i in 2:3)
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
                              , " -TARGET_USER_SIZE="
                              , unique(proj.res.tempCHELSA)
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

### CHELSA Temperature : mean, max, min : PAST
new.folder.name = paste0("../", zone_name.tempCHELSA, "_", proj.name, "_resolution", proj.res.tempCHELSA, "_PAST/")
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name))
}

for(ye in past.years)
{
  for (mm in list.mm)
  {
    for (i in 2:3)
    {
      cat("\n ==> Reproject CHELSA ", c("MEAN","MAX","MIN")[i], "temperature for year ", ye, " and month ", mm, "\n")

      input.name = paste0("TEMPERATURE/RAW_TS_PAST/CHELSA_", c("temp","tmax","tmin")[i], "_", ye, "_", mm, "_V1.2.1.tif")
      new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone_name.tempCHELSA, "_"
                             , proj.name, "_resolution", proj.res.tempCHELSA, "_", ye, "_", mm, ".sgrd")
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
}

### CHELSA Temperature : mean, max, min : FUTURE
new.folder.name = paste0("../", zone_name.tempCHELSA, "_", proj.name, "_resolution", proj.res.tempCHELSA, "_FUTURE/")
if (!dir.exists(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "TEMPERATURE/RAW/", new.folder.name))
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
          new.file.name = paste0("TEMP_", c("MEAN","MAX","MIN")[i], "_", zone_name.tempCHELSA, "_", proj.name
                                 , "_resolution", proj.res.tempCHELSA, "_", sce, "_rcp", rcp, "_", mm, "_", ye, ".sgrd")
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
    }
  }
}


###################################################################

# new.folder.name.1 = paste0("../", zone_name.tempERA, "_longlat_resolution0.75/")
# if (!dir.exists(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name.1)))
# {
#   dir.create(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name.1))
# }
new.folder.name.2 = paste0("../", zone_name.tempERA, "_", proj.name,"_resolution", proj.res.tempERA, "/")
if (!dir.exists(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name.2)))
{
  dir.create(paste0(path.to.data, "LAPSE_RATE/RAW/", new.folder.name.2))
}

# ### ERA-interim temperature
# for (mm in 1:12)
# {
#   input.name.nc = paste0("LAPSE_RATE/RAW/ERAinterim_modelLevels_Temperature_2017_", mm, ".nc")
#   
#   for (lev in 1:28)
#   {
#     cat("\n ==> Extract and average ERA-interim temperature for month ", mm, " level ", lev, "\n")
#     
#     setwd(path.to.data)
#     input.ras = brick(input.name.nc, level=lev)
#     input.ras = stack(input.ras)
#     input.ras = mean(input.ras, na.rm = T)
#     new.file.name = paste0("ERAinterim_TEMP_MEAN_", zone_name.tempERA, "_longlat_resolution", unique(res(input.ras)),"_", mm,"_LEVEL", lev, ".img")
#     output.name = sub(
#       basename(input.name.nc),
#       paste0(new.folder.name.1, new.file.name),
#       input.name.nc
#     )
#     if (!file.exists(paste0(path.to.data, output.name)))
#     {
#       writeRaster(input.ras, filename = output.name)
#     }
#     
#     ### NOT WORKING : gives a 2 column band
#     # setwd(path.to.SAGA)
#     # 
#     # cat("\n ==> Reproject ERA-interim temperature for month ", mm, " level ", lev, "\n")
#     # 
#     # input.name = output.name
#     # new.file.name = paste0("ERAinterim_TEMP_MEAN_", zone_name.tempERA, "_", proj.name, "_resolution", proj.res.tempERA,"_", mm, "_LEVEL", lev, ".sgrd")
#     # output.name = sub(
#     #   basename(input.name.nc),
#     #   paste0(new.folder.name.2, new.file.name),
#     #   input.name.nc
#     # )
#     # 
#     # if (!file.exists(paste0(path.to.data, output.name)))
#     # {
#     #   system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
#     #                           , paste0("\"", proj.value, "\"")
#     #                           , " -SOURCE="
#     #                           , paste0("\"", path.to.data, input.name, "\"")
#     #                           , " -GRIDS="
#     #                           , paste0("\"", path.to.data, output.name, "\"")
#     #                           , " -RESAMPLING=3") ## B-spline interpolation
#     #   
#     #   system(system.command)
#     # }
#   }
# }

### ERA5 temperature
past.years = 1979:2019
levels.pressure = rev(c(seq(250, 750, 50), seq(775, 1000, 25)))
for (lev in 1:length(levels.pressure))
{
  input.name.nc = paste0("LAPSE_RATE/RAW/ERA5_modelLevels_Temperature_1979-2019_LEVEL", lev, "_", levels.pressure[lev], "hPa.nc")
  
  for (ye in past.years)
  {
    cat("\n ==> Extract and average ERA-interim temperature for level ", lev, " and year ", ye, "\n")
    
    setwd(path.to.data)
    input.ras = brick(input.name.nc)
    input.ras = stack(input.ras)
    input.ras = input.ras[[grep(paste0("^X", ye), names(input.ras))]]
    new.mm = sapply(names(input.ras), function(x) strsplit(x, "[.]")[[1]][2])
    # new.file.name = paste0("ERA5_TEMP_MEAN_", zone_name.tempERA, "_longlat_resolution", unique(res(input.ras))
    #                        , "_", as.numeric(new.mm) ,"_LEVEL", lev, ".img")
    # output.name = sapply(new.file.name
    #                      , function(x) sub(
    #                        basename(input.name.nc),
    #                        paste0(new.folder.name.1, x),
    #                        input.name.nc
    #                      ))
    new.file.name = paste0("ERA5_TEMP_MEAN_", zone_name.tempERA, "_", proj.name, "_resolution", proj.res.tempERA
                           , "_", as.numeric(new.mm) ,"_LEVEL", lev, ".img")
    output.name = sapply(new.file.name
                         , function(x) sub(
                           basename(input.name.nc),
                           paste0(new.folder.name.2, x),
                           input.name.nc
                         ))
    if (!file.exists(paste0(path.to.data, output.name[1])))
    {
      hop = projectRaster(input.ras, crs = CRS(proj.value))
      writeRaster(input.ras, filename = output.name, bylayer = TRUE)
    }
    
    ### NOT WORKING : gives a 2 column band
    # setwd(path.to.SAGA)
    # 
    # input.name = output.name
    # new.file.name = paste0("ERA5_TEMP_MEAN_", zone_name.tempERA, "_", proj.name, "_resolution", proj.res.tempERA
    #                        , "_", as.numeric(new.mm) ,"_LEVEL", lev, ".sgrd")
    # output.name = sapply(new.file.name
    #                      , function(x) sub(
    #                        basename(input.name.nc),
    #                        paste0(new.folder.name.2, x),
    #                        input.name.nc
    #                      ))
    # 
    # for (i in 1:length(input.name))
    # {
    #   cat("\n ==> Reproject ERA-interim temperature for month ", mm, " level ", lev, "\n")
    #   
    #   if (!file.exists(paste0(path.to.data, output.name[i])))
    #   {
    #     system.command = paste0("saga_cmd pj_proj4 3 -CRS_PROJ4="
    #                             , paste0("\"", proj.value, "\"")
    #                             , " -SOURCE="
    #                             , paste0("\"", path.to.data, input.name[i], "\"")
    #                             , " -GRIDS="
    #                             , paste0("\"", path.to.data, output.name[i], "\"")
    #                             , " -RESAMPLING=3") ## B-spline interpolation
    #     
    #     system(system.command)
    #   }
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
