
rm(list=ls())

# machine = "leca"
machine = "luke"
# machine = "froggy"


## LOAD R packages
if (machine == "luke")
{
  .libPaths("/bettik/emabio/R_PKG_NIX/")
} else if (machine == "froggy")
{
  .libPaths("/home/mayagueguen/R_PKG_NIX/")
}
library(raster)
library(rgdal)


## DEFINE paths to data and to SAGA executable file
if (machine == "leca")
{
  # path.to.data = "C:/Users/gueguen/Documents/CLIMATE_DOWNSCALING/"
  # path.to.SAGA = "C:/Program Files (x86)/SAGA-GIS/"
  path.to.data = "/media/gueguen/equipes/macroeco/GIS_DATA/CHELSA_DOWNSCALING/"
} else if (machine == "luke")
{
  path.to.data = "/bettik/mayagueguen/CHELSA_DOWNSCALING/"
} else if (machine == "froggy")
{
  path.to.data = "/scratch/mayagueguen/CHELSA_DOWNSCALING/"
}
path.to.SAGA = path.to.data


setwd(path.to.data)
zone_name.clouds = "World" ## DO NOT CHANGE !
zone_name.tempCHELSA = "World" ## DO NOT CHANGE !
zone_name.tempERA = "World" ## DO NOT CHANGE !
zone_name.GMTED = "FID30"
proj.res.clouds = 6000 ## DO NOT CHANGE !
proj.res.tempCHELSA = 4000 ## DO NOT CHANGE !
proj.res.tempERA = "100000" ## DO NOT CHANGE !
proj.res.GMTED = 310 ## DO NOT CHANGE !

# zone_name = "Bauges"
# DEM_name = "DEM/RAW/DEM_Bauges.img"
# zone_name = "Lautaret"
# DEM_name = "DEM/RAW/DEM_Lautaret.img"
zone_name = "Alps"
DEM_name = "DEM/RAW/DEM_Alps_ETRS89_resolution25.img"
# zone_name = "FrenchAlps"
# DEM_name = "DEM/RAW/DEM_FrenchAlps.img"

DEM_ras = raster(DEM_name)
proj.res = unique(res(DEM_ras))[1]
proj.name = "Mercator"
proj.value = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
setwd(path.to.SAGA)

# system("export SAGA_MLB=/home/gueguen/Documents/_SOFTWARES/SAGA_LST_Dirk/Release/")


###################################################################
### CREATE PATHWAYS and FILENAMES
###################################################################

zone.file.name = paste0(zone_name, "_", proj.name, "_resolution", proj.res)
zone.folder.name = paste0(zone.file.name, "/")

clouds.file.name = paste0(zone_name, "_", proj.name, "_resolution", proj.res.clouds)

###################################################################
### REPROJECT INPUT data (must be a conserving angle projection !!)
###################################################################

### DEM
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


###################################################################
### SOLAR RADIATION
###################################################################

new.folder.name = paste0("SOLAR_RADIATION/", zone_name, "_", proj.name,"_resolution", proj.res, "/")
if (!dir.exists(paste0(path.to.data, new.folder.name)))
{
  dir.create(paste0(path.to.data, new.folder.name))
}

args = commandArgs(trailingOnly=T)
mm = as.numeric(args[1])
dd = as.numeric(args[2])
for (VAR in c(DEM_name, input.name.DEM.flat))
{
  
  input.name.dem = VAR
  input.name.svf = sub(extension(input.name.dem), "_SVF.sgrd", input.name.dem)
  if (VAR == input.name.DEM.flat) input.name.svf = VAR
  
  cat("\n ==> Calculate solar radiation for month ", mm, " from day ", dd, "\n")
  mm = as.numeric(mm)
  mm.end = mm + 1
  if (mm == 12) mm.end = 1
  # mm.end = mm
  # if (dd == 1) dd.end = 10
  # if (dd == 10) dd.end = 20
  # if (dd == 20)
  # {
  #   dd.end = 1
  #   mm.end = mm + 1
  # }
  # if (dd == 20 && mm == 12) mm.end = 1
  yy.end = ifelse(mm == 12, 2019, 2018)
  dd.end = 1
  cat("STARTING POINT : ",ISOdate(2018, mm, dd), "\n")
  cat("ENDING POINT : ",ISOdate(yy.end, mm.end, dd.end), "\n")
  
  nb.days = nrow(as.data.frame(seq.POSIXt(from = ISOdate(2018, mm, dd),
                                          to = ISOdate(yy.end, mm.end, dd.end),
                                          by = "day"))) - 1
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
    if (mm <= 9 && nchar(mm) == 1) mm = paste0("0",mm)
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
                            , " -LOCATION=1 -PERIOD=2 -DAY=2018-", mm, "-1 -DAY_STOP=2018-", mm, "-", nb.days#dd + nb.days -1
                            , " -DAYS_STEP=1")
    
    system(system.command)
  }
}
