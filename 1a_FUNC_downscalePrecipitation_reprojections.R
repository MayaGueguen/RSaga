
### TO BE DONE ONLY ONCE !!

rm(list=ls())
library(raster)
library(rgdal)

# path.to.data = "C:/Users/gueguen/Documents/CLIMATE_DOWNSCALING/"
# path.to.SAGA = "C:/Program Files (x86)/SAGA-GIS/"
path.to.data = "/media/gueguen/equipes/emabio/GIS_DATA/CHELSA_DOWNSCALING/"
path.to.SAGA = path.to.data
path.to.data = "/bettik/mayagueguen/CHELSA_DOWNSCALING/"
path.to.SAGA = path.to.data

zone_name.precip = "World"
proj.res.precip = 1200
proj.name = "ETRS89"
proj.value = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

fut.scenarios = c("CESM1-BGC", "IPSL-CM5A-MR", "MPI-ESM-LR")
fut.rcp = c("45", "85")
fut.years = c("2061-2080")

setwd(path.to.SAGA)


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
### REPROJECT INPUT data (must be EQUAL AREA projection !!)
### FUTURE
###################################################################


new.folder.name = paste0("../", zone_name.precip, "_", proj.name, "_resolution", proj.res.precip, "_FUTURE/")
if (!dir.exists(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name)))
{
  dir.create(paste0(path.to.data, "PRECIPITATION/RAW/", new.folder.name))
}

for (sce in fut.scenarios)
{
  for (rcp in fut.rcp)
  {
    for (ye in fut.years)
    {
      ### Monthly precipitations
      for (mm in 1:12)
      {
        cat("\n ==> Reproject CHELSA precipitations for month ", mm, "\n")
        
        input.name = paste0("PRECIPITATION/RAW/FUTURE/CHELSA_pr_mon_", sce, "_rcp", rcp, "_r1i1p1_g025.nc_", mm, "_", ye, ".tif")
        new.file.name = paste0("PRECIP_", zone_name.precip, "_", proj.name, "_resolution", proj.res.precip, "_"
                               , sce, "_rcp", rcp, "_", mm, "_", ye, ".sgrd")
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
    }
  }
}