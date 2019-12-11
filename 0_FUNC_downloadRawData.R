
### TO BE DONE ONLY ONCE !!

rm(list=ls())
library(foreach)
library(doParallel)
registerDoParallel(cores = 7)

path.to.data = "/home/gueguema/Documents/CHELSA_DOWNSCALING/"
setwd(path.to.data)

list.mm = c(paste0("0", 1:9), 10:12)
past.years = 1979:2013

fut.scenarios = c("CESM1-BGC", "IPSL-CM5A-MR", "MPI-ESM-LR")
fut.rcp = c("45", "85")
fut.years = c("2041-2060", "2061-2080")

fut.ts.scenarios1 = c("cmcc_cmcc-cm_ar5", "csiro-bom_access1-3_ar5", "miroc_miroc5_ar5", "nsf-doe-ncar_cesm1-bgc_ar5")
fut.ts.scenarios2 = c("CESM1-BGC", "IPSL-CM5A-MR", "MIROC5", "CESM1-BGC")
fut.ts.rcp = c("45", "85")
fut.ts.years = c("2041-2060", "2061-2080")


###################################################################
### CHELSA DATA (1979 - 2013)
###################################################################

path.CHELSA = "https://www.wsl.ch/lud/chelsa/data/"

path.CHELSA.climatologies = paste0(path.CHELSA, "climatologies/")
path.CHELSA.climatologies.prec = paste0(path.CHELSA.climatologies, "prec/")
path.CHELSA.climatologies.tmin = paste0(path.CHELSA.climatologies, "temp/integer/tmin/")
path.CHELSA.climatologies.tmax = paste0(path.CHELSA.climatologies, "temp/integer/tmax/")
path.CHELSA.climatologies.fut = paste0(path.CHELSA, "cmip5/")

path.CHELSA.timeseries = paste0(path.CHELSA, "timeseries/")
path.CHELSA.timeseries.prec = paste0(path.CHELSA.timeseries, "prec/")
path.CHELSA.timeseries.tmin = paste0(path.CHELSA.timeseries, "tmin/")
path.CHELSA.timeseries.tmax = paste0(path.CHELSA.timeseries, "tmax/")
path.CHELSA.timeseries.fut = paste0(path.CHELSA, "cmip5_ts/")

dir.create(paste0(path.to.data, "PRECIPITATION/RAW/"), recursive = TRUE)
dir.create(paste0(path.to.data, "PRECIPITATION/RAW_TS_PAST/"), recursive = TRUE)
dir.create(paste0(path.to.data, "PRECIPITATION/RAW_TS_FUTURE/"), recursive = TRUE)
dir.create(paste0(path.to.data, "TEMPERATURE/RAW/"), recursive = TRUE)
dir.create(paste0(path.to.data, "TEMPERATURE/RAW_TS_PAST/"), recursive = TRUE)
dir.create(paste0(path.to.data, "TEMPERATURE/RAW_TS_FUTURE/"), recursive = TRUE)



## CLIMATOLOGIES
foreach(mm = list.mm) %dopar%
  {
    system(paste0("wget -P ", path.to.data, "PRECIPITATION/RAW/ ", path.CHELSA.climatologies.prec, "CHELSA_prec_", mm, "_V1.2_land.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW/ ", path.CHELSA.climatologies.tmin, "CHELSA_tmin10_", mm, "_1979-2013_V1.2_land.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW/ ", path.CHELSA.climatologies.tmax, "CHELSA_tmax10_", mm, "_1979-2013_V1.2_land.tif"))
  }

## TIMESERIES
combi = expand.grid(mm = list.mm, ye = past.years)
foreach(mm = combi$mm, ye = combi$ye) %dopar%
  {
    cat("\n ==> Year ", ye, " month ",  mm)
    system(paste0("wget -P ", path.to.data, "PRECIPITATION/RAW_TS_PAST/ ", path.CHELSA.timeseries.prec, "CHELSA_prec_", ye, "_", mm, "_V1.2.1.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW_TS_PAST/ ", path.CHELSA.timeseries.tmin, "CHELSA_tmin_", ye, "_", mm, "_V1.2.1.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW_TS_PAST/ ", path.CHELSA.timeseries.tmax, "CHELSA_tmax_", ye, "_", mm, "_V1.2.1.tif"))
  }

## FUTURE CLIMATOLOGIES
combi = expand.grid(mm = list.mm, sce = fut.scenarios, rcp = fut.rcp, ye = fut.years)
foreach(mm = combi$mm, sce = combi$sce, rcp = combi$rcp, ye = combi$ye) %do%
  {
    system(paste0("wget -P ", path.to.data, "PRECIPITATION/RAW/ ", path.CHELSA.climatologies.fut, ye
                  , "/prec/CHELSA_pr_mon_", sce, "_rcp", rcp, "_r1i1p1_g025.nc_", mm, "_", ye, ".tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW/ ", path.CHELSA.climatologies.fut, ye
                  , "/prec/CHELSA_tasmin_mon_", sce, "_rcp", rcp, "_r1i1p1_g025.nc_", mm, "_", ye, "_V1.2.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW/ ", path.CHELSA.climatologies.fut, ye
                  , "/prec/CHELSA_tasmax_mon_", sce, "_rcp", rcp, "_r1i1p1_g025.nc_", mm, "_", ye, "_V1.2.tif"))
  }

###################################################################
### EarthEnv CLOUD COVER (2000 - 2014)
###################################################################

path.EarthEnv = "https://data.earthenv.org/cloud/"

dir.create(paste0(path.to.data, "CLOUDS/RAW/"), recursive = TRUE)

foreach(mm = list.mm) %dopar%
  {
    system(paste0("wget -P ", path.to.data, "CLOUDS/RAW/ ", path.EarthEnv, "MODCF_monthlymean_", mm, ".tif"))
  }


###################################################################
### ERA5 (previously ERA interim)
###################################################################

path.ERA5 = "https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels-monthly-means?tab=form"

dir.create(paste0(path.to.data, "LAPSE_RATE/RAW/"), recursive = TRUE)

## Monthly averaged reanalysis
## Temperature
## 250 to 1000 hPa (one by one)
## From 1979 to 2019
## 12 months
## 00:00

###################################################################
### GMTED2010
###################################################################

## GLOBAL
system(paste0("wget -P ", path.to.data, "DEM/RAW/ http://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/Grid_ZipFiles/mn30_grd.zip"))

## TILES
path.GMTED = "https://www.usgs.gov/land-resources/eros/coastal-changes-and-impacts/gmted2010?qt-science_support_page_related_con=0#qt-science_support_page_related_con"

## Mean 7.5 arc sec