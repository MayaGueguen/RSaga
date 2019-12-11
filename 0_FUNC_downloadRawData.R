
### TO BE DONE ONLY ONCE !!

rm(list=ls())
library(foreach)
library(doParallel)
registerDoParallel(cores = 7)

path.to.data = "/home/gueguema/Documents/CHELSA_DOWNSCALING/"
setwd(path.to.data)

list.mm = c(paste0("0", 1:9), 10:12)
past.years = 1979:2013

fut.scenarios = c("ACCESS1-0", "BNU-ESM", "CCSM4"
                  , "CESM1-BGC", "CESM1-CAM5", "CMCC-CMS", "CMCC-CM", "CNRM-CM5"
                  , "CSIRO-Mk3-6-0", "CSIRO-Mk3L-1-2", "CanESM2"
                  , "EC-EARTH", "FGOALS-g2", "FIO-ESM", "GFDL-CM3", "GFDL-ESM2G", "GFDL-ESM2M"
                  , "GISS-E2-H-CC", "GISS-E2-H", "GISS-E2-R-CC", "GISS-E2-R"
                  , "HadGEM2-AO", "HadGEM2-CC", "HadGEM2-ES"
                  , "IPSL-CM5A-LR", "IPSL-CM5A-MR", "MIROC-ESM-CHEM", "MIROC-ESM", "MIROC5"
                  , "MPI-ESM-LR", "MPI-ESM-MR", "MRI-CGCM3", "MRI-ESM1"
                  , "NorESM1-M", "bcc-csm1-1", "inmcm4", "")
fut.scenarios = c("CESM1-BGC", "IPSL-CM5A-MR", "MPI-ESM-LR")
fut.scenarios = c("CMCC-CM", "ACCESS1-0", "MIROC5", "CESM1-BGC")
fut.rcp = c("45", "85")
fut.years = c("2041-2060", "2061-2080")

fut.ts.scenarios1 = c("cmcc_cmcc-cm_ar5", "csiro-bom_access1-3_ar5", "miroc_miroc5_ar5", "nsf-doe-ncar_cesm1-bgc_ar5")
fut.ts.scenarios2 = c("CMCC-CM", "ACCESS1-3", "MIROC5", "CESM1-BGC")
fut.ts.rcp = c("45", "85")
fut.ts.years = c(2010:2019, seq(2020, 2100, 10))
fut.ts.years = seq(2010, 2100, 10)


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
dir.create(paste0(path.to.data, "PRECIPITATION/RAW_FUTURE/"), recursive = TRUE)
dir.create(paste0(path.to.data, "PRECIPITATION/RAW_TS_PAST/"), recursive = TRUE)
dir.create(paste0(path.to.data, "PRECIPITATION/RAW_TS_FUTURE/"), recursive = TRUE)
dir.create(paste0(path.to.data, "TEMPERATURE/RAW/"), recursive = TRUE)
dir.create(paste0(path.to.data, "TEMPERATURE/RAW_FUTURE/"), recursive = TRUE)
dir.create(paste0(path.to.data, "TEMPERATURE/RAW_TS_PAST/"), recursive = TRUE)
dir.create(paste0(path.to.data, "TEMPERATURE/RAW_TS_FUTURE/"), recursive = TRUE)



## CLIMATOLOGIES
foreach(mm = list.mm) %dopar%
  {
    system(paste0("wget -P ", path.to.data, "PRECIPITATION/RAW/ ", path.CHELSA.climatologies.prec, "CHELSA_prec_", mm, "_V1.2_land.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW/ ", path.CHELSA.climatologies.tmin, "CHELSA_tmin10_", mm, "_1979-2013_V1.2_land.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW/ ", path.CHELSA.climatologies.tmax, "CHELSA_tmax10_", mm, "_1979-2013_V1.2_land.tif"))
  }

## PAST TIMESERIES
combi = expand.grid(mm = list.mm, ye = past.years)
foreach(mm = combi$mm, ye = combi$ye) %dopar%
  {
    cat("\n ==> Year ", ye, " month ",  mm)
    system(paste0("wget -P ", path.to.data, "PRECIPITATION/RAW_TS_PAST/ ", path.CHELSA.timeseries.prec, "CHELSA_prec_", ye, "_", mm, "_V1.2.1.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW_TS_PAST/ ", path.CHELSA.timeseries.tmin, "CHELSA_tmin_", ye, "_", mm, "_V1.2.1.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW_TS_PAST/ ", path.CHELSA.timeseries.tmax, "CHELSA_tmax_", ye, "_", mm, "_V1.2.1.tif"))
  }

## FUTURE CLIMATOLOGIES
combi = expand.grid(mm = as.numeric(list.mm), sce = fut.scenarios, rcp = fut.rcp, ye = fut.years)
foreach(mm = combi$mm, sce = combi$sce, rcp = combi$rcp, ye = combi$ye) %do%
  {
    system(paste0("wget -P ", path.to.data, "PRECIPITATION/RAW_FUTURE/ ", path.CHELSA.climatologies.fut, ye
                  , "/prec/CHELSA_pr_mon_", sce, "_rcp", rcp, "_r1i1p1_g025.nc_", mm, "_", ye, ".tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW_FUTURE/ ", path.CHELSA.climatologies.fut, ye
                  , "/tmin/CHELSA_tasmin_mon_", sce, "_rcp", rcp, "_r1i1p1_g025.nc_", mm, "_", ye, "_V1.2.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW_FUTURE/ ", path.CHELSA.climatologies.fut, ye
                  , "/tmax/CHELSA_tasmax_mon_", sce, "_rcp", rcp, "_r1i1p1_g025.nc_", mm, "_", ye, "_V1.2.tif"))
  }

## FUTURE TIMESERIES
combi = expand.grid(mm = list.mm, i.sce = length(fut.ts.scenarios1), rcp = fut.ts.rcp, ye = fut.ts.years)
foreach(mm = combi$mm, i.sce = combi$i.sce, rcp = combi$rcp, ye = combi$ye) %do%
  {
    cat("\n ==> Scenario ", fut.ts.scenarios2[i.sce], " year ", ye, " month ",  mm, " rcp ", rcp)
    
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW_TS_FUTURE/ ", path.CHELSA.timeseries.fut
                  , "/", fut.ts.scenarios1[i.sce], "/CHELSA_", fut.ts.scenarios2[i.sce]
                  , "_rcp", rcp, "_", ye, "_", mm, "_tmin.tif"))
    system(paste0("wget -P ", path.to.data, "TEMPERATURE/RAW_TS_FUTURE/ ", path.CHELSA.timeseries.fut
                  , "/", fut.ts.scenarios1[i.sce], "/CHELSA_", fut.ts.scenarios2[i.sce]
                  , "_rcp", rcp, "_", ye, "_", mm, "_tmax.tif"))
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