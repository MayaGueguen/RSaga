
rm(list=ls())

machine = "leca"
# machine = "luke"
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
  # path.to.data = "/media/gueguen/equipes/macroeco/GIS_DATA/CHELSA_DOWNSCALING/"
  # path.to.data = "K:/LECA/macroeco/GIS_DATA/CHELSA_DOWNSCALING/"
  # path.to.data = "/home/gueguen/Bureau/CHELSA_DOWNSCALING/"
  path.to.data = "/run/user/30241/gvfs/smb-share:server=129.88.191.70,share=equipes/macroeco/GIS_DATA/CHELSA_DOWNSCALING/"
} else if (machine == "luke")
{
  path.to.data = "/bettik/mayagueguen/CHELSA_DOWNSCALING/"
} else if (machine == "froggy")
{
  path.to.data = "/scratch/mayagueguen/CHELSA_DOWNSCALING/"
}
path.to.SAGA = path.to.data
# path.to.SAGA = "C:/Users/renaud/Downloads/saga-7.2.0_x64/"


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

###################################################################
### LAND SURFACE TEMPERATURE computation
###################################################################

### WITH MAP VALUES

new.folder.name = paste0(zone_name, "_", proj.name,"_resolution", proj.res, "/")
lst.folder.name = paste0("LAND_SURFACE_TEMPERATURE/", zone_name, "_", proj.name,"_resolution", proj.res, "/")

library(foreach)
library(viridis)
# library(doParallel)
# registerDoParallel(cores = 8)

combi = expand.grid(mm = 1:12, i = 1:3)
cairo_pdf("COMPARISON_Chelsa_LST.pdf", width = 8.5, height = 7.7)
COR = foreach (mm = combi$mm, i = combi$i, .combine = "rbind") %do%
{
    cat("\n ==> Calculate ", c("MIN","MEAN","MAX")[i], " Land Surface Temperature for month ", mm, "\n")
    
    ## Temperature REF (CHELSA)
    input.name.temp = paste0("TEMP_", c("MIN","MEAN","MAX")[i], "_", zone.file.name, "_", mm, ".sdat")
    input.name.temp = paste0("TEMPERATURE/", new.folder.name, input.name.temp)
    
    input.name.temp.ori = paste0("CHELSA_", c("tmin", "temp", "tmax")[i], "10_", mm, "_1979-2013_V1.2_land.tif")
    input.name.temp.ori = paste0("TEMPERATURE/RAW/", input.name.temp.ori)
    
    ## Land Surface Temperature
    output.name = paste0("LST_", c("MIN_","MEAN_","MAX_")[i], zone.file.name, "_", mm, ".sdat")
    output.name = paste0(lst.folder.name, output.name)
    
    ##
    TEMP.ori = raster(input.name.temp.ori) ## ~ 1000m - world - long/lat
    TEMP = raster(input.name.temp) ## 25m - alps - mercator
    LST = raster(output.name) ## 25m - alps - mercator

    ##
    no.pix = 500
    
    ##
    TEMP.tmp1 = crop(TEMP, extent(TEMP, 9961, 11000, 34961, 36000)) ## 25m - crop - mercator
    new_extent = projectExtent(TEMP.tmp1, crs = CRS(projection(TEMP.ori)))
    TEMP.ori = crop(TEMP.ori, extent(new_extent)) ## ~ 1000m - crop - long/lat
    TEMP.tmp0 = projectRaster(TEMP.ori, res = 1000, crs = CRS(projection(TEMP.tmp1))) ## 1000m - crop - mercator
    TEMP.tmp0 = TEMP.tmp0 / 10
    TEMP.tmp0 = crop(TEMP.tmp0, extent(TEMP.tmp0, 4, 29, 4, 29))
    
    LST.tmp1 = crop(LST, extent(LST, 9961, 11000, 34961, 36000)) ## 25m - crop - mercator
    LST.tmp2 = aggregate(LST.tmp1, fact = 40, fun = mean) ## 1000m - crop - mercator
    
    ##
    PIX = sample(1:ncell(TEMP.tmp1), no.pix)
    XY = xyFromCell(TEMP.tmp1, PIX)
    
    PIX_TEMP_ori = TEMP.tmp0[cellFromXY(TEMP.tmp0, XY)]
    PIX_TEMP = TEMP.tmp1[PIX]
    PIX_LST = LST.tmp1[PIX]
    PIX_LST_agg = LST.tmp2[cellFromXY(LST.tmp2, XY)]
    
    ##
    PIX_agg = sample(1:ncell(LST.tmp2), no.pix)

    LST.tmp3 = LST.tmp2 ## 1000m - crop - mercator
    LST.tmp3[] = NA
    LST.tmp3[PIX_agg] = PIX_agg
    LST.tmp3 = disaggregate(LST.tmp3, fact = 40) ## 25m - crop - mercator
    
    zon.cv = zonal(LST.tmp1, LST.tmp3, fun = 'cv')
    zon.sd = zonal(LST.tmp1, LST.tmp3, fun = 'sd')
    
    ##
    MIN = min(c(min(TEMP.tmp0[], na.rm = TRUE)
                , min(TEMP.tmp1[], na.rm = TRUE)
                , min(LST.tmp1[], na.rm = TRUE)
                , min(LST.tmp2[], na.rm = TRUE)))
    MAX = max(c(max(TEMP.tmp0[], na.rm = TRUE)
                , max(TEMP.tmp1[], na.rm = TRUE)
                , max(LST.tmp1[], na.rm = TRUE)
                , max(LST.tmp2[], na.rm = TRUE)))
    
    par(mfrow = c(3,2), mar = c(2, 2, 4, 0))
    plot(TEMP.tmp0
         , zlim = c(MIN, MAX)
         , box = FALSE
         , axes = FALSE
         , col = viridis(100)
         , main = "1000 m"
         , legend.width = 2
         , legend.mar = 10
    )
    plot(TEMP.tmp1
         , zlim = c(MIN, MAX)
         , box = FALSE
         , axes = FALSE
         , col = viridis(100)
         , main = "25 m"
         , legend = FALSE
    )
    plot(LST.tmp2
         , zlim = c(MIN, MAX)
         , box = FALSE
         , axes = FALSE
         , col = viridis(100)
         , legend.width = 2
         , legend.mar = 10
    )
    plot(LST.tmp1
         , zlim = c(MIN, MAX)
         , box = FALSE
         , axes = FALSE
         , col = viridis(100)
         , legend = FALSE
    )
    mtext(text = paste0(c("MIN","MEAN","MAX")[i], " - "
                        , c("January", "February", "March", "April"
                            , "May", "June", "July", "August",
                            "September", "October", "November", "December")[mm])
          , side = 3, line = -2, at = 0.46, font = 2, outer = TRUE)
    mtext(text = "CHELSA", side = 2, line = -5, at = 0.82, font = 2, outer = TRUE)
    mtext(text = "LST", side = 2, line = -5, at = 0.485, font = 2, outer = TRUE)
    mtext(text = "<---->", side = 2, line = -13.5, at = 0.65, col = "orange", font = 2, outer = TRUE)
    mtext(text = "<---->", side = 4, line = -18, at = 0.65, col = "darkblue", font = 2, outer = TRUE)
    mtext(text = "------------------------------------------>", side = 3, line = -6, at = 0.46, outer = TRUE)
    mtext(text = "<-------------", side = 4, line = -14, at = 0.65, outer = TRUE)
    mtext(text = "<------------------------------------------", side = 1, line = -23, at = 0.46, outer = TRUE)
    
    plot(PIX_TEMP, PIX_LST
         , col = "darkblue", pch = 20
         , xlim = c(MIN, MAX), ylim = c(MIN, MAX)
         , xlab = "CHELSA", ylab = "LST")
    points(PIX_TEMP_ori, PIX_LST_agg, col = "orange", pch = 20)
    abline(v = 0, h = 0, lty = 2)
    abline(lm(PIX_LST ~ PIX_TEMP), lty = 1, lwd = 2, col = "darkblue")
    abline(lm(PIX_LST_agg ~ PIX_TEMP_ori), lty = 1, lwd = 2, col = "orange")
    legend("topleft"
           , legend = c("CHELSA (25m) ~ LST (25m)", "CHELSA (1km) ~ LST (1km)")
           , col = c("darkblue", "orange")
           , pch = 15
           , bty = 'n')
    
    hist(zon.sd[,2], breaks = 1000, main = "Standard deviation\n LST (25m) values within 1000m pixel\n (degrees)")
    
    
    ##
    corr1 = cor.test(PIX_TEMP, PIX_LST)
    corr2 = cor.test(PIX_TEMP, PIX_LST_agg)
    corr3 = cor.test(PIX_TEMP_ori, PIX_LST_agg)
    
    ##
    return(data.frame(MONTH = mm
                      , STAT = c("MIN","MEAN","MAX")[i]
                      , P_VALUE_1 = corr1$p.value
                      , COR_1 = corr1$estimate
                      , P_VALUE_2 = corr2$p.value
                      , COR_2 = corr2$estimate
                      , P_VALUE_3 = corr3$p.value
                      , COR_3 = corr3$estimate
                      ))
}
dev.off()

save(COR, file = "COR.RData")

par(mfrow = c(1,1))
plot(COR[, c(4,6,8)], col = c("darkblue", "darkgreen", "orange"))

library(reshape2)
tmp = melt(COR[, c(2,4,6,8)])
tmp = tmp[which(tmp$variable %in% c("COR_1", "COR_3")), ]
boxplot(value ~ variable * STAT, tmp)

# plot(COR$COR_1, COR$COR_2, col = "coral", pch = 20
#      , xlim = range(COR$COR_1), ylim = range(COR[, c(6,8)]))
# points(COR$COR_1, COR$COR_3, col = "brown", pch = 20)
# abline(a = 0, b = 1, lty = 2, lwd = 2, col = "brown")

