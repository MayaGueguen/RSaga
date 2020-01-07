

###########################################################################
### LAPSE RATE : ERA5 data
###########################################################################

# equation to find Geopotential altitude (m) in function of Pressure levels (hPa) according to ERA-interim values :

``` R
ff = function(x){ 0.0093*x*x - 24.9763*x + 15850.56 }
xx = c(seq(250, 750, 50), seq(775, 1000, 25))
length(xx)
sapply(xx, ff)
```

 [1] 10187.735  9194.670  8248.105  7348.040  6494.475  5687.410  4926.845
 [8]  4212.780  3545.215  2924.150  2349.585  2079.740  1821.520  1574.925
[15]  1339.955  1116.610   904.890   704.795   516.325   339.480   174.260


###########################################################################
## IF SPACE IS NEEDED !
## Some files can be erased and quite quickly reproduced.
## 	- CLOUDS/[user-defined]_Mercator_resolution[user-defined]/
## 	- PRECIPITATION/[user-defined]_ETRS89_resolution[user-defined]/
## 	- SOLAR_RADIATION/[user-defined]_Mercator_resolution[user-defined]/Quotient...
## 	- TEMPERATURE/[user-defined]_Mercator_resolution[user-defined]/
