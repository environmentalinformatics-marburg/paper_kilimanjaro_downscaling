inpath1 <- "C:/Prozessierte_Daten/2_metre_temperature/Output_eot4"
inpath2 <- "C:/Prozessierte_Daten/2_metre_temperature"

#load libraries
library(raster)
if(!require(remote)){install.packages('remote')}
library(remote)

setwd(inpath2)
temp <- "agg_kili_1600pix.tif"
agg_netcdf_kili <- stack(temp)

temp <- "predict_stack.nc"
predict_stack <- stack(temp)

# make EOT, the predictor is the raster with 40x40 pixels, the predicant is the raster we created
# with the data from the Kili-Stations
setwd(inpath1)
eot(agg_netcdf_kili, predict_stack, write.out = TRUE, path.out=inpath1)

