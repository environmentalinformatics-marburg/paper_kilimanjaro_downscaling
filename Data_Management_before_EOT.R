inpath <- "C:/2_metre_temperature/"
setwd(inpath)

# install required Packages 
if (!require(CDFt)){install.packages('CDFt')}
library(CDFt)
if (!require(raster)){install.packages('raster')}
library(raster)
if (!require(ncdf)){install.packages("ncdf")}
library(ncdf4)
library(ncdf)
#type = "source", configure.args="--with-netcdf-include=/usr/include"
if (!require(rgdal)){install.packages('rgdal')}
library(rgdal)
if (!require(chron)){install.packages('chron')}
library(chron)
if(!require(caret)){install.packages('caret')}
library(caret)



temp <- "netcdf-web224-20150223144140-8541-12767.nc"
netcdf <- stack(temp)

# Data Management: Here we need to clean up the raster stack - since there is information which is not needed
# => 3 years with 4 Measurments a day => 4384 Rasters - after those there are different timestamps saved in the original 
# Raster Stack => they are not needed anymore

netcdf <- netcdf[[1:4384]]

df <- data.frame(names(netcdf),seq(1:4384))
df$time <- as.character(strptime(names(netcdf), "X%Y.%m.%d.%H.%M.%S"))
df$timeday <- substr(df$time, 1, 10)
df$timehour <- substr(df$time, 12,19)

df$timehour <- times(df$timehour)
badtimes <- times(c(df$timehour[1977:1980]))
a <- times("01:00:00")

for (i in 1:length(df$timehour)){
  if (df$timehour[i] == badtimes[1]) {
 
    df$timehour[i] <- (df$timehour[i] - a)
  } else if (df$timehour[i] == badtimes[2]){
 
    df$timehour[i] <- (df$timehour[i] - a)
  } else if (df$timehour[i] == badtimes[3]) {
 
    df$timehour[i] <- (df$timehour[i] - a) 
  } else if (df$timehour[i] == badtimes[4]) {
   
    df$timehour[i] <- (df$timehour[i] - a)
  }
}

df$time <- paste0(df$timeday,".",df$timehour)
colnames(df)[1] <- "names.netcdf"
df[,2] <- NULL
df$names.netcdf <- paste0("X",df$timeday,"-",df$timehour)

names(netcdf) <- df$names.netcdf

row <- rowFromY(netcdf, -3.076475)
col <- colFromX(netcdf, 37.353205)
# create Rastermask for the Kili Region. 
raster_mask <- netcdf[[1]]
#write out the original Raster for control
#writeRaster(raster_mask, "where_is_kili", format="GTiff", overwrite=TRUE)
vals <- 0
raster_mask <- setValues(raster_mask, vals)

# get position of thes Cells around the Kili Pixel
cellntl <- cellFromRowCol(raster_mask, row-20, col-20)
cellnbl <- cellFromRowCol(raster_mask, row+20, col-20)
cellntr <- cellFromRowCol(raster_mask, row-20, col+20)
cellnbr <- cellFromRowCol(raster_mask, row+20, col+20)
ext <- data.frame(xyFromCell(raster_mask, cellntl)) #topleft
ext <- rbind(ext,xyFromCell(raster_mask, cellnbl))  #bottomleft
ext <- rbind(ext, xyFromCell(raster_mask, cellntr)) #topright
ext <- rbind(ext, xyFromCell(raster_mask, cellnbr)) # bottomright


# moving edgepoints to middle cell to get extent
ext[1,1] <- ext[1,1]-0.375
ext[2,1] <- ext[2,1]-0.375
ext[3,1] <- ext[3,1]+0.375
ext[4,1] <- ext[4,1]+0.375
ext[1,2] <- ext[1,2]+0.375
ext[2,2] <- ext[2,2]-0.375
ext[3,2] <- ext[3,2]+0.375
ext[4,2] <- ext[4,2]-0.375

edgepoints2 <- SpatialPointsDataFrame(ext, ext)
crs(edgepoints2) <- crs(raster_mask)
writeOGR(edgepoints2, inpath, "edgepoints", driver="ESRI Shapefile", overwrite=TRUE)

ext <- extent(ext)
# extract 20x20 matrix around kili
netcdf_kili <- crop(netcdf, ext)

len <- (length(df$names.netcdf)/4)


# here we aggregate the pixels of the netcdf files to days, there was a 6 h rhythm before
agg_4_times <- netcdf_kili[[1]]
agg_netcdf_kili <- netcdf_kili[[1:1096]]
for (i in (seq(1:len))){
  print(i)
  agg_4_times <- mean(netcdf_kili[[i:(i+3)]])
  agg_netcdf_kili[[i]] <- agg_4_times 
}

df_names <- aggregate(df,by=list(df$timeday), FUN=mean)
df_names <- df_names[,1]
names(agg_netcdf_kili) <- df_names

writeRaster(agg_netcdf_kili, "agg_kili_1600pix", format="CDF", overwrite=TRUE)

# take the plots from the kili and make a rasterstack out of it so it can be used in the EOT
temp_hourly <- read.csv("plots.csv")
temp_hourly$day <- substr(temp_hourly$datetime, 1, 10)

# aggregate on the daily temperature
temp_daily <- aggregate(temp_hourly$Ta_200~temp_hourly$day+temp_hourly$plotID,temp_hourly,mean, na.action = na.pass)

# transform the stations .shp file to wgs 84, the netcdf from emwcf is in wgs 84 too. 
station_shp <- readOGR("C:/2_metre_temperature/plots.shp","plots")
station_wgs <- spTransform(station_shp, crs(agg_netcdf_kili))

colnames(temp_daily) <- c("day", "ID", "Ta_200")
write.csv(temp_daily, "test.csv")

# produce a matrix/raster which is used to print the stations-datasets on
xy <- matrix(NA,200,200)
z <- matrix(NA,200,200)
rast <- raster(xy)
rast1 <- raster(z)
extent(rast) <- extent(station_wgs)
projection(rast) <- crs(agg_netcdf_kili)
# extract the extent on the extent of the stations
extent(rast1) <- extent(station_wgs)
projection(rast1) <- crs(agg_netcdf_kili)

# create a raster with NA's (1096 Layers)
predict_stack <- stack(rast,rast1)
for (i in (1:len)){
  print(i)
  predict_stack[[i]] <- rast
}


# transform the stations shp in a data.frame
# We have to match the PlotID'S from the shp to the ID's from the data.frame
station_wgs_df <- data.frame(station_wgs)
station_wgs_df$PlotID <- as.character(station_wgs_df$PlotID)
station_wgs_df <- station_wgs_df[station_wgs_df$PlotID %in% levels(temp_daily$ID),]

calccelln <- function (x,y, raster){
  # this function creates col and row numebers from a given x and y in a raster
  row <- rowFromY(raster,  y)
  col <- colFromX(raster, x)
  celln <- cellFromRowCol(raster, row, col)
  return(celln)
}  


# here we create new cols in the data.frame where we put in the raster col/rownumbers of the stations 
for (i in 1:(length(station_wgs_df$PlotID))){
  station_wgs_df$celln[i] <- calccelln(station_wgs_df$coords.x1[i], station_wgs_df$coords.x2[i], rast)
}



# Now we put the stationsdata in the rasterstack. The for-loops are created this way:
# take a station, take the raster col/rownumber of the station, put all 1096 in the right rasterstack
# LONG COMPUTATION TIME!!!!
i <- 1
m <- 1
for (i in (1:length(station_wgs_df$PlotID))){
  print(i)
  act_station <- station_wgs_df$PlotID[[i]]
  tmp_df <- temp_daily[temp_daily$ID==act_station,]
  celln <- station_wgs_df$celln[i]
  for (m in(1:len)){
    act_raster <- predict_stack[[m]]
    act_raster[celln] <- tmp_df$Ta_200[m]
    predict_stack[[m]] <- act_raster
    }
}


# write out raster - we can start with the eot now.

test <- (predict_stack[[1]])
extent(predict_stack) <- extent(rast1) 
writeOGR(station_wgs, inpath, "station_wgs", driver="ESRI Shapefile", overwrite = TRUE)
writeRaster(predict_stack, "predict_stack", format="GTiff", overwrite=TRUE)
writeRaster(agg_netcdf_kili, "agg_kili_1600pix", format="GTiff", overwrite=TRUE)



# Probleme mit dem foreach Package - deshalb erst dann verwenden wenn es gebraucht wird.
#if(!require(remote)){install.packages('remote')}
#library(remote)
#??remote
#?eot
#eot(predict_stack, agg_netcdf_kili)
