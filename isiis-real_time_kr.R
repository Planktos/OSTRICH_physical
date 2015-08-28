#
#      Read data from the ISIIS computer and Walton Smith on the fly and plot it
#
#  (c) Copyright 2013-2014 Jean-Olivier Irission and Jessica Luo
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library("plyr")
library("stringr")
library("ggplot2")
library("reshape2")
library("grid")
library("gridExtra")

source("lib_process.R")
source("lib_plot.R")

`%ni%` <- Negate(`%in%`) 

# TODO:
# 1. clean up document so that real.time.plot calls another function instead of trying to process all isiis data itself
# 2. see about switching from plyr library to dplyr, might be faster
# 3. create a proper data processing file where all data is processed, combined (if necessary) with shipboard gps data, 
# and then written to one physical_data_processed.csv file
# 4. add meaningful comments
# 5. modify binning code so can bin physical data
# 6. figure out how to bin the fixed depth transects
# 7. make code ready for OSTRICH 2015 cruise


real.time.plot <- function(vars=c("Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm")) {
  # get data
  # hydroFiles <- list.files("/Volumes/ISIIShydro", pattern=glob2rx("ISIIS2014*.txt"), full=T)
  hydroFiles <- list.files("raw_physical_data", pattern=glob2rx("ISIIS2014*.txt"), full=T)
  
  transects <- read.csv("raw_physical_data/transects.csv")
  transects$start <- as.POSIXct(transects$start)
  transects$end <- as.POSIXct(transects$end)
 
  # read in gps files
  gpsFiles <- list.files("raw_physical_data/gps", pattern=glob2rx("*.dat"), full=T)
  
  gps <- adply(gpsFiles, 1, function(file){
    t <- read.gps(file)
    return(t)
  }, .progress="text")
  
  gps <- gps[,-1]
  
  # round gps dateTime
  gps$dateTime <- as.POSIXct(round(gps$dateTime))
  
  
  d <- adply(hydroFiles, 1, function(file){
    d <- read.isiis(file)
    d$dateTime <- as.POSIXct(round(d$dateTimeMsec))
    
    date <- as.numeric(str_sub(file, 30,31))
    
    if (date < 4 | date >= 28){
      d <- d[,which(names(d) %ni% c("Lat.decimals", "Long.decimals"))]
      d <- join(d, gps, by="dateTime")
      
      # interpolate data points
      d$lat <- approx(x=as.numeric(d$dateTime), y=d$lat, xo=as.numeric(d$dateTime))$y
      d$long <- approx(x=as.numeric(d$dateTime), y=d$long, xo=as.numeric(d$dateTime))$y
    } else {
      d$lat <- to.dec(d$Lat.decimals)
      d$long <- to.dec(d$Long.decimals)
      # we are in the western hemisphere so longitude should be negative
      d$long <- -d$long
    }
    
    d <- d[,c("dateTime", "Pressure.dbar", "Depth.m", "Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm", "lat", "long")]
    return(d)
  }, .progress="text")

  d <- d[,-1]

  # remove some erroneous values
  d <- d[-which(d$Salinity.PPT==0),]
  
  
  # compute distance from a reference point
  d$distanceFromMiami <- dist.from.start(d$lat, d$long)
  
  # detect yos
  casts <- detect.casts(d$Depth.m, order=50)
  d <- cbind(d, casts)
  # check
  # ggplot(d) + geom_point(aes(x=dateTime, y=-Depth.m, colour=down.up))
  
  # interpolate all variables
  dm <- melt(d, id.vars=c("Depth.m", "down.up", "distanceFromMiami"), measure.vars=vars)  
  
  di <- ddply(dm, ~variable, function(x) {
    x <- na.omit(x[which(x$down.up=="up"),])
    xi <- interp.dist(x=x$distanceFromMiami, y=x$Depth.m, z=x$value, duplicate="mean", x.step=300, y.step=1, anisotropy=1000)
  })
  di <- rename(di, c("x"="distance", "y"="Depth.m"))
  
  plots <- dlply(di[di$variable != "Irrandiance.UE.cm",], ~variable, function(x) {
    ggplot(x, aes(x=distance, y=-Depth.m)) +
       # geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
       geom_tile(aes(fill=value), na.rm=T) +
       stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
      scale_fill_gradientn(colours=spectral(), na.value=NA) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0))
  })
  
t <- ggplot(di[di$variable == "Temp.C",], aes(x=distance, y=-Depth.m)) +
  geom_tile(aes(fill=value), na.rm=T) +
  stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
  scale_fill_gradientn(colours=spectral(), na.value=NA) +
  scale_x_continuous("distance (km)", expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) + facet_grid(variable~.)

s <- ggplot(di[di$variable == "Salinity.PPT",], aes(x=distance, y=-Depth.m)) +
  geom_tile(aes(fill=value), na.rm=T) +
  stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
  scale_fill_gradientn(colours=spectral(), na.value=NA) +
  scale_x_continuous("distance (km)", expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) + facet_grid(variable~.)

f <- ggplot(di[di$variable == "Fluoro.volts",], aes(x=distance, y=-Depth.m)) +
  geom_tile(aes(fill=value), na.rm=T) +
  stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
  scale_fill_gradientn(colours=spectral(), na.value=NA) +
  scale_x_continuous("distance (km)", expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) + facet_grid(variable~.)

o <- ggplot(di[di$variable == "Oxygen.ml.l",], aes(x=distance, y=-Depth.m)) +
  geom_tile(aes(fill=value), na.rm=T) +
  stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
  scale_fill_gradientn(colours=spectral(), na.value=NA) +
  scale_x_continuous("distance (km)", expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) + facet_grid(variable~.)

grid.arrange(t, s, f, o, ncol=1)

do.call(grid.arrange, c(plots,list(ncol=1)))
}

ggplot(d) + geom_point(aes(x=distanceFromMiami, y=-Depth.m, colour=Temp.C), size=8) + scale_color_spectral() + theme_bw()

ggplot(d) + geom_point(aes(x=distanceFromMiami, y=-Depth.m, colour=Salinity.PPT), size=8) + scale_color_spectral() + theme_bw()

ggplot(d) + geom_point(aes(x=distanceFromMiami, y=-Depth.m, colour=Fluoro.volts), size=8) + scale_color_spectral() + theme_bw()

ggplot(d) + geom_point(aes(x=distanceFromMiami, y=-Depth.m, colour=Oxygen.ml.l), size=8) + scale_color_spectral() + theme_bw()

ggplot(d) + geom_point(aes(y=Temp.C, x=Salinity.PPT, colour=-Depth.m)) + scale_color_spectral() + theme_bw()


gpseveryhr <- gps[seq(1,nrow(gps),by=3600),]
gpseveryhr$times <- str_sub(as.character(gpseveryhr$dateTime),12,16)

gpsevery2hr <- gps[seq(1,nrow(gps),by=17000),]
gpsevery2hr$days <- trunc(gpsevery2hr$dateTime, units="days")

ggplot() + geom_point(aes(x=lon, y=lat), size=0.7, alpha=0.4, data=coast) + 
  geom_point(aes(x=long, y=lat), size=1, colour="red", alpha=0.7, data=gps) + 
  geom_text(aes(x=long, y=lat, label=times), size=4,  data=gpseveryhr) + 
  # geom_text(aes(x=long, y=lat, label=days), size=3, colour="red", data=gpsevery2hr) + 
  scale_x_continuous(limits=c(-80.3,-79.25)) + scale_y_continuous(limits=c(25,26.25)) + theme_bw()


while ( 1 == 1 ) {
  real.time.plot(c("Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l"))
  Sys.sleep(180)
}


write.csv(d, "data/combined_phys_201406021736.csv", row.names=F)

