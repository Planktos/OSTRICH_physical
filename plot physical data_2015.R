
#Plot physical data for each transect

library("plyr")
library("stringr")
library("ggplot2")
library("reshape2")
library("grid")
library("gridExtra")
library("scales")

#2015 data:
load("ost15_phy_t.R")


#Call functions:
#-----------------
# Compute the straight line distance (km) from the starting point of a lat,lon trajectory
  dist.from.start <- function(lat, lon) {
    library("oce")
    geodDist(lat1=lat, lon1=lon, lat2=na.omit(lat)[1], lon2=na.omit(lon)[1]) # / 1.852 # use if you want to convert from km to nautical miles
  }

# Spectral colour map from ColorBrewer
  spectral <- function(n=6) {
    library("RColorBrewer")
    rev(brewer.pal(name="Spectral", n=n))
  }
  
  scale_fill_spectral <- function(...) {
    scale_fill_gradientn(colours=spectral(...))
  }
  scale_colour_spectral <- function(...) {
    scale_colour_gradientn(colours=spectral(...))
  }

# Interpolate a slice of data for which the x-axis is a distance in nautical miles
interp.dist <- function(x, y, z, anisotropy=1000, x.step=500, y.step=2.5, smooth=FALSE, theta=0.2, ...) {
    #
    # Interpolate data over a distance coordinate
    #
    # x   vector of distance *IN KILOMETERS*
    # y   vector of depth in m
    # z   vector of measured variable
    # anisotropy  anisotropy ratio between x and y
    # x/y.step    interpolation grid steps in m
    # smooth      boolean, wether to smooth the first interpolation using fields::image.smooth
    # x/y.step.smooth   interpolation grid step for the smoothing
    # grid.smooth intepolation grid for the smoothing, overrides x/y.step.smooth
    # theta       bandwidth for the kernel smoother in fields::image.smooth
  
    library("akima")
    library("reshape2")
  
    # correct x-axis for anisotropy between horizontal and vertical
    x <- x / anisotropy # 'x <- x*1852/ anisotropy' if x unit is nautical miles
  
    # interpolate
    i <- interp(x=x, y=y, z=z, xo=seq(0, max(x), by=x.step/anisotropy), yo=seq(0, max(y), by=y.step), ...)
  
    # smooth
    if ( smooth ) {
      library("fields")
      i <- image.smooth(i, grid=list(x=i$x, y=i$y), theta=theta)
    }
  
    # extract a data.frame
    out <- melt(i$z, varnames=c("x","y"))
    out$x <- i$x[out$x] * anisotropy / 1852
    out$y <- i$y[out$y]
  
    return(out)
  }

## Detect up and down casts in a depth yo
  detect.casts <- function(depth, order=200) {
    # smoothing the depth profile using a moving average and find the turning points
    
    # smooth depths
    library("pastecs")
    depth_avg <- decaverage(-depth, times=3, weights=c(seq(1, order), order+1, seq(order, 1, -1)))
    # plot(depth_avg)
    depth_avg <- as.numeric(pastecs::extract(depth_avg, component="filtered"))
    
    # detect turning points
    TP <- suppressWarnings(turnpoints(depth_avg))
    
    # set cast numbers (different for up and down casts)
    cast <- cumsum(TP$peaks | TP$pits) + 1
    
    # detect which are up and which are down casts:
    # if the first turning point is a peak, then the first cast (and all odd casts) are upcasts
    if ( TP$firstispeak ) {
      # these are the types for
      #              even  & odd   cast numbers
      castTypes <- c("down", "up")
    } else {
      castTypes <- c("up", "down")
    }
    down.up <- castTypes[cast %% 2 + 1]
    
    return(data.frame(cast, down.up))
  }

#rowShift
rowShift <- function(x, shiftLen = 1L) {
  r <- (1L + shiftLen):(length(x) + shiftLen)
  r[r<1] <- NA
  return(x[r])
  }
  #-------------------

#load physical data frame
#2014 data: 
load("ost14_physM.R")

vars <- c("temp", "salinity", "density", "fluoro", "oxygen", "irradiance")

# Plots using dateTime as x-axis variable
# ---------------------------------------
phys.timeM <- physM[,c("depth", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "density", "haul", "dateTime", "transect.id")]
  # if metadata character strings give you trouble, then exclude transect.id from subset
    # phys.time <- physM[,c("depth", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "density", "haul", "dateTime", "transect.id")]

#Run ddply function in subsets because hauls 36 (OST14-L3-Und) and 37 (OST-L3S) irradiance values are all NAs.
dt <- subset(phys.timeM, haul >= 1 & haul =< 35, select=c(depth:transect.id))

  #if metadata character strings give you trouble...
  #dt <- subset(phys.time, haul >= 1 & haul < 36, select=c(depth:dateTime))

#comment out plot i for hauls 36-61 because no irradiance values will cause ggplot to throw an error
dt <- subset(phys.timeM, haul >= 36 & haul <= 61, select=c(depth:transect.id))

#adjust date range values for 62 since it was performed over a 48-h time period
dt <- subset(phys.timeM, haul == 62, select=c(depth:transect.id))
  #replace ("15 min") with ("6 hour") and ("5 min") with ("1 hour")

#adjust date range values for 63 since it was performed over a 5-h time period
dt <- subset(phys.timeM, haul == 63, select=c(depth:transect.id))
  # make date_breaks ("1 hour") and minor breaks ("15 min")

  #exclude erroneous observations at the beginning of haul 36
  #dt <- subset(phys.time, haul == 36, select=c(depth:dateTime))
  #dt <- subset(dt, dateTime >"2014-06-06 00:00:00")

ddply(.data = dt, .variables = "transect.id", function(x){
    
  na.omit(x) 
  
  # x$cast <- detect.casts(x$depth)
  
  #detect temperature gradients (e.g. thermoclines and thermal fronts)
  
#   x$temp.grad <- ifelse(x$temp[i]>1, yes = (x$temp[i]-x$temp[i-1])/(x$depth[i]-x$depth[i-1]), NA)
#   MAX.temp.grad <- max(x$temp.grad, na.rm=TRUE)
#   MIN.temp.grad <- min(x$temp.grad, na.rm=TRUE)
  
   dm <- melt(x, id.vars=c("dateTime", "depth"), measure.vars=vars)
  
  t <- ggplot(dm[dm$variable == "temp",], aes(x=(dateTime), y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA)+
    scale_x_datetime("Time", labels = date_format("%H:%M"), breaks = date_breaks("15 min"), minor_breaks = "5 min") +  
    scale_y_continuous("depth", expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))

  s <- ggplot(dm[dm$variable == "salinity",], aes(x=dateTime, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_datetime("Time", labels = date_format("%H:%M"), breaks = date_breaks("15 min"), minor_breaks = "5 min") +
    scale_y_continuous("depth", expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
  
  d <- ggplot(dm[dm$variable == "density",], aes(x=dateTime, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_datetime("Time", labels = date_format("%H:%M"), breaks = date_breaks("15 min"), minor_breaks = "5 min") +
    scale_y_continuous("depth", expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
  
  f <- ggplot(dm[dm$variable == "fluoro",], aes(x=dateTime, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_datetime("Time", labels = date_format("%H:%M"), breaks = date_breaks("15 min"), minor_breaks = "5 min") +
    scale_y_continuous("depth", expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))

  o <- ggplot(dm[dm$variable == "oxygen",], aes(x=dateTime, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_datetime("Time", labels = date_format("%H:%M"), breaks = date_breaks("15 min"), minor_breaks = "5 min") +
    scale_y_continuous("depth", expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
#  
#   i <- ggplot(dm[dm$variable == "irradiance",], aes(x=dateTime, y=-depth)) +
#     geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
#     scale_colour_gradient(high=spectral(), na.value=NA) +
#     scale_x_datetime("Time", labels = date_format("%H:%M"), breaks = date_breaks("15 min"), minor_breaks = "5 min") +
#     scale_y_continuous("depth", expand=c(0.01,0.01)) + facet_grid(variable~.) +
#     theme(strip.text.y = element_text(size = 10)) +
#     theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
  
g <- grid.arrange(t, s, d, f, o, ncol=1) #remove 'i' for hauls 36-63

#print image files to a directory
png(file = paste0("plots_time/",unique(x$transect.id), ".png"), width = 8.5, height = 14, units = "in", res = 300)
plot(g)
dev.off()

}, .progress = "text")

# Geotile plots (interpolated) - work in progress -kr 4 sept 2015
#-------------------------------------
# Tile plots for undulation transects - work in progress KR
#-------------------------------------
physMund <- physM[,c("depth", "lat", "lon", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "density", "haul", "dateTime", "transect.id", "tow")]

und <- subset(physMund, tow == "und", select=c(depth:tow))

#Run ddply function in subsets because hauls 36 (OST14-L3-Und) and 37 (OST-L3S) irradiance values are all NAs.
d.und <- subset(und, haul >= 1 & haul <= 35, select=c(depth:transect.id))

d.und <- subset(und, haul == 12, select=c(depth:transect.id))

#if metadata character strings give you trouble...
#dt.und <- subset(phys.time, haul >= 1 & haul < 36, select=c(depth:dateTime))

#comment out plot i for hauls 36-61 because no irradiance values will cause ggplot to throw an error
d.und <- subset(und, haul >= 36 & haul <= 61, select=c(depth:transect.id))

# #adjust date range values for 62 since it was performed over a 48-h time period
# dt.und <- subset(phys.timeM, haul == 62, select=c(depth:transect.id))
# #replace ("15 min") with ("6 hour") and ("5 min") with ("1 hour")
# 
# #adjust date range values for 63 since it was performed over a 5-h time period
# dt.und <- subset(phys.timeM, haul == 63, select=c(depth:transect.id))
# # make date_breaks ("1 hour") and minor breaks ("15 min")


ddply(.data = d.und, .variables = "transect.id", function(x){
  
  na.omit(x) 
  
  d.und$distance <- dist.from.start(d.und$lat, d.und$lon)
  
  dm <- melt(d.und, id.vars=c("distance", "depth"), measure.vars=vars)
  
    di <- ddply(dm, ~variable, function(x) {
      #x <- na.omit(x[which(x$down.up=="up"),])
      xi <- interp.dist(x=x$distance, y=x$depth, z=x$value, duplicate="mean", x.step=300, y.step=1, anisotropy=1000)
    })
    di <- rename(di, c("x"="distance", "y"="depth"))
    
    plots <- dlply(di[di$variable != "Irrandiance.UE.cm",], ~variable, function(x) {
      ggplot(x, aes(x=distance, y=-depth)) +
        # geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
        geom_tile(aes(fill=value), na.rm=T) +
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(colours=spectral(), na.value=NA) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0))
    })
   

    t <- ggplot(di[di$variable == "temp",], aes(x=distance, y=-depth)) +
      geom_tile(aes(fill=value), na.rm=T) +
      #stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
      scale_fill_gradient(high = spectral(), na.value=NA) +
      scale_x_continuous(expand=c(0,0), "distance (km)", breaks = 1, minor_breaks = 0.5) +
      scale_y_continuous(expand=c(0,0), "depth (m)") + facet_grid(variable~.)
  
    t <- ggplot(di[di$variable == "temp",], aes(x=dateTime, y=-depth)) +
      geom_tile(aes(fill=value), na.rm=T) +
      #stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
      scale_fill_gradient(high = spectral(), na.value=NA) +
      scale_x_continuous("dateTime", expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) + facet_grid(variable~.) +
      theme(strip.text.y = element_text(size = 10)) +
      theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
  
    s <- ggplot(di[di$variable == "salinity",], aes(x=distance, y=-Depth.m)) +
      geom_tile(aes(fill=value), na.rm=T) +
      stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
      scale_fill_gradient(colours=spectral(), na.value=NA) +
      scale_x_continuous("distance (km)", expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) + facet_grid(variable~.)
  
  grid.arrange(t, s, f, o, ncol=1)
  
  do.call(grid.arrange, c(plots,list(ncol=1)))
}

# Plots by distance from start
# -----------------------------
# Line plots by distance from transect start (km)
# -----------------------------
#drop character fields
physdM <- physM[,c("depth", "lat", "lon", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "density", "haul", "dateTime", "transect.id")]

#Run ddply function in subsets because haul 36 (OST14-L3-Und) has no irradiance
dd <- subset(physdM, haul >= 1 & haul <= 35, select=c(depth:transect.id))

dd <- subset(physdM, haul == 12, select=c(depth:transect.id))


#comment out plot i for hauls 36-63 because no irradiance values will cause ggplot to throw an error
dd <- subset(physdM, haul >= 36 & haul <= 63, select=c(depth:transect.id))

dd <- subset(physdM, haul >= 60 & haul <= 63, select=c(depth:transect.id))

#Create and save plots to a directory
ddply(.data = dd, .variables = "transect.id", function(x){
  
  dd$distance <- dist.from.start(dd$lat, dd$lon)
  
  # the lat gets stuck from time to time and that results in jumps afterwards. Remove those stuck points and reinterpolate the depth linarly using time.
    
  # interpolation
  # dd$distance <- approdd(dd$dateTime, dd$distance, dd$dateTime, method="linear")$y
 
   # detect yos
   casts <- detect.casts(dd$depth, order=50)
   dd <- cbind(dd, casts)
   # check
   #ggplot(dd) + geom_point(aes(x=dateTime, y=-depth, colour=down.up))
  
   dd <- dd[order(dd$dateTime),]
  
  #detect temperature gradients (e.g. thermoclines and thermal fronts)
  dd$temp.gradient <- abs((dd$temp-rowShift(dd$temp, -1))/(dd$depth-rowShift(dd$depth,-1)))
  dd$density.gradient <- abs((dd$density-rowShift(dd$density, -1))/(dd$depth-rowShift(dd$depth,-1)))
  dc <- dd[ which(is.finite(dd$temp.gradient)),]
   
# Need to fix function to find max by cast and bin temperature gradient to a larger depth bin
#   dy <- ddply(.data = dc, .variables = ~cast, function(z) {
#     z$temp.grad <- (z$temp-rowShift(z$temp, -1))/(z$depth-rowShift(z$depth,-1))
#     max.temp.grad <- max(z$temp.grad, na.rm=TRUE)
#     min.temp.grad <- min(z$temp.grad, na.rm=TRUE)
#     z$s.temp.grad <- ifelse(z$temp.grad == max.temp.grad, 1, 0)
#   })
#   
#   dy <- rename(dy, c("temp.grad"="temp.gradient", "s.temp.grad" = "thermocline"))

  
  dm <- melt(dc, id.vars=c("depth", "distance"), measure.vars=vars)
    
# Create plots using 'distance.km' as the x-axis value
  t <- ggplot(dm[dm$variable == "temp",], aes(x=distance, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance (km)", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))

  p <- ggplot(dm[dm$variable == "density.gradient",], aes(x=distance, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance (km)", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))

  g <- ggplot(dm[dm$variable == "temp.gradient",], aes(x=distance, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance (km)", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
  
  s <- ggplot(dm[dm$variable == "salinity",], aes(x=distance, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance (km)", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))

  d <- ggplot(dm[dm$variable == "density",], aes(x=distance, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance (km)", expand=c(0.01,0.01)) +
    scale_y_continuous("depth", expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
  
  f <- ggplot(dm[dm$variable == "fluoro",], aes(x=distance, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance(km)", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))

  o <- ggplot(dm[dm$variable == "oxygen",], aes(x=distance, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance (km)", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))
# 
#   i <- ggplot(dm[dm$variable == "irradiance",], aes(x=distance, y=-depth)) +
#     geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
#     scale_colour_gradient(high=spectral(), na.value=NA) +
#     scale_x_continuous("distance (km)", expand=c(0.01,0.01)) +
#     scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
#     theme(strip.text.y = element_text(size = 10)) +
#     theme(axis.text = element_text(size = 10)) + theme(axis.title = element_text(size = 12))

g <- grid.arrange(t, s, d, f, o, ncol=1) #remove i for hauls 36-63

#print image files to a directory
png(file = paste0("plots_distance/",unique(x$transect.id), ".png"), width = 8.5, height = 14, units = "in", res = 300)
plot(g)
dev.off()

}, .progress = "text")
