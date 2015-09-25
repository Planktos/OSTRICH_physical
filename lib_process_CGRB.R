#
#      Functions to process raw data
#
#  (c) Copyright 2013-2014 Jean-Olivier Irisson
#       and Jessica Luo
#      GNU General Public License v3
#   
#   Updated: 25 August 2015 - klr
#
#--------------------------------------------------------------------------
# This section performs the following tasks:
# 1) read in the physical data files from a folder
# 2) cleans up the names by removing the extra characters, units and making everything lowercase
# 3) from the time format in the physical data, converts it to the universal time format
# 4) corrects the time zone and adds in transect number, and converts lat/long into decimal degrees

library("lubridate")
library("plyr")
library("stringr")
library("reshape2")
library("pastecs")
library("ggplot2")
library("oce")

options("digits.secs"=3)

# Functions
# -----------
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


#Current code to process ISIIS Physical data files
#------------------------------------------
# create the final data repository

dir.create("data", showWarnings=FALSE)

message("Read and process physical data")


# list all the physical data files in a given directory
phyFiles <- list.files("raw_physical_data_2015", full=TRUE) #use directory name of where physical data files are stored

# read them all
phy <- adply(phyFiles, 1, function(file) {
  
  # read the data
  d <- read.table(file, sep="\t", skip=10, header=TRUE, fileEncoding="ISO-8859-1", stringsAsFactors=FALSE, quote="\"", check.names=FALSE, encoding="UTF-8", na.strings="9999.99")
  
  # clean names
  head <- names(d)
  head <- str_replace(head, "\\(.*\\)", "")
  head <- str_trim(head)
  head <- make.names(head)
  head <- tolower(head)
  head <- str_replace(head, fixed(".."), ".")
  # assign names
  names(d) <- head
  
  # create a proper date + time format
  date <- scan(file, what="character", skip=1, nlines=1, quiet=TRUE)
  date <- date[2]
  mm <- str_sub(date,1,2)
  dd <- str_sub(date,4,5)
  dd <- as.numeric(dd)
  yy <- str_sub(date,7,8)
  dateNextDay <- str_c(mm,as.character(dd+1),yy, sep="/")
  
  # shift by one day when we cross midnight
  d$hour <- as.numeric(str_sub(d$time,1,2))
  d$date <- date
  d$dateTime <- str_c(d$date, d$time, sep=" ")
  d$dateTime <- as.POSIXct(strptime(d$dateTime, format="%m/%d/%y %H:%M:%OS", tz="America/New_York"))
  # Set time zone of collection location
  # NB: we say it is America/New York when it is in fact local time, just to avoid having to deal with time zones
  # shift time of all physical data if not collected in Eastern Time Zone by appropriate number of hourse
  # if collected in the Eastern Time Zone, then set change value to 0
  d$dateTime <- d$dateTime - 0 * 3600 
  
  # detect midnight shifts
  midnightShift <- which(diff(d$dateTime) < 0)
  if (length(midnightShift) > 0) {
    d$dateTime[midnightShift:nrow(d)] <- d$dateTime[midnightShift:nrow(d)] + 24 * 3600
  }
  
  d$date <- NULL
  
  # code in a transect number
  # use the file name as a dummy variable for transect number. Will assign proper transect number later in the pipeline.
  d$transect <- basename(file)
  
  # reformat the lat and long in decimal degrees
  
  # KR: Modifying because JL original wasn't quite working with structure of 2014 physical data files. This may not be robust for data
  # collected in regions with single digit lat and longitude coordinates
  to.dec <- function(x) {
    # split in degree, minute, second
    pieces <- str_split_fixed(x, "°|'| ",3) #Can't make R split at the degree symbol
    deg <- substr(pieces[,1], 1, 2)
    min <- pieces[,2]
    # extract orientation (S/N and E/W)
    orientation <-substrRight(pieces[,3], 1) #Changed from pieces[,3] to sec
    # remove orientation to only keep numbers
    sec <- str_replace(pieces[,3], pattern = "[NSEW]", replacement = "")
    # convert to decimal degrees
    dec <- as.numeric(deg) + as.numeric(min) / 60 + as.numeric(sec) / 3600 #replaced 'pieces' with 'deg', 'min', and 'sec'
    # orient the coordinate
    ifelse(orientation %in% c("S", "W"), -dec, dec)
    
    return(dec)
  }
  
  d$lat <- to.dec(d$lat)
  d$long <- to.dec(d$long)
  # we are in the western hemisphere so longitude should be negative
  d$long <- -d$long
  
  # columns that are all zero are not possible
  # they are actually missing data
  # detect them
  totCol <- colSums(d[llply(d, class) == "numeric"])
  allZeroCols <- names(totCol)[which(totCol == 0)]
  # replace the content with NA
  d[,allZeroCols] <- NA
  
  # rename some columns
  d <- rename(d, c("horizontal.vel.in.water"="horizontal.vel",
                   "irrandiance"="irradiance"))
  
  # keep only interesting data
  d <- d[,c("transect", "dateTime","depth", "lat", "long", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "vol.imaged")]
  
  return(d)
  head
}, .progress="text")

# remove adply crap
phy <- phy[,-1]

##{ Check and correct the physical data -----------------------------------

summary(phy)

# the depth gets stuck from time to time and that results in jumps afterwards. Remove those stuck points and reinterpolate the depth linarly using time.
# assign the depths in which the difference before the previous depth is 0 to be NA
phy$depth[which(diff(phy$depth)==0)+1] <- NA

# interpolation
phy$depth <- approx(phy$dateTime, phy$depth, phy$dateTime, method="linear")$y

#calculate seawater density
phy$density <- swRho(salinity = phy$salinity, temperature = phy$temp, pressure = phy$pressure, eos = "unesco")

# remove some erroneous values
phy$oxygen <- ifelse(phy$oxygen < 0, NA, phy$oxygen)
phy$temp <- ifelse(phy$temp <= 0, NA, phy$temp)
phy$salinity <- ifelse(phy$salinity <= 0, NA, phy$salinity)
phy$irradiance <- ifelse(phy$irradiance <= 0, NA, phy$irradiance)
phy$pressure <- ifelse(phy$pressure < 0, NA, phy$presssure)
phy$density <- ifelse(phy$density < 1000, NA, phy$density)

#calculate vol.imaged
#convert horizontal velocity in (mm/s) to (cms/s)
phy$t0 <- phy$horizontal.vel/10
phy$vel_cms <- ifelse(phy$t0 < 0, NA, phy$t0)
phy$t0 <- NULL

#set constants
image.width <- 2048
scan.rate <- 35000
dof <- 50
fov <- 13.5
#calculate px/image
px.image <- scan.rate/image.width

#calculate distanced imaged per second
phy$dist.image <- phy$vel_cms/px.image
#calculate volume imaged in cm3.
phy$vol.image <- dof*fov*phy$dist.image
#calculate volume rate(L/s)
phy$vol.rate <- (px.image*phy$vol.image)/1000
#1m3/sec
phy$m3.sec <- phy$vol.rate/1000

#sec to image 1-m3
print(1/mean(phy$m3.sec, na.rm=T)) #this value should be around 4-7 sec if DPI was towed at 5 knots (2.5 m/s)

summary(phy)

    #1. If lat and long have lots of NAs, then round physical data to the nearest second so it can be merged with the gps data
    #   phyt.sec <- phy
    #   phyt.sec$dateTime <- round_date(phyt.sec$dateTime, "second")

#Read in transect IDs from log sheets (exported from OSTRICH MS Access database table "ISIIS_Table")
transect.names <- read.csv(file = "transect file names.csv", sep=",", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, na.strings="9999.99")
transect.names <- as.data.frame(transect.names)
transect.names$haul <- as.numeric(transect.names$haul)

phy_t <- merge(x=phy, y=transect.names, by.x = "transect", by.y = "physicaldatafilename", all.x=T)

# remove redundant file names
# if GPS on ISIIS is working at all, then use: 
phy_t <- phy_t[,c("dateTime", "depth", "lat", "long", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "density", "vel_cms", "vol.rate", "dist.image","vol.image", "m3.sec","cruise","transect.id","haul","area","tow","region")]

    # 2. if GPS on ISIIS is NOT working, then use: 
    # phyt.sec <- phyt.sec[,c("dateTime", "depth", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "vel_cms", "vol.rate", "dist.image","vol.image", "m3.sec","density", "haul")]

    # 3. Average parameter values by dateTime (i.e. to the nearest whole second)
    # phy.sec <- aggregate(cbind(depth, temp, salinity, pressure, fluoro, oxygen, irradiance, heading, horizontal.vel, vertical.vel, pitch, vel_cms, vol.rate, dist.image, vol.image, m3.sec, density, haul)~dateTime, data = phyt.sec, FUN = mean, na.action = na.pass)

    # summary(phy.sec)

    # 4. If there are NAs in lat and long, then merge GPS data frmae with physical data frame 
    # Add latitude and longitude using ship's GPS data stream. 
    # phys <- merge(x = gps.sec, y =  phy.sec, by = "dateTime", all.y = T)

#rename lat.gps and long.gps columns to be lat and lon so that they will work with the distance.from.start function in 'plot physical data.R'
names(phy_t)[names(phy_t)=="lat"] <- "lat"
names(phy_t)[names(phy_t)=="long"] <- "lon"

#check phys data
summary(phy_t)

#if needed: fix 4 hour time offset if present (start and end times should match phy data set)
phys$dateTime <- phys$dateTime - 4 * 3600

# inspect water mass data. Replace object "phy_t" with "phys" if needed.
phyM <- melt(phy_t, id.vars=c("dateTime"), measure.vars=c("depth", "temp", "salinity", "density","fluoro", "oxygen", "irradiance"))
ggplot(data=phyM) + geom_histogram(aes(x=value)) + facet_wrap(~variable, scales="free")


# look at profiles
ggplot(phy_t) + geom_path(aes(x=temp, y=-depth), alpha=0.5) + facet_wrap(~transect.id)
ggplot(phy_t) + geom_path(aes(x=salinity, y=-depth), alpha=0.5) + facet_wrap(~transect.id)
ggplot(phy_t) + geom_path(aes(x=density, y=-depth), alpha=0.5) + facet_wrap(~transect.id)
ggplot(phy_t) + geom_path(aes(x=fluoro, y=-depth), alpha=0.5) + facet_wrap(~transect.id)
ggplot(phy_t) + geom_path(aes(x=oxygen, y=-depth), alpha=0.5) + facet_wrap(~transect.id)
ggplot(phy_t) + geom_path(aes(x=irradiance, y=-depth), alpha=0.5) + facet_wrap(~transect.id)

 
  #   5. If needed to pull ship GPS: Add transect metadata to phys data frame
  #   transect.names2 <- read.csv(file = "transect file names2.csv", sep=",", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, na.strings="9999.99")
  #   transect.names2 <- as.data.frame(transect.names2)
  #   transect.names2$haul <- as.numeric(transect.names2$haul)
  #   temp <- merge(x=phys, y=transect.names2, by = "haul", all.x=T)
  #   
  #   physM <- temp

#save phy frame (non-averaged data) as R object
save(phy_t, file = "ost15_phy_t.R")

#print summary statistics
options(digits=2)
options(scipen=100)
attach(phy_t) #or physM

phy_t_stats <- cbind(dateTime, depth, lat, lon, temp, salinity, pressure, fluoro, oxygen, irradiance, heading, horizontal.vel, vertical.vel, pitch, density, vel_cms, vol.rate, dist.image, vol.image, m3.sec)
summary <- stat.desc(phy_t_stats)
write.csv(summary, "ost15_physical_data_summary_stats.csv")


# Functions we use sometimes.
#--------------------------------------------
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

# Compute the straight line distance (km) from the starting point of a lat,lon trajectory
dist.from.start <- function(lat, lon) {
  library("oce")
  geodDist(lat1=lat, lon1=lon, lat2=na.omit(lat)[1], lon2=na.omit(lon)[1]) # / 1.852 # use if you want to convert from km to nautical miles
}

# Compute the distance from rsmas
dist.from.rsmas <- function(lat, lon) {
  # TODO should compensate for the length of cable put out and the angle of the cable (i.e. depth of ISIIS)
  library("oce")
  geodDist(lat, lon, lat2=25.731384, lon2=-80.1621017) / 1.852
}
