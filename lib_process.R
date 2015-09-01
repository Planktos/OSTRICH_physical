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
# 1) read in the GPS data files from a folder
# 2) cleans up the names by removing the extra characters, units and making everything lowercase
# 3) from the time format in the gps data, converts it to the universal time format
# 4) corrects the time zone and converts lat/long into decimal degrees

library("lubridate")
library("plyr")
library("stringr")
library("reshape2")
library("pastecs")
library("ggplot2")

options("digits.secs"=3)

#list GPS files from ship
gps.files <- list.files("gps_2014", recursive = TRUE, full=TRUE)

# reformat the lat and long in decimal degrees
to.dec.gps <- function(y) {
  # split in degree, minute, second
  pieces <- str_split_fixed(x, "° |'", 3)
  # extract orientation (S/N and E/W)
  orientation <- str_sub(pieces[,3], -1)
  # remove orientation to only keep numbers
  pieces[,3] <- str_replace(pieces[,3], "[NSEW]", "")
  # convert to decimal degrees
  dec <- as.numeric(pieces[,1]) + as.numeric(pieces[,2]) / 60 + as.numeric(pieces[,3]) / 3600
  # orient the coordinate
  ifelse(orientation %in% c("S", "W"), -dec, dec)
  
  return(dec)
}

# Read ship's GPS data
# replaced read.gps with an adply call - kr
gps <- adply(gps.files, 1, function(file) {

# read.gps <- function(file) {
  library(stringr)
  
  options(digits.secs=3)  # allow split seconds
  
  # read table
  t <- read.table(file, header=F, skip=2, sep="\t",fileEncoding="UTF-8", skipNul = TRUE)
  #name the fields
  colnames(t) <- c("date", "time", "gpsTime", "latdeg", "latminsec", "N-S", "londeg", "lonminsec","E-W", "altitude","height","dilution","satellites", "fix.quality","model", "checksum")
  # use the following instead of the above if reading in GPSVTG String .DAT GPS files
    # split the first column
      # tmp <- colsplit(t$V1, pattern="\t", names=c("date", "time", "model"))
    # bind together
      #   t <- cbind(tmp, t[2:ncol(t)])
    #   t <- t[,names(t) %in% c("date", "time", "V3", "V4", "V5", "V6")]
    #     
    #   names(t) <- c("date", "time", "lat", "N-S", "long", "E-W")
        
    #   t$latdeg <- str_sub(t$lat,1,2)
    #   t$latminsec <- str_sub(t$lat,3,-1)
    #   t$londeg <- str_sub(t$long,1,2)
    #   t$lonminsec <- str_sub(t$long,3,-1)
    
  # convert to degrees, take into account the N/S/E/W orientation
  t$lat.gps <- as.numeric(t$latdeg) + as.numeric(t$latminsec)/60
  t$lat.gps <- ifelse(t$"N-S" == "S", -t$lat.gps, t$lat.gps)
  t$long.gps <- as.numeric(t$londeg) + as.numeric(t$lonminsec)/60
  t$long.gps <- ifelse(t$"E-W" == "W", -t$long.gps, t$long.gps)
    
  # combine date and time and convert dateTime to POSIX
  t$dateTime <- str_c(t$date, " ",t$time)
  t$dateTime <- as.POSIXct(t$dateTime, format="%m/%d/%Y %H:%M:%OS", tz="GMT")
  # convert from GMT to local time (Eastern Standard Time)
  t$dateTime <- t$dateTime - 4*3600
 
 # keep only relevant columns & reorder them
  t <- t[,names(t) %in% c("dateTime", "lat.gps", "long.gps")]
  t <- t[c("dateTime", "lat.gps", "long.gps")]

  return(t)
}, .progress="text")

# remove adply crap
gps <- gps[,-1]

gps.temp <- gps

#round dateTime to the nearest whole second
# binSize <- 0.01

#calculate the number of bins
# Try this later to keep subsecond accuracy using computer with more RAM 8/27/2015 kr
#   gps.temp.maxT <- max(gps.temp$dateTime, na.rm=T)
#   gps.temp.minT <- min(gps.temp$dateTime, na.rm=T)
#   gps.bins=seq(gps.temp.minT, gps.temp.maxT, by = binSize)
#   gps.temp$dateTimeB <- cut(gps.temp$dateTime, breaks=gps.bins, labels=1:(length(gps.bins)-1))
#   gps.temp$dateTimeB <- as.numeric(gps.temp$dateTimeB)
#   gps.dateTimeBin <- aggregate(dateTime~dateTimeB, data = gps.temp, FUN = mean)

gps.temp$dateTimeR <- round_date(gps.temp$dateTime, "second")

gps.sec <- aggregate(cbind(dateTime, lat.gps, long.gps)~dateTimeR, data = gps.temp, FUN = mean)


# # Read hydrological data from ISIIS
#-------------------------------------
# read.isiis <- function(file) {
#   library("stringr")
# 
# 	options(digits.secs=3)  # allow split seconds #changed from 2 to 3
# 
# 	# read the data
#   d <- read.delim(file, skip=10, fileEncoding="ISO-8859-1", encoding="UTF-8", stringsAsFactors=FALSE)
# 
#   # clean names
# 	# remove double dots
#   names(d) <- str_replace_all(names(d), fixed(".."), ".")
#   names(d) <- str_replace_all(names(d), fixed(".."), ".")
# 	# remove dots at end of names
#   names(d) <- str_replace(names(d), "\\.$", "")
# 
#   # extract date from file name
#   year <- str_sub(file, -16, -13)
#   month <- str_sub(file, -12, -11)
#   day <- str_sub(file, -10, -9)
# 
#   # compute date and time
#   d$dateTimeMsec <- as.POSIXct(str_c(year, "-", month, "-", day, " ", d$Time), tz="America/New_York")
#   # detect midnight shifts
#   midnightShift <- which(diff(d$dateTimeMsec) < 0)
#   if (length(midnightShift) > 0) {
#     d$dateTimeMsec[midnightShift:nrow(d)] <- d$dateTimeMsec[midnightShift:nrow(d)] + 24 * 3600
#   }
#   
#   # keep important columns
#   d <- d[,c("dateTimeMsec", "Pressure.dbar", "Depth.m", "Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm", "Lat.decimals", "Long.decimals")]
# 
# 	return(d)
# } #original "lib_process" code to process ISIIS physical data
# -----------------------------------------

# create the final data repository
dir.create("data", showWarnings=FALSE)

message("Read and process physical data")

##{ Read and reformat physical data ---------------------------------------

# This section performs the following tasks:
# 1) read in the physical data files from a folder
# 2) cleans up the names by removing the extra characters, units and making everything lowercase
# 3) from the time format in the physical data, converts it to the universal time format
# 4) corrects the time zone and adds in transect number, and converts lat/long into decimal degrees

# list all the physical data files in a given directory
phyFiles <- list.files("raw_physical_data_2014", full=TRUE)

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
  # d$date <- ifelse(d$hour >= 18 & d$hour <= 23, date, dateNextDay) #NEED to ask Jessica to explain this line
  d$dateTime <- str_c(d$date, d$time, sep=" ")
  d$dateTime <- as.POSIXct(strptime(d$dateTime, format="%m/%d/%y %H:%M:%OS", tz="America/New_York"))
  # Set time zone of collection location
  # NB: we say it is America/New York when it is in fact local time, just to avoid having to deal with time zones
  
  # shift time of all physical data if not collected in Eastern Time Zone by appropriate number of hourse
  # if collected in the Eastern Time Zone, then set change value to 0
  d$dateTime <- d$dateTime - 0 * 3600  
  
  # code in a transect number
  # use the file name as a dummy variable for transect number. Will assign proper transect number later in the pipeline.
  d$transect <- basename(file)
    # San Diego transect assignment
      # this is not robust for all physical data but is necessary here
      #d$transect <- dd-14
  
  # reformat the lat and long in decimal degrees
  
  # KR: Modifying because JL original wasn't quite working with structure of 2014 physical data files. This may not be robust for data
  # collected in regions with single digit lat and longitude coordinates
  to.dec <- function(x) {
    # split in degree, minute, second
    #pieces <- str_split_fixed(x, "°|'",2) #Can't make R split at the degree symbol
    deg <- substr(x, 1, 2)
    min <- substr(x, 5, 6)
    sec <- substr(x, 8, 12)
        # extract orientation (S/N and E/W)
    orientation <- substr(x, 13,13) #Changed from pieces[,3] to sec
    # remove orientation to only keep numbers
      #pieces[,2] <- str_replace(pieces[,3], pattern = "[NSEW]", replacement = "")
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
  d <- d[,c("transect", "dateTime", "depth", "lat", "long", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "vol.imaged")]
  
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

#Read in transect IDs from log sheets (exported from OSTRICH MS Access database table "ISIIS_Table")
transect.names <- read.csv(file = "transect file names.csv", sep=",", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, na.strings="9999.99")
transect.names <- as.data.frame(transect.names)
phyt <- merge(x=phy, y=transect.names, by.x = "transect", by.y = "physicaldatafilename", all.x=T)

# remove redundant file names
  # if GPS on ISIIS is working well, then use: 
  # phyt <- phyt[,c("dateTime", "depth", "lat", "long", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "vol.imaged", "cruise", "haul.no", "transect.id")]
  # if GPS on ISIIS is NOT working, then use: 
  phyt <- phyt[,c("dateTime", "depth", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "haul")]
  # Add latitude and longitude using ship's GPS data stream

#round physical data to the nearest second so it can be merged with the gps data
phyt$dateTime <- round_date(phyt$dateTime, "second")
phyt.sec <- aggregate(cbind(depth, temp, salinity, pressure, fluoro, oxygen, irradiance, heading, horizontal.vel, vertical.vel, pitch, haul)~dateTime, data = phyt, FUN = mean)

summary(phyt.sec)

#merge gps $ physical data
phys <- merge(x = gps.sec, y =  phyt.sec, by.x = "dateTimeR", by.y = "dateTime", all.y = T)

#check phys data
summary(phys)

#fix 4 hour time offset if present (start and end times should match phy data set)
phys$dateTimeR <- phys$dateTimeR - 4 * 3600


#check to make sure lat and lon from ship matched up correctly with physical data via timestamps
# wm1 <- subset(phys, haul == 3, select=c(dateTimeR:haul))
# eund1 <- subset(phys, haul == 4, select=c(dateTimeR:haul))
# ws1 <- subset(phys, haul == 2, select=c(dateTimeR:haul))
# cs2 <- subset(phys, haul == 9, select=c(dateTimeR:haul))
# cu2 <- subset(phys, haul == 8, select=c(dateTimeR:haul))
# ed1 <- subset(phys, haul == 7, select=c(dateTimeR:haul))


#Figure out why there are NAs in lat and long fields
gps.na <- phys[is.na(phys$lat),]


# inspect water mass data
phyM <- melt(phy, id.vars=c("dateTime"), measure.vars=c("depth", "temp", "salinity", "fluoro", "oxygen", "irradiance"))
ggplot(data=phyM) + geom_histogram(aes(x=value)) + facet_wrap(~variable, scales="free")


# look at profiles
ggplot(phyt) + geom_path(aes(x=temp, y=-depth), alpha=0.5) + facet_wrap(~transect.id)
ggplot(phyt) + geom_path(aes(x=salinity, y=-depth), alpha=0.5) + facet_wrap(~transect.id)
ggplot(phyt) + geom_path(aes(x=fluoro, y=-depth), alpha=0.5) + facet_wrap(~transect.id)
ggplot(phyt) + geom_path(aes(x=oxygen, y=-depth), alpha=0.5) + facet_wrap(~transect.id)
ggplot(phyt) + geom_path(aes(x=irradiance, y=-depth), alpha=0.5) + facet_wrap(~transect.id)
# Not needed
  # ggplot(phyt) + geom_path(aes(x=irradiance, y=-depth), alpha=0.5) + facet_wrap(~transect.id) + scale_x_continuous(limits=c(-1.7E-6, -7.5E-7)) 



#save phys frame as R object
save(phys, "ost14_phy.R")

write.table(phys, file = "ost14-phys.txt", sep = "\t", row.names = F, col.names = T)

# Need to get distance from start individually for each transect for plotting
# subset phys by transect
cs2 <- subset(phys, haul == 9, select=c(dateTimeR:haul))

ed1 <- subset(phys, haul == 7, select=c(dateTimeR:haul))

#get distance
ed1$distance.km <- dist.from.start(ed1$lat, ed1$long)

# Functions
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

