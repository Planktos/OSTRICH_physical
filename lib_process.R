#
#      Functions to process raw data
#
#  (c) Copyright 2013-2014 Jean-Olivier Irisson
#       and Jessica Luo
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

# reformat the lat and long in decimal degrees
to.dec <- function(x) {
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

# Read GPS data from ship
read.gps <- function(file) {
  library(stringr)
  
  options(digits.secs=3)  # allow split seconds
  
  # read table
  t <- read.table(file, header=F, skip=2, sep=",")
  # split the first column
  tmp <- colsplit(t$V1, pattern="\t", names=c("date", "time", "model"))
  # bind together
  t <- cbind(tmp, t[2:ncol(t)])
  t <- t[,names(t) %in% c("date", "time", "V3", "V4", "V5", "V6")]
    
  names(t) <- c("date", "time", "lat", "N-S", "long", "E-W")
    
  t$latdeg <- str_sub(t$lat,1,2)
  t$latminsec <- str_sub(t$lat,3,-1)
  t$londeg <- str_sub(t$long,1,2)
  t$lonminsec <- str_sub(t$long,3,-1)
    
  # convert to degrees, take into account the N/S/E/W orientation
  t$lat <- as.numeric(t$latdeg) + as.numeric(t$latminsec)/60
  t$lat <- ifelse(t$"N-S" == "S", -t$lat, t$lat)
  t$long <- as.numeric(t$londeg) + as.numeric(t$lonminsec)/60
  t$long <- ifelse(t$"E-W" == "W", -t$long, t$long)
    
  # combine date and time and convert dateTime to POSIX
  t$dateTime <- str_c(t$date, " ",t$time)
    
  t$dateTime <- as.POSIXct(t$dateTime, format="%m/%d/%Y %H:%M:%OS", tz="America/New_York")
    
  # convert from GMT to local time
  t$dateTime <- t$dateTime - 4*3600
  # keep only relevant columns
  t <- t[,names(t) %in% c("dateTime", "lat", "long")]
    
  return(t)  
  
}

# Read hydrological data from ISIIS
read.isiis <- function(file) {
  library("stringr")

	options(digits.secs=3)  # allow split seconds #changed from 2 to 3

	# read the data
  d <- read.delim(file, skip=10, fileEncoding="ISO-8859-1", encoding="UTF-8", stringsAsFactors=FALSE)

  # clean names
	# remove double dots
  names(d) <- str_replace_all(names(d), fixed(".."), ".")
  names(d) <- str_replace_all(names(d), fixed(".."), ".")
	# remove dots at end of names
  names(d) <- str_replace(names(d), "\\.$", "")

  # extract date from file name
  year <- str_sub(file, -16, -13)
  month <- str_sub(file, -12, -11)
  day <- str_sub(file, -10, -9)

  # compute date and time
  d$dateTimeMsec <- as.POSIXct(str_c(year, "-", month, "-", day, " ", d$Time), tz="America/New_York")
  # detect midnight shifts
  midnightShift <- which(diff(d$dateTimeMsec) < 0)
  if (length(midnightShift) > 0) {
    d$dateTimeMsec[midnightShift:nrow(d)] <- d$dateTimeMsec[midnightShift:nrow(d)] + 24 * 3600
  }
  
  # keep important columns
  d <- d[,c("dateTimeMsec", "Pressure.dbar", "Depth.m", "Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm", "Lat.decimals", "Long.decimals")]

	return(d)
}

# Detect up and down casts in a depth yo
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

# Compute the straight line distance from the starting point of a lat,lon trajectory
dist.from.start <- function(lat, lon) {
	library("oce")
	geodDist(lat1=lat, lon1=lon, lat2=na.omit(lat)[1], lon2=na.omit(lon)[1]) / 1.852
}

# Compute the distance from rsmas
dist.from.rsmas <- function(lat, lon) {
	# TODO should compensate for the length of cable put out and the angle of the cable (i.e. depth of ISIIS)
	library("oce")
	geodDist(lat, lon, lat2=25.731384, lon2=-80.1621017) / 1.852
}

