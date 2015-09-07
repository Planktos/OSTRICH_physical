#
#      Functions to process raw GPS data
#
#  (c) Copyright 2013-2014 Jean-Olivier Irisson
#       and Jessica Luo
#      GNU General Public License v3
#   
#   Updated: 6 September 2015 - Kelly Robinson
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

options("digits.secs"=3)

#To process ship's GPS files with file nomenclature "[datetime] Bridge GPSGGA String"
#--------------------------

#list GPS files from ship
gps.files <- list.files("gps_2014_string", recursive = TRUE, full=TRUE)

# Read in functions
# --------------------------------------------
#1 reformat the lat and long in decimal degrees
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

#2 Function to read GPS data from ship
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


# Run functions
  gps <- adply(gps.files, 1, function(file){
    t <- read.gps(file)
    return(t)
  }, .progress="text")
  
  gps <- gps[,-1]
  
  gps.temp <- gps

#round dateTime to the nearest whole second
gps.temp$dateTime <- round_date(gps.temp$dateTime, "second")

gps.sec <- aggregate(cbind(lat, long)~dateTime, data = gps.temp, FUN = mean)

#To process ship's GPS files with file nomenclature "[datetime] Bridge GPSGGA"
# ----------------------------------------------------------------------------
#list GPS files from ship
gps.files <- list.files("gps_2014", recursive = TRUE, full=TRUE)

# function reformat the lat and long in decimal degrees
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

gps.temp$dateTime <- round_date(gps.temp$dateTime, "second")

gps.sec <- aggregate(cbind(lat.gps, long.gps)~dateTime, data = gps.temp, FUN = mean)
