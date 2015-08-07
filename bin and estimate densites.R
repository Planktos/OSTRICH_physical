#
# Calculate taxa densities (#/m3) for OSTRICH REU
#
##The purpose of this script is to estimate the taxa densities (individuals per m3) in the 1E and 1W sites for OSTRICH 2014.

# 1. reads in physical data files (code from JL) for sites and R objects that have all the biological and physical data fields
# 2. subset phyiscal data to get just horizontal velocity
# 3. estimate the average number of seconds it took to image 1-m3 of water in 1E and 1W
# 4. create time bins for each 1E and 1W data sets using the average number of seconds as the interval
# 5. calculate the total number of each taxon in a given time bin (each representing 1-m3), so individuals/m3
# 6. merge taxon-specific densities with more closely associated physical data using time stamp

# Created by: Kelly Robinson
# Created on: 5 August 2015
# Date last modified: 6 August 2015
#
#---------------------------------------------------------------------------

library("plyr")
library("stringr")
library("reshape2")
library("pastecs")
library("lubridate")

load("all_joined_1E.robj")
load("all_joined_1W.robj")

e <- as.data.frame(bio_1e)
w <- as.data.frame(bio_1w)

e$dateTime <- format(e$dateTime, format = "%m-%d-%y %H:%M:%OS", tz="GMT")
w$dateTime <- format(w$dateTime, format = "%m-%d-%y %H:%M:%OS", tz="GMT")

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

##{ Read and reformat physical data ---------------------------------------

# This section performs the following tasks:
# 1) read in the physical data files from a folder
# 2) cleans up the names by removing the extra characters, units and making everything lowercase
# 3) from the time format in the physical data, converts it to the universal time format
# 4) corrects the time zone and adds in transect number, and converts lat/long into decimal degrees

# list all the physical data files in a given directory
phyFiles <- list.files("raw_physical_data_REU", full=TRUE)

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
  d$dateTime <- as.POSIXct(strptime(d$dateTime, format="%m/%d/%y %H:%M:%OS", tz="GMT"))
#   d$dateTime <- as.POSIXct(strptime(d$dateTime, format="%m/%d/%y %H:%M:%OS", tz="America/New_York"))
  # Set time zone of collection location
  # NB: we say it is America/New York when it is in fact local time, just to avoid having to deal with time zones
  
  # shift time of all physical data if not collected in Eastern Time Zone by appropriate number of hourse
  # if collected in the Eastern Time Zone, then set change value to 0
  d$dateTime <- d$dateTime - 0 * 3600  
  
  # code in a transect number
  # this is not robust for all physical data but is necessary here
  # use the file name as a dummy variable for transect number
  d$transect <- file
  #d$transect <- dd-14
  
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
                   "irrandiance"="irradiance"
  ))
  
  # keep only interesting data
  d <- d[,c("transect", "dateTime", "depth", "lat", "long", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "vol.imaged")]
  
  return(d)
  head
}, .progress="text")

# remove adply crap
phy <- phy[,-1]

# }

#KR edits to get data to estimate average volume imaged
#calcuate average number of seconds to image 1 m3

  #pull out horizontal velocity from physical data
  h <- phy[,c("dateTime", "horizontal.vel")]
  
  #convert horizontal velocity in (mm/s) to (cms/s)
  h$vel.cms <- h$horizontal.vel/10
  h <- h[ which(h$vel.cms > 0),]

  #Keep horizontal velocity in cm/s
  h <- h[,c("dateTime", "vel.cms")]

  #set constants
  image.width <- 2048
  scan.rate <- 35000
  dof <- 50
  fov <- 13.5
  #calculate px/image
  px.image <- scan.rate/image.width
  
  #calculate distanced imaged per second
  h$dist.image <- h$vel.cms/px.image
  #calculate volume imaged in cm3.
  h$vol.image <- dof*fov*h$dist.image
  #calculate volume rate(L/s)
  h$vol.rate <- (px.image*h$vol.image)/1000
  #1m3/sec
  h$m3.sec <- h$vol.rate/1000
  #sec to image 1-m3
  h$sec.image.m3 <- 1/h$m3.sec
  
  print(mean(h$sec.image.m3)) #this value should be 25,5803 sec


#bin taxon counts into 1-m3 bins using time it takes to image 1-m3 (25.5803 sec)

  e$dateTime <- as.POSIXct(e$dateTime, format = "%m-%d-%y %H:%M:%OS", tz="GMT")
  w$dateTime <- as.POSIXct(w$dateTime, format = "%m-%d-%y %H:%M:%OS", tz="GMT")

  binSize <- '25.5803 sec'
  
  #calculate the number of bins and label them in each 'e' and 'w' data frame
  e.maxT <- max(e$dateTime, na.rm=T)
  e.minT <- min(e$dateTime, na.rm=T)
  e.bins=seq(e.minT, e.maxT, by=binSize)
  e$dateTimeB <- cut(e$dateTime, breaks=e.bins, labels=1:(length(e.bins)-1))

  w.maxT <- max(w$dateTime, na.rm=T)
  w.minT <- min(w$dateTime, na.rm=T)
  w.bins=seq(w.minT, w.maxT, by=binSize)
  w$dateTimeB <- cut(w$dateTime, breaks=w.bins, labels=1:(length(w.bins)-1))

  # convert from factor to numeric
  e$dateTimeB <- as.numeric(e$dateTimeB)
  w$dateTimeB <- as.numeric(w$dateTimeB)
  
  #1E binning
  e.count.bin <- aggregate(count~dateTimeB+taxon, data = e, FUN = sum)
  e.dateTimeBin <- aggregate(dateTime~dateTimeB+taxon, data = e, FUN = mean)
  e.dateTimeBin$dateTime <- e.dateTimeBin$dateTime + 7*3600
  
  e.taxa <- merge(e.dateTimeBin, e.count.bin, by = c("dateTimeB", "taxon"))
  
  e.taxa$dateTime <- round_date(e.taxa$dateTime, "second")

  e.taxa$dateTime <- format(e.taxa$dateTime, format = "%m-%d-%y %H:%M:%OS", tz="GMT")
  e.taxa$dateTime <- as.POSIXct(e.taxa$dateTime, format = "%m-%d-%y %H:%M:%OS", tz="GMT")
  e.taxa$dateTime <- e.taxa$dateTime - 7*3600

  #1W binning
  #bin animals (total number of taxa in a bin)
  w.count.bin <- aggregate(count~dateTimeB+taxon, data = w, FUN = sum)
  
  #estimate average dateTime for each bin and taxon combination
  w.dateTimeBin <- aggregate(dateTime~dateTimeB+taxon, data = w, FUN = mean)
  
  #add 7 hours to get back to GMT because computer kept setting the tz for dateTime to be Pacific Standard Time
  w.dateTimeBin$dateTime <- w.dateTimeBin$dateTime + 7*3600
  
  #associate the average dateTime for a given bin with the total count per taxon in that bin
  w.taxa <- merge(w.dateTimeBin, w.count.bin, by = c("dateTimeB", "taxon"))
  
  w.taxa$dateTime <- round_date(w.taxa$dateTime, "second")
  
  w.taxa$dateTime <- format(w.taxa$dateTime, format = "%m-%d-%y %H:%M:%OS", tz="GMT")
  w.taxa$dateTime <- as.POSIXct(w.taxa$dateTime, format = "%m-%d-%y %H:%M:%OS", tz="GMT")
  w.taxa$dateTime <- w.taxa$dateTime - 7*3600

  #Function to try to integrate later
  #   closest.time <- function(dateTime, e.taxa){
  #   return(which.min(abs(difftime(e.taxa$datetime, e.phy$dateTime, units="secs"))))
  #   }   

#round "dateTime" to the nearest second for all frames
# h$dateTime <- round_date(h$dateTime, "second")
e$dateTime <- round_date(e$dateTime, "second")
w$dateTime <- round_date(w$dateTime, "second")

#bin physical data as a mean per second
  #Get the physical data associated with each second
  e.phy <- e[ ,c("dateTime", "Pressure.dbar", "Depth.m", "Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm", "lat","long")]
  w.phy <- w[ ,c("dateTime", "Pressure.dbar", "Depth.m", "Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm", "lat","long")]

  #physical data are averaged
  e.physec <- aggregate(cbind(Pressure.dbar, Depth.m, Temp.C, Salinity.PPT, Fluoro.volts, Oxygen.ml.l, Irrandiance.UE.cm, lat, long)~dateTime, data = e.phy, FUN = mean)
  w.physec <- aggregate(cbind(Pressure.dbar, Depth.m, Temp.C, Salinity.PPT, Fluoro.volts, Oxygen.ml.l, Irrandiance.UE.cm, lat, long)~dateTime, data = w.phy, FUN = mean)


#Merge the binned biological data with the physical data, matching by the "dateTime" fields and keeping all the biological records
e.sec <- merge(e.taxa, e.physec, by = "dateTime", all.x = T)
w.sec <- merge(w.taxa, w.physec, by = "dateTime", all.x = T)


#Add column with density of animals per 1-m3. Equals count because count is the total number of animals imaged in
#25 seconds, the average time it took the instrument to sample 1-m3
e.sec$density.1m3 <- e.sec$count
w.sec$density.1m3 <- w.sec$count

#check data
summary(e.sec)
summary(w.sec)

#save data
save(e.sec, file="all_joined_1E_density.robj")
save(w.sec, file="all_joined_1W_density.robj")
