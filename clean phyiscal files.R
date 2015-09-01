

library("plyr")
library("stringr")
library("reshape2")
library("pastecs")
library("lubridate")

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

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
  d$dateTime <- as.POSIXct(strptime(d$dateTime, format="%m/%d/%y %H:%M:%OS", tz="American/New_York"))
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

