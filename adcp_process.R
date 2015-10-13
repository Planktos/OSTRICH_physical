#process ADCP data

library(stringr)
library(plyr)
library(dplyr)


options("digits.secs" = 3)

u <- read.csv("u_test.csv",stringsAsFactors = F, header = F)
v <- read.csv("v_test.csv",stringsAsFactors = F, header = F)
xyt <- read.csv("current_xyt.csv",stringsAsFactors = F, header = F)
z <- read.csv("depth_adcp_test.csv",stringsAsFactors = F, header = F)

dir.create("temp")

#function-----
adcp <- function(u,v,z,yxt){
  for (i in 1:length(u)){
    u_temp <- u[,i]
    v_temp <- v[,i]
    uv <- cbind(u_temp, v_temp,z)
    
    xyt1 <- xyt[,c(i)]
    lon <- xyt1[1:1]
    lat <- xyt1[2:2]
    doy <- xyt1[3:3]
    d <- data.frame(uv,lat,lon,doy)
    write.csv(d,file=paste0("temp/d",i,".csv"),row.names = FALSE)
  }
}
#------
#run function
adcp(u,v,z,xyt)

setwd("C:/Users/kelly.robinson/Dropbox/Cowen_Sponaugle_share/OSTRICH/ISIIS data/OSTRICH_physical/temp")
files <- list.files("C:/Users/kelly.robinson/Dropbox/Cowen_Sponaugle_share/OSTRICH/ISIIS data/OSTRICH_physical/temp")


adcp_data <- adply(files, 1, function(file) {
  t <- read.table(file, header=T, sep=",")
  #name the fields
  colnames(t) <- c("u","v","z_bin","lat","lon","doy")
    return(t)
}, .progress="text")

setwd("C:/Users/kelly.robinson/Dropbox/Cowen_Sponaugle_share/OSTRICH/ISIIS data/OSTRICH_physical")

# remove adply crap
adcp_data <- adcp_data[,-1]

d<-adcp_data

  library(lubridate)

  #create proper date
  d$date<-strptime(x = d$doy, format ="%j", tz="GMT")
  d$date<-as.POSIXct(d$date)
  year <- year(d$date)
  month <- month(d$date)
  day <- day(d$date)
  
  #create proper time
  t <- trunc(d$doy,3)
  dec.day <- d$doy-t
  dec.hh <- dec.day*24
  hour <- as.integer(dec.hh)
  dec.min <- dec.hh-hour
  min <- as.integer(dec.min*60)
  dec.sec <- (dec.min*60)-min
  sec <- dec.sec*60
  
  #combine date and time
  d$dateTime<-ISOdatetime(year,month,day,hour,min,sec, tz="GMT")
 
  #create time field in format yyyymmddhhminss.s for ArcMap
  d$dateTime <- as.character(d$dateTime)
  t <- str_split_fixed(d$dateTime, "-|:| ",6)
  yy <- substring(t[,1],1)
  mm <- substring(t[,2],1)
  dd <- substring(t[,3],1)
  hh <- substring(t[,4],1)
  min <- substring(t[,5],1)
  sec <- substring(t[,6],1)
  d$time <- as.numeric(paste0(yy,mm,dd,hh,min,sec))
  
  #calculate current speed
  d$current_speed <- sqrt(d$u^2+d$v^2)
  
  #calculate rotation
  d$rotation <- (180/3.14)*atan2(d$u,d$v)
  
  #change sign of depth bin for ArcGIS
  d$z_bin <- d$z_bin*-1
  
  #sort data by day of year (doy)
  sort(d$doy)
  
  d <- na.omit(d)
  
  write.table(d, file="ost15_WS_currents.txt", sep="\t", row.names = F, col.names = T)
  







