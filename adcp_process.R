library("chron")
library("lattice")
library("RNetCDF")
library("plyr")
library("lubridate")
library("stringr")

options("digits.secs" = 3)

#functions
#-----
adcp_bin <- function(u,v,z,yxt){ #used for ADCP_test (currently stored in bin folder)
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

adcp <- function(u,v,z,lat,lon,time,flag){
  for (i in 1:length(u)){
    u_temp <- u[,i]
    v_temp <- v[,i]
    z_temp <- z[,i]
    flag_temp <- flag[,i]
    uvz <- cbind(z_temp, u_temp, v_temp, flag_temp)
    
    lat_temp <- lat[i,]
    lon_temp <- lon [i,]
    doy <- time[i,]
    
    c <- cbind(doy, lat_temp, lon_temp, uvz)
    c <- as.data.frame(c)
    c <- c[order(z_temp), ]
    
    #rename columns
    names(c)[names(c)=="lat_temp"] <- "lat"
    names(c)[names(c)=="lon_temp"] <- "lon"
    names(c)[names(c)=="z_temp"] <- "depth[m]"
    names(c)[names(c)=="u_temp"] <- "u.velocity[m/s]"
    names(c)[names(c)=="v_temp"] <- "v.velocity[m/s]"
    names(c)[names(c)=="flag_temp"] <- "qual.flag"
    
    write.csv(c, file=paste0("temp/c",i,".csv"), row.names = FALSE)
  }
}

#--------
#2014 data
#--------
setwd("C:/Users/kelly.robinson/Dropbox/Cowen_Sponaugle_share/OSTRICH/ISIIS data/OSTRICH_physical/ADCP_2014/")

#create temporary directory
dir.create("temp")
#create dictory for processed files
dir.create("ship_currents")

#--------
#2015 data
#--------
setwd("C:/Users/kelly.robinson/Dropbox/Cowen_Sponaugle_share/OSTRICH/ISIIS data/OSTRICH_physical/ADCP_2015/")

#create temporary directory
dir.create("temp")
#create dictory for processed files
dir.create("ship_currents")

#----
#code
#----

#list sensor netCDF files
f <- list.files(pattern = ".nc")

#loop through sensor netCDF files
for (i in 1:length(f)){
  fname <- f[i]
  
  fid <- open.nc(fname)
  
  #read in data fields from netCDF file
  dat<-read.nc(fid)
    
    u <- as.data.frame(dat$u)
    v <- as.data.frame(dat$v)
    lat <- as.data.frame(dat$lat)
    lon <- as.data.frame(dat$lon)
    time <- as.data.frame(dat$time)
    z <- as.data.frame(dat$depth)
    flag <- as.data.frame(dat$pflag)
  
  #run function to create individual temp files for current velocities at each lat, lon, and time observation
  adcp(u,v,z,lat,lon,time,flag)

  files <- list.files("temp", full=T)
  
  adcp_data <- adply(files, 1, function(file) {
    
    d <- read.table(file, header=T, sep=",")
   
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
    d$current_speed <- sqrt(d$u.velocity.m.s^2+d$v.velocity.m.s.^2)
  
    #calculate rotation
    d$rotation <- (180/3.14)*atan2(d$u.velocity.m.s, d$v.velocity.m.s.)
  
    #change sign of depth bin for ArcGIS
    d$depth.m. <- d$depth.m.*-1
    
    #sort data by day of year (doy)
    sort(d$doy)
    
  return(d)
  
  }, .progress="text")
  
  # remove adply crap
  adcp_data <- adcp_data[,-1]
  
  #remove obserations with non-zero quality flags
  adcp_data_good <- subset(adcp_data, qual.flag == 0)
  
  name <- str_split_fixed(string = f[i], pattern = ".nc", n = 2)
  filename <- substring(name[,1],1)
  
  write.table(adcp_data_good, file=paste0("ship_currents/",filename,".txt"), sep="\t", row.names = F, col.names = T)
}

#clean-up: delete temporary directory
unlink("temp", recursive = T)

            