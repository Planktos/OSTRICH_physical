library("data.table")
library("dplyr")

#functions
#Compute the straight line distance (km) from the starting point of a lat,lon trajectory
dist.from.start <- function(lat, lon) {
  library("oce")
  geodDist(lat1=lat, lon1=lon, lat2=na.omit(lat)[1], lon2=na.omit(lon)[1]) # use geodDist(lat1=lat, lon1=lon, lat2=na.omit(lat)[1]/1.852 if you want to convert from km to nautical miles
}


#2015 data
#-----

## os75bb sensor @18.5m
##-------
c <- fread("os75bb_currents18m_isiis_2015.txt", sep = ",", header = T, stringsAsFactors = F)
c <- as.data.frame(c)

c$distanceFromStart <- dist.from.start(c$lat, c$lon)
c$depth_m_ <- c$depth_m_*-1

t <- group_by(c, transect_id)
cvel <- summarise(t, avg_current_speed = mean(current_speed), sd = sd(current_speed))
ncvel <- cvel[order(cvel$avg_current_speed),]

t0 <- str_split_fixed(ncvel$transect_id, "-", 4)
ncvel$study <- t0[,2]
ncvel$transect <- t0[,3]

s <- subset(ncvel, study == "Eddy")
s <- s[order(s$avg_current_speed),]

g <- ggplot(s, aes(x = reorder(transect, avg_current_speed), y = avg_current_speed)) + geom_bar(stat = "identity", fill = "grey55") +
  geom_errorbar(aes(ymax = avg_current_speed + sd, ymin=avg_current_speed - sd), width=0.5) +
  scale_y_continuous("18.5m: mean current velocity (m/s) ± 1 SD") +
  scale_x_discrete("transect")+
  theme_bw()+theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), axis.title = element_text(size=14), axis.title.y = element_text(vjust = 1, face = "bold"),axis.title.x = element_text(vjust = -0.5, face = "bold"))+
  ggtitle("os75bb") + theme(plot.title = element_text(lineheight=1.0, vjust = 1.0, face="bold"))
                                                                                                                                                                                                           
png(file = paste0("data/2015_eddy_currents_os75bb.png"), width = 11, height = 8, units = "in", res = 300)
plot(g)
dev.off()




