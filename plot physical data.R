
#Plot physical data for each transect

library("plyr")
library("stringr")
library("ggplot2")
library("reshape2")
library("grid")
library("gridExtra")

load("ost14_phy.R")




vars <- c("temp", "salinity", "fluoro", "oxygen", "irradiance")

dist.from.start(d$lat, d$long)

dm <- melt(ed1, id.vars=c("dateTimeR", "depth", "distance.km"), measure.vars=vars)  


# For constant depth transects (e.g. shallow, mid, deep)

## Create plots using 'dateTime' as the x-axis value

  t <- ggplot(dm[dm$variable == "temp",], aes(x=(dateTimeR), y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 15))
  
  s <- ggplot(dm[dm$variable == "salinity",], aes(x=dateTimeR, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 15))
  
  f <- ggplot(dm[dm$variable == "fluoro",], aes(x=dateTimeR, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 15))

  o <- ggplot(dm[dm$variable == "oxygen",], aes(x=dateTimeR, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 15))

  i <- ggplot(dm[dm$variable == "irradiance",], aes(x=dateTimeR, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 15))
  
  
  # Create plots using 'distance.km' as the x-axis value
  t <- ggplot(dm[dm$variable == "temp",], aes(x=distance.km, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance.km", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 15))
  
  s <- ggplot(dm[dm$variable == "salinity",], aes(x=distance.km, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance.km", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 15))
  
  f <- ggplot(dm[dm$variable == "fluoro",], aes(x=distance.km, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance.km", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 15))

  o <- ggplot(dm[dm$variable == "oxygen",], aes(x=distance.km, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance.km", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 15))

  i <- ggplot(dm[dm$variable == "irradiance",], aes(x=distance.km, y=-depth)) +
    geom_line(aes(colour=value, size = value), na.rm=T, show_guide = FALSE) +
    scale_colour_gradient(high=spectral(), na.value=NA) +
    scale_x_continuous("distance.km", expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) + facet_grid(variable~.) +
    theme(strip.text.y = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 15))