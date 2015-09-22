library("lubridate")
library("plyr")
library("stringr")
library("reshape2")
library("pastecs")
library("ggplot2")
library("oce")

#load data
load("ost15_phy_t.R")

#historgrams of 2015 data
phy <- melt(phy_t, id.vars=c("dateTime"), measure.vars=c("depth", "temp", "salinity", "density","fluoro", "oxygen", "irradiance"))
ggplot(data=phy) + geom_histogram(aes(x=value)) + facet_wrap(~variable, scales="free")

# check normality
##temp
  qqnorm(phy_t$temp)
  ks.test(x = phy_t$temp,y = pnorm, alternative = "two.sided")
  
  phy_t$trans.temp <- log10(phy_t$temp)
  qqnorm(phy_t$trans.temp)
  ks.test(x = phy_t$trans.temp,y = pnorm, alternative = "two.sided")

##salinity
  qqnorm(phy_t$salinity)
  ks.test(x = phy_t$salinity,y = pnorm, alternative = "two.sided")
  
  phy_t$trans.salinity <- log10(phy_t$salinity)
  qqnorm(phy_t$trans.salinity)
  ks.test(x = phy_t$trans.salinity,y = pnorm, alternative = "two.sided")
  
##density
  qqnorm(phy_t$density)
  ks.test(x = phy_t$density,y = pnorm, alternative = "two.sided")
  
  phy_t$trans.density <- log10(phy_t$density)
  qqnorm(phy_t$trans.density)
  ks.test(x = phy_t$trans.density,y = pnorm, alternative = "two.sided")
  
##fluoro
  qqnorm(phy_t$fluoro)
  ks.test(x = phy_t$fluoro,y = pnorm, alternative = "two.sided")
  
  phy_t$trans.fluoro <- log10(phy_t$fluoro)
  qqnorm(phy_t$trans.fluoro)
  ks.test(x = phy_t$trans.fluoro,y = pnorm, alternative = "two.sided")

##oxygen
  qqnorm(phy_t$oxygen)
  ks.test(x = phy_t$oxygen,y = pnorm, alternative = "two.sided")
  
  phy_t$trans.oxygen <- log10(phy_t$oxygen)
  qqnorm(phy_t$trans.oxygen)
  ks.test(x = phy_t$trans.oxygen,y = pnorm, alternative = "two.sided")
  
##irradiance
  qqnorm(phy_t$irradiance)
  ks.test(x = phy_t$irradiance,y = pnorm, alternative = "two.sided")
  
  phy_t$trans.irradiance <- log10(phy_t$irradiance)
  qqnorm(phy_t$trans.irradiance)
  ks.test(x = phy_t$trans.irradiance,y = pnorm, alternative = "two.sided")

phyt <- phyt[,c("dateTime", "depth", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "density", "haul")]
