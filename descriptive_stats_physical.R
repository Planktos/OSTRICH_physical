library("plyr")
library("stringr")
library("reshape2")
library("pastecs")
library("ggplot2")
library("ggbiplot")
library("caret")
library("e1071")
library("dplyr")

#load data
load("ost15_phy_t.R")

#historgrams of 2015 data
phy <- melt(phy_t, id.vars=c("dateTime"), measure.vars=c("depth", "temp", "salinity", "density","fluoro", "oxygen", "irradiance"))
ggplot(data=phy) + geom_histogram(aes(x=value)) + facet_wrap(~variable, scales="free")

# ---------------
# check normality
# ---------------
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

  
# Transform the data -- not used since no transformations were found to be helpful 22 Sept 2015 KR

  # FUNCTIONS:
  box.cox <- function(x, parms=c(1,0)) {
    lambda <- parms[1]
    offset <- parms[2]
    if (lambda==0) log(x+offset) else ((x+offset)^lambda - 1)/lambda
  }
  threepoint <- function(x, y, ladder=c(1, 1/2, 1/3, 0, -1/2, -1)) {
    # x and y are length-three samples from a dataset.
    dx <- diff(x)
    f <- function(parms) (diff(diff(box.cox(y, parms)) / dx))^2
    fit <- nlm(f, c(1,0))
    parms <- fit$estimate #$
    lambda <- ladder[which.min(abs(parms[1] - ladder))]
    if (lambda==0) offset = 0 else {
      do <- diff(range(y))
      offset <- optimize(function(x) f(c(lambda, x)), 
                         c(max(-min(x), parms[2]-do), parms[2]+do))$minimum    
    }
    c(lambda, offset)
  }
  
  n <- dim(p)[1]
  i3 <- c(2, floor((n+1)/2), n-1)
  
  parms <- threepoint(p$depth[i3], p$temp[i3])
  y <- box.cox(p$temp, parms)
  transformed.temp1 <- 1/sqrt(p$temp)
  
  boxcox(lm(p$temp~1))
  transformed.temp2 <- p$temp^2
  ks.test(transformed.temp2, y = pnorm, alternative = "two.sided")
  
# ---------------
# subset spatial data
# ---------------
  # subset spatial data fields
  p <- phy_t[,c("depth", "temp", "salinity", "pressure", "fluoro", "oxygen", "density","region")]
  p <- na.omit(p)
  
#   #examine variance in each parameter by depth bins
#   depth_bin = 5
#   max <- max(p$depth, na.rm=T)
#   min <- min(p$depth, na.rm=T)
#   bins=seq(from = min, to = max, by=depth_bin)
#   p$depth_bin <- as.numeric(cut(p$depth, breaks=bins, labels=1:(length(bins)-1)))
#   
#   stat.desc(p_grouped)
#   p_melt <- melt(p, id.vars = c("region", "depth_bin"))
#   p_grouped <- group_by(p_melt, depth_bin, variable)
#   summarise(p_grouped, mean = mean(value), sd=sd(value))
  
  # p$depth_bin <- NULL


p_ssub <- subset(p, p$depth >= 5 & p$region != "eddy", depth:region)
p_region <- p_ssub$region
  
p.pca <- p_ssub[,c("temp", "salinity", "fluoro", "oxygen", "density")]
  
  
#   # Check for inifite values. PCA will not run with them in the data frame
#     is.finite.data.frame <- function(obj){
#       sapply(obj,FUN = function(x) all(is.finite(x)))
#     }
#     is.finite.data.frame(p)
#     ptemp_inf <- is.infinite(p$temp)
#     psal_inf <- is.infinite(p$salinity)
#     ppress_inf <- is.infinite(p$pressure)
#     pfluoro_inf <- is.infinite(p$fluoro)
#     poxy_inf <- is.infinite(p$oxygen)
#     pden_inf <- is.infinite(p$density)

# principle component analysis
# Run PCA. Set scale=TRUE to run on a correlation matrix because variables are at different units and variances.
pca <- prcomp(p.pca, center = TRUE, scale = TRUE)
print(pca)
    
#Evaluate the varianace explained by each component  
summary(pca <- prcomp(p.pca, center = TRUE, scale = TRUE))

#correlations between PCs and variables
loadings <- pca$rotation
scores <- pca$x
correlations <- cor(scores,p.pca)
print(correlations)

  #Scree place
  plot(pca, type = "lines")

#plot Spatial PCA
g <- ggbiplot(pca, choices = 1:2, obs.scale = 1,var.scale = 1,groups = p_region, ellipse = TRUE, circle = TRUE)
g <- g + scale_color_discrete(name ="") + theme(strip.text.y = element_text(size = 12)) +
     theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 12))+
     ggtitle("2015 Spatial: PCA") + theme(plot.title = element_text(lineheight=1.5, face="bold"))

png(file = paste0("data/spatial_pca.png"), width = 8.5, height = 14, units = "in", res = 300)
plot(g)
dev.off()

# ---------------
# subset spatial data fields with only east and west observations
# ---------------
p <- phy_t[,c("depth", "temp", "salinity", "pressure", "fluoro", "oxygen", "density","area", "region")]
p <- na.omit(p)
ps <- subset(p, p$depth >= 5 & p$region != "eddy" & p$region != "c", depth:region)
pew <- subset(ps, ps$area == "6"| ps$area == "7", depth:region)

p_region <- pew$region

pew.pca <- pew[,c("temp", "salinity", "fluoro", "oxygen", "density")]

# principle component analysis
# Run PCA. Set scale=TRUE to run on a correlation matrix because variables are at different units and variances.
pca <- prcomp(pew.pca, center = TRUE, scale = TRUE)
print(pca)

#Evaluate the varianace explained by each component  
summary(pca <- prcomp(pew.pca, center = TRUE, scale = TRUE))

#Scree place
plot(pca, type = "lines")

#correlations between PCs and variables
loadings <- pca$rotation
scores <- pca$x
correlations <- cor(scores,pew.pca)
print(correlations)

#plot Spatial PCA
g <- ggbiplot(pca, choices = 4:5, obs.scale = 1,var.scale = 1,groups = p_region, ellipse = TRUE, circle = TRUE)
g <- g + scale_color_discrete(name ="") + theme(strip.text.y = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 12))+
  ggtitle("PCA score plot: 2015 Spatial E-W ") + theme(plot.title = element_text(lineheight=1.5, face="bold"))

png(file = paste0("data/spatial_EW_PC4-PC5.png"), width = 8.5, height = 14, units = "in", res = 300)
plot(g)
dev.off()

