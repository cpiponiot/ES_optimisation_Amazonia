
library(data.table)
library(rgdal)
library(automap)
library(raster)
library(ggplot2)
# library(foreign)
# library(BIOMASS)
# library(pgirmess)

load("C:/Users/camille.piponiot/Google Drive/radam/volumeRadam.Rdata")

##############################################################################
##################                 KRIGING                  ##################
##############################################################################

dsn="C:/Users/camille.piponiot/Google Drive/maps/amazonia"
amazonia <- readOGR(dsn=dsn, layer="BassinAmazonien")

### problem: we have repeated coordinates ###
dupl = unique(vol_data[duplicated(cbind(vol_data$long,vol_data$lat)),c("long","lat")])
##average volume values on those plots ###
dupl2 = apply(dupl, 1, function(coords) {
  c(vol_data$plot[vol_data$long==coords[1] & vol_data$lat==coords[2]][1],
    coords,
    mean(vol_data$vol[vol_data$long==coords[1] & vol_data$lat==coords[2]]),
    mean(vol_data$om0[vol_data$long==coords[1] & vol_data$lat==coords[2]]), 
    mean(vol_data$WDcom[vol_data$long==coords[1] & vol_data$lat==coords[2]]))
})
dupl2 = matrix(unlist(dupl2), ncol=6, byrow=T)
dupl2 = data.frame(dupl2)
colnames(dupl2) = colnames(vol_data)

not_dupl = vol_data[!(duplicated(cbind(vol_data$long,vol_data$lat)))]
vol_data = rbind(not_dupl, dupl2)
vol_data = subset(vol_data, !is.na(WDcom))

# vol_data -> spatial object
coordinates(vol_data) = ~long+lat
# reproject to krige
proj4string(vol_data) <- CRS("+proj=longlat +datum=WGS84")
vol_data <- spTransform(vol_data, proj4string(amazonia))
plot(vol_data)


# create grid
grd <- expand.grid(x=seq(from=-3e6, to=1.5e6, by=1e5), y=seq(from=-2.1e6, to=1.2e6, by=1e5))
coordinates(grd) <- ~ x+y
gridded(grd) <- TRUE
proj4string(grd) <- proj4string(amazonia)
grd = grd[!is.na(over(spTransform(grd,proj4string(amazonia)), amazonia))]
plot(grd)

# choose variogram model
variog <- autofitVariogram(formula = WDcom~1, input_data = vol_data)
plot(variog)

# WDext
krig_wd <- autoKrige(WDcom~1, vol_data, new_data=grd)
krig.output_wd <- krig_wd$krige_output[1:3]
names(krig.output_wd) <- c("x","y","WDcom")
plot(krig.output_wd)
borders <- readOGR(dsn = "C:/Users/camille.piponiot/Google Drive/maps/borders",layer="TM_WORLD_BORDERS-0.3")
brasil <- spTransform(borders[borders$ISO3=="BRA",],proj4string(krig.output_wd))
output_wd <- raster(krig.output_wd)
output_wd <- projectRaster(output_wd,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
names(output_wd) <- "wd"

writeRaster(output_wd,file="maps/map_WDcom_radam.tif", overwrite=TRUE)
