library(sp)
library(rgdal)
library(raster)

#### create Amazonian 1deg grid #### 
amazonia <- readOGR("C:/Users/camille.piponiot/Google Drive/maps/amazonia","BassinAmazonien")
amazonia <- spTransform(amazonia, CRS("+proj=longlat"))
grd <- expand.grid(long=seq(from=-79.5, to=-44.5, by=1), lat=seq(from=-19.5, to=9.5, by=1))
coordinates(grd) <- ~ long+lat
gridded(grd) <- TRUE
proj4string(grd) <- proj4string(amazonia)
grd = grd[!is.na(over(spTransform(grd,proj4string(amazonia)), amazonia))]


####  proportion of harvestable area #### 

## within 25 km of a road
load("C:/Users/camille.piponiot/Google Drive/volume_recovery/AmazonianHarvestableAreas.Rdata")
pHarv <- aggregate(harvestableAreas, fact=1/res(harvestableAreas))
pHarv <- projectRaster(pHarv,raster(grd))
grd$pHarv <- raster::extract(pHarv, grd)
# grd$pHarv[grd$pHarv==0] <- NA

##  all unprotected areas 
load("C:/Users/camille.piponiot/Google Drive/maps/forest cover/FC.Rdata")
load("C:/Users/camille.piponiot/Google Drive/volume_recovery/WDPA_AmazonRaster.Rdata")
mask_FC_AR <- mask_FC
mask_FC_AR[wdpa_am_raster,] <- NA
mask_FC_AR[is.na(mask_FC_AR),] <- 0
pHarvAR <- aggregate(mask_FC_AR, fact=1/res(mask_FC_AR))
pHarvAR <- projectRaster(pHarvAR,raster(grd))
grd$pHarvAR <- raster::extract(pHarvAR, grd)

## total area of a pixel (km2)
grd$area <- raster::extract(area(raster(grd)), grd)

#### environmental variables ####

e = extent(-80,-44,-20,10)

prec <- raster("C:/Users/camille.piponiot/Google Drive/maps/wc2.0/wc2.0_bio_2.5m_12.tif")
prec <- aggregate(crop(prec,e), fact=24)
grd$prec <- raster::extract(prec, grd)

seas <- raster("C:/Users/camille.piponiot/Google Drive/maps/wc2.0/wc2.0_bio_2.5m_15.tif")
seas <- aggregate(crop(seas,e), fact=24)
grd$seas <- raster::extract(seas, grd)

# bkd <- raster("C:/Users/camille.piponiot/Google Drive/maps/soilgrids/bulk density -- 100cm/BLDFIE_M_sl6_250m.tif")
# bkd <- raster::aggregate(bkd, fact = 480)
# writeRaster(bkd,"maps/bkd_aggr.tif")
# bkd <- raster("maps/bkd_aggr.tif")
# grd$bkd <- raster::extract(bkd, grd)

bkd_HWSD <- raster("C:/Users/camille.piponiot/Google Drive/maps/HWSD/topBD.tif")
bkd_HWSD <- raster::aggregate(bkd_HWSD, fact = 120)
grd$bkd <- raster::extract(bkd_HWSD, grd)

## !!! make variable maps only on forested pixels! mask with forest cover map
# acs <- raster("C:/Users/camille.piponiot/Google Drive/maps/carbon/Avitabile_AGB_Map.tif")
# acs <- crop(acs, e)/2
# ## forest cover map
# load("C:/Users/camille.piponiot/Google Drive/maps/forest cover/FC.Rdata")
# mask_FC <- resample(mask_FC, acs)
# acs[is.na(mask_FC),] <- NA
# acs[acs<20,] <- NA
# acs <- aggregate(acs, fact = 120)
# writeRaster(acs,"maps/acs_aggr.tif", overwrite=TRUE)
acs <- raster("maps/acs_aggr.tif")
grd$acs <- raster::extract(acs, grd)

#### volume parameters ####
load("C:/Users/camille.piponiot/Google Drive/volume_recovery/maps.Rdata")
load("C:/Users/camille.piponiot/Google Drive/volume_recovery/map_aG_FORMIND.Rdata")
map_om0 <- raster("C:/Users/camille.piponiot/Google Drive/volume_recovery/map_omega0.tif")
grd$om0 <- raster::extract(map_om0, grd)
grd$aG <- raster::extract(aG_FORMIND, grd)
grd$vmax <- raster::extract(vmax_Am, grd)
grd$stem_mort <- raster::extract(stem_mort, grd)

#### richness ####
load("maps/richness_maps.Rdata")
grd$mammals <- raster::extract(richness$mammals.all_spp, grd)
grd$amphi <- raster::extract(richness$amphibians.all_spp, grd)

#### ter Steege ecoregions ####
ecoregions <- readOGR("C:/Users/camille.piponiot/Google Drive/maps/ecoregions","regions1")
ecoregions$ID = c("GS","SEA","CA","EA","NWA","SWA")
grd$region <- raster::extract(rasterize(ecoregions,raster(grd)), grd)
grd$region <- c("GS","SEA","CA","EA","NWA","SWA")[grd$region]

#### proportion of IFL ####
propIFL <- raster("maps/propIFL.tif")
grd$pIFL <- raster::extract(propIFL, grd)

#### commercial trees' wood density ####
grd$WDext <- raster::extract(raster("maps/map_WDcom_radam.tif"), grd)
## 1 NA value
isna = which(is.na(grd$WDext))
coordsna = coordinates(grd)[isna,]
# neighbors
neighna = which((coordinates(grd)[,2] >= coordsna[2]-1    & 
                   coordinates(grd)[,2] <= coordsna[2]+1  & 
                   coordinates(grd)[,1] >= coordsna[1]-1  & 
                   coordinates(grd)[,1] <= coordsna[1]+1 )) 
grd$WDext[isna]<- mean(grd$WDext[neighna], na.rm=TRUE)
