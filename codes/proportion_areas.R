area_FC = area(forest_cover_90)
area_FC[is.na(forest_cover_90), ] <- NA
areaFC = sum(values(area_FC), na.rm=T)

protected_forest = protected
protected_forest[is.na(forest_cover_90),] <- NA
area_protected = area(protected_forest)
area_protected[is.na(protected_forest),] <- NA
areaProtected = sum(values(area_protected), na.rm=T)

area_avail = area(current_available)
area_avail[is.na(current_available),] <- NA
areaAvail = sum(values(area_avail), na.rm=T)

Pwdpa = round(areaProtected/areaFC*100)
Pharv = round(areaAvail/areaFC*100)
Punacc = 100-Pwdpa-Pharv


### calculate PPF area for each grid cell

grd$areaTot = extract(area(raster(grd)), grd) ## in km2 
grd$areaTot = grd$areaTot * 100 ## in ha
grd$area = grd$areaTot * 0.58 ## multiplication by coefficient pi (proportion of areas not suitable for logging)

# proportion of forested areas (FC > 90%)
forest_cover_90b = forest_cover_90
forest_cover_90b[is.na(forest_cover_90b),] = 0
grd$pAreaForest = extract(aggregate(forest_cover_90b, fact = 1/res(forest_cover_90b)), grd)

## proportion of currently available areas
current_availb = current_available
current_availb[is.na(current_availb),] = 0
grd$pAreaAvail = extract(aggregate(current_availb, fact = 1/res(current_availb)), grd)

## proportion of currently available areas
unprotected_forest = forest_cover_90
unprotected_forest[is.na(forest_cover_90),] = 0
unprotected_forest[!is.na(protected),] = 0
grd$pAreaUnprotect = extract(aggregate(unprotected_forest, fact = 1/res(unprotected_forest)), grd)


### diagram 
png("graphs/map_PPF1.png", height=350)
par(mar=c(0,0,0,0), oma = c(0,0,0,0))
plot(forest_cover_90, legend = FALSE, axes=FALSE, box=FALSE, col="black")
dev.off()
png("graphs/map_PPF2.png", height=350)
par(mar=c(0,0,0,0), oma = c(0,0,0,0))
unprotected_forest[unprotected_forest==0,] = NA
plot(unprotected_forest, legend = FALSE, axes=FALSE, box=FALSE, col="black")
dev.off()
png("graphs/map_PPF3.png", height=350)
par(mar=c(0,0,0,0), oma = c(0,0,0,0))
plot(current_available, legend = FALSE, axes=FALSE, box=FALSE, col="black")
dev.off()
png("graphs/map_PPF4.png", height=200, width = 320)
par(mar=c(0,0,0,0), oma = c(0,0,0,0))
plot(raster(grd, "pAreaUnprotect")*raster(grd, "area")*1e-6, box = FALSE,
     axes=FALSE, legend.args=list(text='PPF area (Mha)', side = 4, line = 2.5, cex = 1.5))
dev.off()
png("graphs/map_PPF5.png", height=200, width = 320)
par(mar=c(0,0,0,0), oma = c(0,0,0,0))
plot(raster(grd, "pAreaAvail")*raster(grd, "area")*1e-6, 
     axes=FALSE, box=FALSE, legend.args=list(text='PPF area (Mha)', side = 4, line = 2.5, cex = 1.5))
dev.off() 
png("graphs/map_PPF6.png", height=200, width = 320)
par(mar=c(0,0,0,0), oma = c(0,0,0,0))
plot(raster(grd, "pAreaForest")*raster(grd, "area")*1e-6, 
     axes=FALSE, box=FALSE, legend.args=list(text='PPF area (Mha)', side = 4, line = 2.5, cex = 1.5))
dev.off() 