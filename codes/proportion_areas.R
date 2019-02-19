area_FC = area(mask_FC)
area_FC[is.na(mask_FC),]<-NA
areaFC = sum(values(area_FC), na.rm=T)

wdpa_am_forest = wdpa_am_raster
wdpa_am_forest[is.na(mask_FC),]<-NA
area_wdpa = area(wdpa_am_forest)
area_wdpa[is.na(wdpa_am_forest),]<-NA
areaWDPA = sum(values(area_wdpa), na.rm=T)

area_harv = area(harvestableAreas)
area_harv[is.na(harvestableAreas),]<-NA
areaHarv = sum(values(area_harv), na.rm=T)


Pwdpa = round(areaWDPA/areaFC*100)
Pharv = round(areaHarv/areaFC*100)
Punacc = 100-Pwdpa-Pharv