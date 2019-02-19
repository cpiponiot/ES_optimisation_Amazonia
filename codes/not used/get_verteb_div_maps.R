
# final resolution of maps in ° (here: 1°)
resol = 1 

#### Amphibians richness ####
amphi <- sapply(list.files(path="biodiversity maps/biodiversitymapping_TIFFs/Amphibians", pattern=".tif$",full.names = T), 
                  function(pth) {
                    mp = raster(pth)
                    mp <- projectRaster(mp, crs = CRS("+proj=longlat"))
                    r1 <- raster(matrix(1, ncol = floor(resol/res(mp)[1])*diff(e[1:2]), 
                                        nrow = floor(resol/res(mp)[2])*diff(e[3:4])))
                    extent(r1) <- e; crs(r1) <- CRS("+proj=longlat")
                    mp <- projectRaster(mp, r1)
                    
                    mp <- aggregate(mp, fact = resol/res(mp))
                  })
names(amphi) = paste("amphibians", gsub(gsub(list.files(path="biodiversity maps/biodiversitymapping_TIFFs/Amphibians", pattern=".tif$"),pattern = "_raster.tif",replacement = ""), 
                                       pattern = "richness_10km_", replacement = ""))

#### Mammals richness ####

mammals <- sapply(list.files(path="biodiversity maps/biodiversitymapping_TIFFs/Mammals", pattern=".tif$",full.names = T), 
                  function(pth) {
                    mp = raster(pth)
                    mp <- projectRaster(mp, crs = CRS("+proj=longlat"))
                    r1 <- raster(matrix(1, ncol = floor(resol/res(mp)[1])*diff(e[1:2]), 
                                        nrow = floor(resol/res(mp)[2])*diff(e[3:4])))
                    extent(r1) <- e; crs(r1) <- CRS("+proj=longlat")
                    mp <- projectRaster(mp, r1)
                    
                    mp <- aggregate(mp, fact = resol/res(mp))
                  })
names(mammals) = paste("mammals", gsub(gsub(list.files(path="biodiversity maps/biodiversitymapping_TIFFs/Mammals", pattern=".tif$"),pattern = "_raster.tif",replacement = ""), 
     pattern = "richness_10km_", replacement = ""))


#### Birds richness ####

birds <- sapply(list.files(path="biodiversity maps/biodiversitymapping_TIFFs/Birds", pattern=".tif$",full.names = T), 
                  function(pth) {
                    mp = raster(pth)
                    mp <- projectRaster(mp, crs = CRS("+proj=longlat"))
                    r1 <- raster(matrix(1, ncol = floor(resol/res(mp)[1])*diff(e[1:2]), 
                                        nrow = floor(resol/res(mp)[2])*diff(e[3:4])))
                    extent(r1) <- e; crs(r1) <- CRS("+proj=longlat")
                    mp <- projectRaster(mp, r1)
                    
                    mp <- aggregate(mp, fact = resol/res(mp))
                  })
names(birds) = paste("Birds", gsub(gsub(list.files(path="biodiversity maps/biodiversitymapping_TIFFs/Birds", pattern=".tif$"),pattern = "_raster.tif",replacement = ""), 
                                       pattern = "richness_10km_", replacement = ""))

richness = stack(c(amphi,mammals,birds))
save(richness, file = "richness_maps.Rdata")
