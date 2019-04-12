
##### Features ##### 
## PPF area: harv = currently available, harvAR = all unprotected areas
harv <- crop(raster(grd,"pAreaAvail")*raster(grd,"area"), extent(cost_carbon[[1]]))
harvAR <- crop(raster(grd,"pAreaUnprotect")*raster(grd,"area"), extent(cost_carbon[[1]])) 


# Intact forest landscape #
feature_ifl = sapply(1:nrow(df_zones), function(i) {
  sapply(unique(grd$region), function(reg) {
    grd$ifl_zr = (grd$area*grd$pIFL)*(grd$region==reg)*(1-grd$pAreaAvail*(df_zones$zname[i]!="NL"))
    crop(raster(grd,"ifl_zr"), extent(cost_carbon))
  })
})

feature_iflAR = sapply(1:nrow(df_zones), function(i) {
  sapply(unique(grd$region), function(reg) {
    grd$ifl_zr = (grd$area*grd$pIFL)*(grd$region==reg)*(1-grd$pAreaUnprotect*(df_zones$zname[i]!="NL"))
    crop(raster(grd,"ifl_zr"), extent(cost_carbon))
  })
})

ones <- feat_prod[[1]][[1]]
ones[!is.na(ones),] <- 1

### base scenario
features <- sapply(1:nrow(df_zones), function(i) stack(feat_prod[[i]]*harv,
                                                       feature_ifl[,i]$CA,
                                                       feature_ifl[,i]$EA,
                                                       feature_ifl[,i]$GS,
                                                       feature_ifl[,i]$NWA,
                                                       feature_ifl[,i]$SEA,
                                                       feature_ifl[,i]$SWA, ones))
z0 <- zones("LS" = features[[1]], "LM" = features[[2]], "LL" = features[[3]], 
            "MS" = features[[4]], "MM" = features[[5]], "ML" = features[[6]],
            "HS" = features[[7]], "HM" = features[[8]], "HL" = features[[9]], 
            "NL" = features[[10]],
            feature_names = c("timber_extr", "IFL_CA", "IFL_EA", "IFL_GS", 
                              "IFL_NWA", "IFL_SEA", "IFL_SWA","allIncluded"))
### no short cutting cycle
zMed <- z0
## remove 'allIncluded' feature where trot != 30 (or 1 : no logging)
for(i in which(!(df_zones$trot %in% c(1,30)))) 
  zMed[[i]][[dim(zMed[[i]])[3]]][!is.na(ones),] <- 0 

### STY scenario
## avoid negative features for timber recovery:
offset_vrec <- max(values(cost_vrec*harvAR), na.rm=TRUE)
## the total recovery when there is no logging (ie target for STY):
zero_vcost <- offset_vrec*nrow(grd)

featuresSTY <- sapply(1:nrow(df_zones), function(i) stack(feat_prod[[i]]*harv,
                                                          -1*cost_vrec[[i]]*harv + offset_vrec,
                                                          feature_ifl[,i]$CA,
                                                          feature_ifl[,i]$EA,
                                                          feature_ifl[,i]$GS,
                                                          feature_ifl[,i]$NWA,
                                                          feature_ifl[,i]$SEA,
                                                          feature_ifl[,i]$SWA, ones) )
zSTY <- zones("LS" = featuresSTY[[1]], "LM" = featuresSTY[[2]], "LL" = featuresSTY[[3]], 
              "MS" = featuresSTY[[4]], "MM" = featuresSTY[[5]], "ML" = featuresSTY[[6]],
              "HS" = featuresSTY[[7]], "HM" = featuresSTY[[8]], "HL" = featuresSTY[[9]], 
              "NL" = featuresSTY[[10]],
              feature_names = c("timber_extr", "timber_rec",
                                "IFL_CA", "IFL_EA", "IFL_GS", 
                                "IFL_NWA", "IFL_SEA", "IFL_SWA","allIncluded"))

### land sharing scenario
featuresAR <- sapply(1:nrow(df_zones), function(i) stack(feat_prod[[i]]*harvAR, 
                                                         feature_iflAR[,i]$CA,
                                                         feature_iflAR[,i]$EA,
                                                         feature_iflAR[,i]$GS,
                                                         feature_iflAR[,i]$NWA,
                                                         feature_iflAR[,i]$SEA,
                                                         feature_iflAR[,i]$SWA, ones) )
zAR <- zones("LS" = featuresAR[[1]], "LM" = featuresAR[[2]], "LL" = featuresAR[[3]], 
             "MS" = featuresAR[[4]], "MM" = featuresAR[[5]], "ML" = featuresAR[[6]],
             "HS" = featuresAR[[7]], "HM" = featuresAR[[8]], "HL" = featuresAR[[9]], 
             "NL" = featuresAR[[10]],
             feature_names = c("timber_extr", 
                               "IFL_CA", "IFL_EA", "IFL_GS", 
                               "IFL_NWA", "IFL_SEA", "IFL_SWA","allIncluded"))

## STY and land sharing

features_AR_STY <- sapply(1:nrow(df_zones), function(i) stack(feat_prod[[i]]*harvAR,
                                                              -1*cost_vrec[[i]]*harvAR + offset_vrec,
                                                              feature_iflAR[,i]$CA,
                                                              feature_iflAR[,i]$EA,
                                                              feature_iflAR[,i]$GS,
                                                              feature_iflAR[,i]$NWA,
                                                              feature_iflAR[,i]$SEA,
                                                              feature_iflAR[,i]$SW, ones) )

zAR_STY <- zones("LS" = features_AR_STY[[1]], "LM" = features_AR_STY[[2]], "LL" = features_AR_STY[[3]], 
                 "MS" = features_AR_STY[[4]], "MM" = features_AR_STY[[5]], "ML" = features_AR_STY[[6]],
                 "HS" = features_AR_STY[[7]], "HM" = features_AR_STY[[8]], "HL" = features_AR_STY[[9]], 
                 "NL" = features_AR_STY[[10]],
                 feature_names = c("timber_extr", "timber_rec",
                                   "IFL_CA", "IFL_EA", "IFL_GS", 
                                   "IFL_NWA", "IFL_SEA", "IFL_SWA","allIncluded"))

featureMaps = list(z0,zAR,zSTY,zAR_STY,zMed)
names(featureMaps) = c("base","sharing","STY","sharingSTY","medium")

##### Costs ##### 
# 
# minV <- min(values(cost_vrec), na.rm = TRUE)
# minC <- min(values(cost_carbon), na.rm = TRUE)
# minB <- min(values(cost_diversity), na.rm = TRUE)

# puV <- (cost_vrec - minV) / mean(values(cost_vrec - minV), na.rm=T)
# puC <- (cost_carbon - minC) / mean(values(cost_carbon - minC), na.rm=T)
# puB <- (cost_diversity - minB) / mean(values(cost_diversity - minB), na.rm=T)
puV <- (cost_vrec - mean(values(cost_vrec), na.rm=T)) / sd(values(cost_vrec), na.rm=T)
puC <- (cost_carbon - mean(values(cost_carbon), na.rm=T)) / sd(values(cost_carbon), na.rm=T)
puB <- (cost_diversity - mean(values(cost_diversity), na.rm=T)) / sd(values(cost_diversity), na.rm=T)
puBal <- puV * coeffs_balanced[1] + puC * coeffs_balanced[2] + puB * coeffs_balanced[3]

costMaps = list(timber = puV * harv, 
                carbon = puC * harv, 
                biodiversity = puB * harv, 
                balanced = puBal * harv, 
                sharing = puBal * harvAR)

areaIFL = tapply(grd$pIFL*grd$area,grd$region,sum)
Nfeature = dim(features[[1]])[3]
NfeatureSTY = dim(featuresSTY[[1]])[3]

#####  Targets ##### 

### base scenario
t0 <- tibble::tibble(feature = c("timber_extr", paste("IFL",names(areaIFL), sep="_"),"allIncluded"), 
                     zone = list(df_zones$zname)[rep(1,Nfeature)],
                     target = c(35e6, 0.8*areaIFL,nrow(grd)), 
                     type=rep("absolute", Nfeature))
### sty scenario
tSTY <- tibble::tibble(feature = c("timber_extr", "timber_rec", 
                                   paste("IFL",names(areaIFL), sep="_"),"allIncluded"), 
                       zone = list(df_zones$zname)[rep(1,NfeatureSTY)],
                       target = c(35e6, zero_vcost, 0.8*areaIFL,nrow(grd)), 
                       type=rep("absolute", NfeatureSTY))

targetsList = list(t0,tSTY)
names(targetsList) = c("base","STY")

