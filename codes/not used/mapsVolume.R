source("C:/Users/camille.piponiot/Google Drive/volume_recovery/codes/omega_t.R")
source("C:/Users/camille.piponiot/Google Drive/volume_recovery/codes/modelPredictions.R")

library(truncnorm)

pred$t0 <- (100/pred$stem_mort) ^ pred$lambda_ti * (1-pred$dti)
pred$aM <- apply(cbind(0, pred$aG - pred$theta*pred$vmax),1,max)

pred$V0 = rtruncnorm(n = nrow(pred), mean = volume(pred$t0, pred$aG, pred$aM, pred$bG, pred$bM, pred$theta, pred$pdef), sd = pred$sigmaV, a = 0)
pred$mu_V0 = volume(pred$t0, pred$aG, pred$aM, pred$bG, pred$bM, pred$theta, pred$pdef)
pred$Vcom0 = pred$V0*pred$om0
pred$vextReal = apply(pred[,c("vext","Vcom0")],1,min)

# deltaV: total volume loss
pred$deltaV = apply(cbind(pred$vext/(pred$om0^(1-pred$rho)),pred$mu_V0), 1, min)
dt0 = pred[,.(t1=t0Prediction(t0, deltaV, aG,aM,bG,bM,theta,pdef)), .(long,lat,vext,trot)]
pred = merge(pred, dt0, by=c("long","lat","vext","trot"))
# post-logging proprotion of commercial species: om1
pred$om1 = (pred$om0*pred$V0-pred$vextReal)/(pred$V0-pred$deltaV)
# proportion of commercial species at the end of the cutting cycle
DTom2 = pred[,.(om2 = omega_t(tvec=c(t1+trot),t1,om1,om0,aG,aM,bG,bM,theta,intercept, slope)),.(vext,trot,long,lat)]
pred = merge(pred, DTom2, by=c("vext","trot","long","lat"))

# recovered timber volume
pred$vrec = volume(t=pred$t1+pred$trot, ag=pred$aG, am=pred$aM,
                   bg=pred$bG, bm=pred$bM,th=pred$theta, pdef=0.2)*pred$om2 - (pred$V0 - pred$deltaV)*pred$om1

# extracted volume maps (ie production)
MvextReal = dcast(pred, long + lat ~ zname , value.var = "vextReal")
MvextReal$NL = 0
coordinates(MvextReal) <- ~ long+lat
gridded(MvextReal) <- TRUE
proj4string(MvextReal) <- proj4string(grd)

feat_prod <- sapply(1:nrow(df_zones), function(i) {
  raster(MvextReal,as.character(df_zones$zname[i]))/df_zones$trot[i]
})

# commercial volume maps 
Mvcom0 = subset(pred, zone==1 )[,c("long","lat","Vcom0")]

# volume recovery maps 
Mvrec = dcast(pred, long + lat ~ zname , value.var = "vrec")
Mvrec$NL = 0
coordinates(Mvrec) <- ~ long+lat
gridded(Mvrec) <- TRUE
proj4string(Mvrec) <- proj4string(grd)

cost_vrec <- stack(sapply(1:nrow(df_zones), function(i) {
  raster(MvextReal,as.character(df_zones$zname[i])) - raster(Mvrec,as.character(df_zones$zname[i]))
}))
