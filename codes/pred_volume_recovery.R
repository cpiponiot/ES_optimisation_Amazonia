
dt = data.table(expand.grid(id=grd$id, zone = df_zones$zone))
dt = merge(dt,data.table(id=grd$id, coordinates(grd)), by="id")
dt = merge(dt, df_zones, by = "zone")
pred = merge(dt, pars, by=c("long","lat"))

pred$V0 = volume(pred$t0, pred$aG, pred$aM, pred$bG, pred$bM, pred$theta, 0.2)
pred$Vcom0 = pred$V0*pred$om0
pred$vextReal = apply(pred[,c("vext","Vcom0")],1,min)
# add actual extracted volume to DTinput (we'll need it for carbon emissions estimation)
if (is.null(DTinput$vextReal)) 
  DTinput = merge(DTinput, pred[,c("id","zone","vextReal")], by=c("id","zone"))
# deltaV: total volume loss
pred$deltaV = apply(cbind(pred$vext/(pred$om0^(1-pred$rho)),pred$V0), 1, min)
dt0 = pred[,.(t1=t0Prediction(t0, deltaV, aG,aM,bG,bM,theta,0.2)), .(long,lat,vext,trot)]
pred = merge(pred, dt0, by=c("long","lat","vext","trot"))

# post-logging proprotion of commercial species: om1
pred$om1 = (pred$om0*pred$V0 - pred$vextReal) / (pred$V0 - pred$deltaV)
pred$om1[pred$deltaV == pred$V0] <- 0
  
# proportion of commercial species at the end of the cutting cycle
DTom2 = pred[,.(om2 = omega_t(tvec=c(t1+trot),t1,om1,om0,aG,aM,bG,bM,theta,intercept, slope)),.(vext,trot,long,lat)]
pred = merge(pred, DTom2, by=c("vext","trot","long","lat"))


# recovered timber volume
pred$vrec = volume(t=pred$t1+pred$trot, ag=pred$aG, am=pred$aM,
                        bg=pred$bG, bm=pred$bM,th=pred$theta, pdef=0.2)*pred$om2 - (pred$V0 - pred$deltaV)*pred$om1

Mvrec = dcast(pred, long + lat ~ zname , value.var = "vrec")
Mvrec$NL = 0
coordinates(Mvrec) <- ~ long+lat
gridded(Mvrec) <- TRUE
proj4string(Mvrec) <- proj4string(grd)

MvextReal = dcast(pred, long + lat ~ zname , value.var = "vextReal")
MvextReal$NL = 0
coordinates(MvextReal) <- ~ long+lat
gridded(MvextReal) <- TRUE
proj4string(MvextReal) <- proj4string(grd)


Mvcom0 = subset(pred, zone==1 )[,c("long","lat","Vcom0")]
