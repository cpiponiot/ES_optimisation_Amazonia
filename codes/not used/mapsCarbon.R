source("codes/functions_carbon.R")

### damage
logist = function(x) 1/(1+exp(-x))
logit = function(p) log(p/(1-p))  

ratioExt <- (pred$vextReal*pred$WDext)/(pred$acs)

ratioDam <- logist(rnorm(nrow(pred), pred$theta_Cdam*logit(ratioExt), pred$sd_Cdam))
pred$dam <-  ratioDam*(pred$acs-pred$vextReal*pred$WDext)

pred$loss <- sapply((pred$dam+pred$WDext*pred$vextReal)/pred$acs, function(x) min(x,1))*100

## scale covariables
standardize <- read.csv2("C:/Users/camille.piponiot/Google Drive/carbon balance/bilanC_version200617/mean_sd_recov.csv")
covar <- do.call(cbind, sapply(1:5, function(i) (pred[,c("loss", "acs", "prec", "seas", "bkd")[i],with=F] - standardize[1,c("loss","acs0", "prec","seas","bd")[i]]) / standardize[2,c("loss","acs0", "prec","seas","bd")[i]] ))

## Survivors parameters ##
pred$bSg <- pred$bSg0[1] + covar %*% t(pred[1,sapply(colnames(covar), function(j) grep(paste("Sg",j,sep="_"), colnames(pred))),with=F]) 
pred$bSg[pred$bSg < 0] <- min(pred$bSg[pred$bSg>0])

pred$bSm <- pred$bSm0[1] + covar %*% t(pred[1,sapply(colnames(covar), function(j) grep(paste("Sm",j,sep="_"), colnames(pred))),with=F]) 
pred$bSm[pred$bSm < 0] <- min(pred$bSm[pred$bSm>0])

## Recruits parameters ## 
## [!] negative beta values... -> recalibrate the model with a log-normal distribution? 
pred$bRr <- pred$bRr0[1] + covar %*% t(pred[1,sapply(colnames(covar), function(i) grep(paste("Rr",i,sep="_"), colnames(pred))),with=F]) 
pred$bRr[pred$bRr < 0] <- min(pred$bRr[pred$bRr>0])

pred$bRg <- pred$bRg0[1] + covar %*% t(pred[1,sapply(colnames(covar), function(i) grep(paste("Rg",i,sep="_"), colnames(pred))),with=F]) 
pred$bRg[pred$bRg < 0] <- min(pred$bRg[pred$bRg>0])

pred$bRm <- pred$bRm0[1] + covar %*% t(pred[1,sapply(colnames(covar), function(i) grep(paste("Rm",i,sep="_"), colnames(pred))),with=F]) 
pred$bRm[pred$bRm < 0] <- min(pred$bRm[pred$bRm>0])


dfC <- pred[,.(Carbon = mean(acs*Pdefor + dam + vext*WDext -
                                     f_recov(t=1:trot, loss = loss, acs0 = acs, 
                                             aSg = aSg, aRr = aRr, aRg = aRg, 
                                             bSg = bSg, bSm = bSm, bRr = bRr, bRg = bRg, bRm = bRm, 
                                             eta = eta)), 
                     Crecov =  mean( f_recov(t=1:trot, loss = loss, acs0 = acs, 
                                             aSg = aSg, aRr = aRr, aRg = aRg, 
                                             bSg = bSg, bSm = bSm, bRr = bRr, bRg = bRg, bRm = bRm, 
                                             eta = eta))),
                  .(trot,vext,long,lat,zone)]

## carbon cost rasters
dfC$Carbon <- apply(cbind(0, dfC$Carbon), 1, max)
cost_carbon <- dcast(dfC, long + lat ~ zone, value.var = "Carbon")
colnames(cost_carbon) <- c("long","lat","LS","LM","LL","MS","MM","ML","HS","HM","HL","NL")
cost_carbon$NL = 0
coordinates(cost_carbon) <- ~ long+lat
gridded(cost_carbon) <- TRUE
proj4string(cost_carbon) <- proj4string(grd)
cost_carbon <- stack(cost_carbon)
names(cost_carbon) <- df_zones$zname
