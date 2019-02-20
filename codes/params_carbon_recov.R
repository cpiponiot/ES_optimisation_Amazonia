dam_pars <- read.csv("data/stanNewDamageModel.csv")
dam_pars <- dam_pars[which.max(dam_pars$lp),]
ratioExt <- (DTinput$vextReal*DTinput$WDext)/(DTinput$acs)

ratioDam <- logist(dam_pars$slope*logit(ratioExt))
DTinput$dam <-  ratioDam*(DTinput$acs-DTinput$vextReal*DTinput$WDext)

## % agb loss (minus roads: no recovery)
DTinput$loss <- sapply( ( DTinput$dam + DTinput$WDext * DTinput$vextReal ) / DTinput$acs, function(x) min(x, 1) ) * 100

## scale covariables
standardize <- read.csv2("data/mean_sd_recov.csv")
covar <- do.call(cbind, sapply(1:5, function(i) (DTinput[,c("loss", "acs", "prec", "seas", "bkd")[i],with=F] - standardize[1,c("loss","acs0", "prec","seas","bd")[i]]) / standardize[2,c("loss","acs0", "prec","seas","bd")[i]] ))

## Survivors parameters ##
load("data/stan_Crecov_surv.Rdata")
SurvC <- SurvC[which.max(SurvC$lp__), -"lp__"]

DTinput$bSg <- SurvC$bSg + covar %*% t(SurvC[,sapply(colnames(covar), function(i) grep(paste("Sg",i,sep="_"), colnames(SurvC))),with=F]) 
DTinput$bSg[DTinput$bSg < 0] <- min(DTinput$bSg[DTinput$bSg>0])

DTinput$bSm <- SurvC$bSm + covar %*% t(SurvC[,sapply(colnames(covar), function(i) grep(paste("Sm",i,sep="_"), colnames(SurvC))),with=F]) 
DTinput$bSm[DTinput$bSm < 0] <- min(DTinput$bSm[DTinput$bSm>0])

## Recruits parameters ## 
load("data/stan_Crecov_recr.Rdata")
RecrC <- RecrC[which.max(RecrC$lp__), -"lp__"]

## [!] negative beta values... -> recalibrate the model with a log-normal distribution? 
DTinput$bRr <- RecrC$bRr + covar %*% t(RecrC[,sapply(colnames(covar), function(i) grep(paste("Rr",i,sep="_"), colnames(RecrC))),with=F]) 
DTinput$bRr[DTinput$bRr < 0] <- min(DTinput$bRr[DTinput$bRr>0])

DTinput$bRg <- RecrC$bRg + covar %*% t(RecrC[,sapply(colnames(covar), function(i) grep(paste("Rg",i,sep="_"), colnames(RecrC))),with=F]) 
DTinput$bRg[DTinput$bRg < 0] <- min(DTinput$bRg[DTinput$bRg>0])

DTinput$bRm <- RecrC$bRm + covar %*% t(RecrC[,sapply(colnames(covar), function(i) grep(paste("Rm",i,sep="_"), colnames(RecrC))),with=F]) 
DTinput$bRm[DTinput$bRm < 0] <- min(DTinput$bRm[DTinput$bRm>0])
