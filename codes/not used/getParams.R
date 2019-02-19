library(rstan)
library(data.table)

#### volume  ####

load("C:/Users/camille.piponiot/Google Drive/volume_recovery/stan models/stan_covar.Rdata")
pars_vol <- rstan::extract(stan_covar)
load("C:/Users/camille.piponiot/Google Drive/volume_recovery/stan models/stan_rho.Rdata")
pars_rho <- rstan::extract(stan_rho)
load("C:/Users/camille.piponiot/Google Drive/volume_recovery/stan models/stan_pR.Rdata")
pars_pR <- rstan::extract(stan_pR)

pars <- data.table(i = sample(200, 100))

pars$lambda_ti <- pars_vol$lambda_ti[pars$i]
pars$dti <- mean(pars_vol$dti[pars$i,])
pars$theta <- pars_vol$theta[pars$i]
pars$sigmaV <- pars_vol$sigma_V[pars$i]
pars$bG <- pars_vol$bP[pars$i]
pars$bM <- pars_vol$bM[pars$i]
pars$rho <- pars_rho$rho[pars$i]
pars$intercept <- pars_pR$intercept[pars$i]
pars$slope <- pars_pR$slope[pars$i]
pars$pdef <- rbeta(nrow(pars), 2, 8)


#### carbon loss ####
pars$Pdefor = runif(nrow(pars),0.12,0.25)*rexp(nrow(pars), 17.95) + 
  rlnorm(nrow(pars), -6.12, 0.58) + 
  runif(nrow(pars), 0.04,0.1)*rbeta(nrow(pars), 15.10,15.32)
dam_pars <- read.csv("stanNewDamageModel.csv")
pars$theta_Cdam = dam_pars$slope[pars$i]
pars$sd_Cdam = dam_pars$sigma[pars$i]
  
#### carbon recovery ####
load("C:/Users/camille.piponiot/Google Drive/recovery/data/surv.Rdata")
SurvC <- data.table((as.matrix(surv, pars=c("alpha","beta1","beta2", names(surv)[grep("lambda",names(surv))]))[pars$i,]))
colnames(SurvC) <- c("aSg","bSg0","bSm0", 
                     paste(rep("l",12), c(rep(c("Sg","Sm"), each=5),c("Sg","Sm")), 
                           c(rep(c("acs","dacs","prec","seas","bkd"),2),rep("loss",2)), sep="_"))
pars = cbind(pars, SurvC)

## Recruits parameters ## 
load("C:/Users/camille.piponiot/Google Drive/recovery/data/newc.Rdata")
RecrC <- data.table((as.matrix(newc, pars=c("alpha1","alpha2","beta1","beta2","beta3","eta", names(newc)[grep("lambda",names(newc))]))[pars$i,]))
colnames(RecrC) <- c("aRr","aRg","bRr0","bRg0","bRm0","eta", 
                     paste(rep("l",18), c(rep(c("Rr","Rg","Rm"), each=5),c("Rr","Rg","Rm")), 
                           c(rep(c("acs","dacs","prec","seas","bkd"),3),rep("loss",3)), sep="_"))
pars = cbind(pars, RecrC)

#### Biodiviersity ####

load("effect_biodiv_stan.Rdata")
pars = cbind(pars, stan_biodiv[pars$i,grep("slope",colnames(stan_biodiv)),with=F])

### maximum likelihood parameters ###

parsmaxL <- data.table(i="maxL")
# volume
parsmaxL$lambda_ti <- pars_vol$lambda_ti[which.max(pars_vol$lp__)]
parsmaxL$dti <- mean(pars_vol$dti[which.max(pars_vol$lp__),])
parsmaxL$theta <- pars_vol$theta[which.max(pars_vol$lp__)]
parsmaxL$sigmaV <- 0
parsmaxL$bG <- pars_vol$bP[which.max(pars_vol$lp__)]
parsmaxL$bM <- pars_vol$bM[which.max(pars_vol$lp__)]
parsmaxL$rho <- pars_rho$rho[which.max(pars_rho$lp__)]
parsmaxL$intercept <- pars_pR$intercept[which.max(pars_pR$lp__)]
parsmaxL$slope <- pars_pR$slope[which.max(pars_pR$lp__)]
parsmaxL$pdef <- 0.2
# carbon emissions
parsmaxL$Pdefor = (0.12+0.25)/2/17.95 + 
  exp(-6.12) + (0.04 + 0.1)/2*15.10/(15.10+15.32)
parsmaxL$theta_Cdam = dam_pars$slope[which.max(dam_pars$lp__)]
parsmaxL$sd_Cdam = 0
# carbon recovery
maxLS = which.max(as.matrix(surv, pars="lp__"))
SurvC <- data.table(t(as.matrix(surv, pars=c("alpha","beta1","beta2", names(surv)[grep("lambda",names(surv))]))[maxLS,]))
colnames(SurvC) <- c("aSg","bSg0","bSm0", 
                     paste(rep("l",12), c(rep(c("Sg","Sm"), each=5),c("Sg","Sm")), 
                           c(rep(c("acs","dacs","prec","seas","bkd"),2),rep("loss",2)), sep="_"))
parsmaxL = cbind(parsmaxL, SurvC)

## Recruits parameters ## 
maxLR = which.max(as.matrix(newc, pars="lp__"))
RecrC <- data.table(t(as.matrix(newc, pars=c("alpha1","alpha2","beta1","beta2","beta3","eta", names(newc)[grep("lambda",names(newc))]))[maxLR,]))
colnames(RecrC) <- c("aRr","aRg","bRr0","bRg0","bRm0","eta", 
                     paste(rep("l",18), c(rep(c("Rr","Rg","Rm"), each=5),c("Rr","Rg","Rm")), 
                           c(rep(c("acs","dacs","prec","seas","bkd"),3),rep("loss",3)), sep="_"))
parsmaxL = cbind(parsmaxL, RecrC)
parsmaxL = cbind(parsmaxL, stan_biodiv[which.max(stan_biodiv$lp__),grep("slope",colnames(stan_biodiv)),with=F])

pars = rbind(pars, parsmaxL)



save(pars, file="pars.Rdata")
