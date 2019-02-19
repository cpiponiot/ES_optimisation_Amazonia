library(data.table)
library(ggplot2)
library(rstan)

## open data
damdata = read.csv("data/02_damdata_mort_TmFO.csv")

## loss = proportion of biomass loss = dAGB / AGB0
# damdata$loss = (damdata$AGBdam+damdata$AGBext)/damdata$AGB0

damdata$RatioExt = damdata$AGBext/damdata$AGB0

## logit and logistic functions
logit = function(p) log(p/(1-p))
logist = function(x) 1/(1+exp(-x))

## new model: logit(loss) = a * log(vext) + b
ggplot(damdata, aes(x=logit(RatioExt), y=logit(RatioDam), colour=Forest)) + 
  geom_point() + geom_smooth(method="lm", colour=1,alpha=0.4)

x=logit(damdata$RatioExt)
y=logit(damdata$RatioDam)
N = length(x)
plotsize = damdata$PlotSurface  # give weight proportional to plot size?

summary(lm(y~x))

# stan_linreg = stan("linear_reg_noIntercept.stan", iter = 500, warmup = 250)
# traceplot(stan_linreg)
# pars = data.table(do.call(cbind,extract(stan_linreg)))
# pars$iter = 1:nrow(pars)
pars = read.csv("data/stanNewDamageModel.csv")

## prediction dataframe with all parameter values for all plots
DTpred = data.table(expand.grid(iter = pars$iter, Plot = damdata$Plot))
DTpred = merge(DTpred, pars, by = "iter")
DTpred = merge(DTpred, damdata, by = "Plot")

## prediction of the proportion of biomass loss (without the error term)
DTpred$Pdam = logist(DTpred$slope*logit(DTpred$RatioExt))#  + DTpred$intercept)
## prediction with the error term 
DTpred$Pdam_err = logist(DTpred$slope*logit(DTpred$RatioExt) + rnorm(nrow(DTpred),0,DTpred$sigma)) #  + DTpred$intercept)

## confidence intervals 
CIpred = DTpred[,.(inf = quantile(Pdam, 0.025),
                   med = quantile(Pdam, 0.5),
                   sup = quantile(Pdam, 0.975),
                   inf_err = quantile(Pdam_err, 0.025),
                   sup_err = quantile(Pdam_err, 0.975)),
                .(RatioExt, RatioDam, Plot, Forest,PlotSurface)]
