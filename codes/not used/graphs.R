library(data.table)
library(ggplot2)

## maps and zones ##
source("codes/getMaps.R")
grd$pu = 1:nrow(grd)
df_zones <- data.table(zone = 1:10, 
                       zname = c("LS","LM","LL","MS","MM","ML","HS","HM","HL","NL"),
                       vext = c(rep(10*1:3, each=3),0), 
                       trot = c(rep(c(15,30,65), 3),1))

## saved outputs ##
op = list.files("outputs")
i1 <- gsub("changeDemand|_|.Rdata", "", op[grep("changeDemand", op)])
i1 <- i1[i1!=""]
i2 <- gsub("changeCost|_|.Rdata", "", op[grep("changeCost", op)])
i2 <- i2[i2!=""]

changeDemandAll = c()

for (i in i1[!is.na(i1)]) {
  load(file = paste("outputs/changeDemand_",i,".Rdata",sep=""))
  changeDemand$i = i
  changeDemandAll = rbind(changeDemandAll, changeDemand)
}

# get initial volume per pixel
demandCosts = merge(changeDemandAll, df_zones, by="zone")
demandCosts = merge(demandCosts, 
                    data.table(coordinates(grd), area=grd$area, 
                               pHarv = grd$pHarv, pHarvAR = grd$pHarvAR, 
                               acs0=grd$acs, Rich0 = grd$mammals+grd$amphi), 
                    by=c("long","lat"))
demandCosts$areaTot = demandCosts$area*demandCosts$pHarvAR  ## total area (accessible & not)
demandCosts$areaLogging = demandCosts$area*demandCosts$pHarv  ## available area
demandCosts$areaLogging[grep("sharing", demandCosts$scenario)] = demandCosts$areaTot[grep("sharing", demandCosts$scenario)]
demandCosts$areaLogging[demandCosts$zone==10] = 0

demandAll = demandCosts[,.(areaTot = sum(areaLogging)*1e-4,
                             logIntens = weighted.mean(x = vextReal*trot, w = areaLogging),
                             cutCycle = weighted.mean(x = trot, w = areaLogging),
                             TtimbLoss=(1-sum(timber*areaLogging)/sum(Vcom0*areaTot))*100,
                             TcarbLoss=(1-sum(carbon*areaLogging)/sum(acs0*areaTot))*100,
                             TdivLoss=(1-sum(biodiversity*areaLogging)/sum(Rich0*areaTot))*100),
                          .(demand,scenario,i)]
demandAll = melt(demandAll, id.vars = c("demand","scenario","i") )

demandFinal = demandAll[,.(inf = quantile(value, 0.025), 
                           med = quantile(value, 0.5), 
                           sup = quantile(value, 0.975),
                           maxL = value[i=="maxL"]), .(demand, scenario,variable)]

levels(demandFinal$variable) <- c("(a) Total area logged (Mha)",
                                  "(b) Mean logging intensity (m3/ha)", 
                                  "(c) Mean cutting cycle length (yr)",
                                  "(d) Timber retained (%)",
                                  "(e) Carbon retained (%)", 
                                  "(f) Biodiversity retained (%)" ) 

demandFinal$scenario <- as.factor(demandFinal$scenario)
levels(demandFinal$scenario) <- c("Balanced","Carbon","Land sharing","STY + Land sharing","STY","Timber") 


ggplot(demandFinal, aes(x=demand, y=maxL, ymin = inf, ymax = sup, colour=scenario, fill=scenario)) + 
  geom_hline(data = data.frame(variable = levels(demandFinal$variable),
                               h = c(rep(c(NA,100), each=3))), aes(yintercept = h), lty=2) + 
  geom_ribbon(alpha=0.1, colour=NA) +
  geom_line() + 
  facet_wrap( ~ variable ,  nrow=3, scale="free_y", dir = "v") + 
  theme_bw() + 
  labs(x=expression("Timber production (M"*m^3*yr^{-1}*")"), y="",colour="Strategies", fill="Strategies") 
ggsave("LaTeX/graphs/increasingDemand_errors.pdf", height=6, width=8)
