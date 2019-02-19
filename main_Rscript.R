library(prioritizr)
library(sp)
library(rgdal)
library(raster)
library(data.table)
library(rstan)
library(ggpubr)
library(parallel)
library(ggtern)

solveProblems <- FALSE
solveProblemsParallel <- FALSE

#### study region & maps ####

load("data/maps.Rdata")
source("codes/proportion_areas.R")

pdf("graphs/harv_areas.pdf", height=5, width=7)
par(mar=c(1,1,1,1))
plot(borders, col = "#F5F5F5", axes=FALSE)
plot(mask_FC, add=TRUE, legend = FALSE, col="#006400")
plot(wdpa_am_forest,legend=FALSE, col="#FF8C00", add=TRUE)
plot(harvestableAreas,legend=FALSE, col="#7FFF00", add=TRUE)
legend(x="bottomright", fill=c("#FF8C00","#006400","#7FFF00"),
       legend = paste(c("Protected forest","Inaccessible forest","Available forest"), " (",c(Pwdpa,Punacc,Pharv),"%)",sep=""), 
       bg = "white",box.col = "white")
dev.off()


#### ES costs & prioritisation ####

### (1) Timber ###

load("data/grd.Rdata")
load("data/parametersTimber.Rdata")

df_zones <- data.table(zone = 1:10, 
                       zname = c("LS","LM","LL","MS","MM","ML","HS","HM","HL","NL"),
                       vext = c(rep(10*1:3, each=3),0), 
                       trot = c(rep(c(15,30,65), 3),1))

DTinput <- data.table(expand.grid(pu = grd$pu, zone = 1:9))
DTinput <- merge(DTinput, data.frame(grd), by="pu")
DTinput <- merge(DTinput, df_zones, by="zone")

source("codes/modelPredictions.R")
source("codes/pred_volume_recovery.R")

## per ha
feat_prod <- sapply(1:nrow(df_zones), function(i) {
  raster(MvextReal,as.character(df_zones$zname[i]))/df_zones$trot[i]
})

cost_vrec <- stack(sapply(1:nrow(df_zones), function(i) {
  raster(MvextReal,as.character(df_zones$zname[i])) - raster(Mvrec,as.character(df_zones$zname[i]))
}))


### (2) Carbon ###

## damage - graph ##

source("codes/new_damage_model.R")
source("codes/params_carbon_recov.R")
source("codes/functions_carbon.R")

## per ha
pdefor_maxL = (0.12+0.25)/2/17.95 + exp(-6.12) + (0.04 + 0.1)/2*15.10/(15.10+15.32)
Carbon <- DTinput[,.(Carbon = mean(acs*pdefor_maxL + dam + vext*WDext -
                                     f_recov(t=1:trot, loss = loss, acs0 = acs, 
                                             aSg = SurvC$aSg, aRr = RecrC$aRr, aRg = RecrC$aRg, 
                                             bSg = bSg, bSm = bSm, bRr = bRr, bRg = bRg, bRm = bRm, 
                                             eta = RecrC$eta)), 
                     Crecov =  mean( f_recov(t=1:trot, loss = loss, acs0 = acs, 
                                             aSg = SurvC$aSg, aRr = RecrC$aRr, aRg = RecrC$aRg, 
                                             bSg = bSg, bSm = bSm, bRr = bRr, bRg = bRg, bRm = bRm, 
                                             eta = RecrC$eta))),
                  .(trot,vext,pu)]
DTinput = merge(DTinput, Carbon, by = c("trot","vext","pu"))
DTinput = DTinput[,c("trot","vext","pu","zone","zname","pHarv","area","acs","long","lat", "Carbon","Crecov")]

## carbon cost rasters
DTinput$Carbon <- apply(cbind(0, DTinput$Carbon), 1, max)
cost_carbon <- dcast(DTinput, long + lat ~ zone, value.var = "Carbon")
colnames(cost_carbon) <- c("long","lat","LS","LM","LL","MS","MM","ML","HS","HM","HL")
cost_carbon$NL = 0
coordinates(cost_carbon) <- ~ long+lat
gridded(cost_carbon) <- TRUE
proj4string(cost_carbon) <- proj4string(grd)
cost_carbon <- stack(cost_carbon)
names(cost_carbon) <- df_zones$zname


### (3) Biodiversity ###

### new burivalova coefficients
richLoss = read.csv("data/Burivalova2014DataNeotropics.csv")
# fit log-log model (vext>0 and pRich > 0 but not always < 1)
# and intercept must be 1 
if (solveProblems) {
  y = 1-richLoss$pRichLoss; N = nrow(richLoss)
  x = richLoss$vext; Groups = as.numeric(richLoss$taxo); K = nlevels(richLoss$taxo)
  stan_biodiv = data.table(do.call(cbind, rstan::extract(stan("codes/linear_reg.stan", chains=1))))
  colnames(stan_biodiv) = c(paste(rep(c("slope","sd"), each=2), rep(c("Amp","Mam"),2), sep=""), "lp__")
  save(stan_biodiv, file="data/effect_biodiv_stan.Rdata")
} else (load("data/effect_biodiv_stan.Rdata"))

## per ha 
cost_diversity <- stack(sapply(df_zones$vext, function(x) {
  grd$rich_mamm = grd$mammals*(0.0144*x)
  grd$rich_mamm[grd$rich_mamm<0] <- 0
  grd$rich_amph = grd$amphi*(0.0153*x)
  grd$rich_amph[grd$rich_amph<0] <- 0
  return(crop(raster(grd, 'rich_mamm') + raster(grd, 'rich_amph'),extent(cost_carbon)))
}))

### (4) Intact forest landscape ###

feature_ifl = sapply(1:nrow(df_zones), function(i) {
  sapply(unique(grd$region), function(reg) {
    grd$ifl_zr = (grd$area*grd$pIFL)*(grd$region==reg)*(1-grd$pHarv*(df_zones$zname[i]!="NL"))
    crop(raster(grd,"ifl_zr"), extent(cost_carbon))
  })
})
## protect 80% of remaining IFL per region 
## total IFL area per ecoregion
areaIFL = tapply(grd$pIFL*grd$area,grd$region,sum)



##### Solving optimisation problems #####

##### (1) Combination of ES costs ####

source("codes/solve_problem.R")
source("codes/createTargetsCostsFeatures.R")

cost_comb = expand.grid(alphaC = seq(0,1,0.1), alphaB = seq(0,1,0.1))
cost_comb = subset(cost_comb, alphaC +alphaB <= 1)
cost_comb$alphaV = round(1 - cost_comb$alphaC - cost_comb$alphaB, digits=1)

if (solveProblemsParallel){
  
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores)
  
  clusterEvalQ(cl, library(prioritizr))
  clusterEvalQ(cl, library(data.table))
  
  clusterExport(cl, varlist = c("costMaps","featureMaps","areaIFL","targetsList","solve_problem","grd","cost_comb"))
  
  changeCost <- parSapply(cl, 1:nrow(cost_comb), function(i){
    
    costsWeighted <- cost_comb[i,1] * costMaps$carbon + cost_comb[i,2] * costMaps$biodiversity + cost_comb[i,3] * costMaps$timber 
    
    ##  optimization
    solZones <- solve_problem(costsWeighted, featureMaps$base, targetsList$base, i, timeLim = 600)
    solZones = data.table(solZones, alphaC = cost_comb$alphaC[i], 
                          alphaB = cost_comb$alphaB[i], alphaV = cost_comb$alphaV[i])
    return(solZones)
    
  })
  stopCluster(cl)
  
  changeCost = lapply(1:dim(changeCost)[2], function(i) { 
    data.table(do.call(cbind, changeCost[,i])) })
  
  changeCost = do.call(rbind, changeCost)
  
  save(changeCost, file="outputs/changeCost.Rdata")
  
  
} else {load("outputs/changeCost.Rdata")}

costs_analysis <- merge(changeCost, data.frame(grd)[,c("long","lat","area","pHarv")], by=c("long","lat"))

costs_analysis = costs_analysis[,.(carbon = raster::extract(costMaps$carbon[[zone]], 
                                                            cbind(long,lat)),
                                   biodiversity = raster::extract(costMaps$biodiversity[[zone]],
                                                                  cbind(long,lat)),
                                   timber = raster::extract(costMaps$timber[[zone]], 
                                                            cbind(long,lat)),
                                   long=long,lat=lat, alphaC,alphaB,alphaV, 
                                   area = area, pHarv=pHarv),.(zone)]

## for now: fix later timber maps ##
for (i in which(is.na(costs_analysis$timber))) {
  costs_analysis$timber[i] <- raster::extract(costMaps$timber[[costs_analysis$zone[i]]], costs_analysis[i, c("long","lat")]+cbind(1,1))
} 

costsTot = costs_analysis[,.(Ccost = sum(carbon*area*pHarv), 
                             Vcost = sum(timber*area*pHarv), 
                             Bcost = sum(biodiversity*area*pHarv)),.(alphaC,alphaB,alphaV)]
costsTot$Ccost_st = costsTot$Ccost/costsTot[alphaC==1]$Ccost*100
costsTot$Bcost_st = costsTot$Bcost/costsTot[alphaB==1]$Bcost*100
costsTot$Vcost_st = costsTot$Vcost/costsTot[alphaV==1]$Vcost*100
costsTot$Ccost_st[costsTot$Ccost_st<100] <- 100 ## fix that

d1 = abs(costsTot$Ccost_st - costsTot$Bcost_st)
d2 = abs(costsTot$Ccost_st - costsTot$Vcost_st)
d3 = abs(costsTot$Vcost_st - costsTot$Bcost_st)
which_balanced = which.min((d1+d2+d3)/3)


### (2) Scenario comparision ###

coeffs_balanced <- unlist(costsTot[which_balanced, c("alphaV","alphaC","alphaB")])
source("codes/createTargetsCostsFeatures.R")

if (solveProblemsParallel){
  
  library(parallel)
  
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores)
  
  clusterEvalQ(cl, library(prioritizr))
  clusterEvalQ(cl, library(data.table))
  
  clusterExport(cl, varlist = c("costMaps","featureMaps","areaIFL","targetsList","solve_problem","grd"))
  
  scenariOptim <- parSapply(cl, c(35, seq(10,80,10)*1e6), function(TD){
    
    for (y in 1:length(targetsList)) 
      targetsList[[y]]$target[1] <- TD
    
    ## increase time limit for sty-landsharing 
    TLim = c(rep(1000,5), 2000, 3000, 5000)[TD*1e-7]
    
    ## carbon optimization
    solZones <- rbind(solve_problem(costMaps$timber, featureMaps$base, targetsList$base, "volume", timeLim = 1000), 
                      solve_problem(costMaps$carbon, featureMaps$base, targetsList$base, "carbon", timeLim = 1000),
                      solve_problem(costMaps$biodiversity, featureMaps$base, targetsList$base, "biodiv", timeLim = 1000),
                      solve_problem(costMaps$balanced, featureMaps$base, targetsList$base, "balanced", timeLim = 1000),
                      solve_problem(costMaps$balanced, featureMaps$medium, targetsList$base, "medium", timeLim = 1000), 
                      solve_problem(costMaps$balanced, featureMaps$STY, targetsList$STY, "STY", timeLim = 1000),
                      solve_problem(costMaps$sharing, featureMaps$sharing, targetsList$base, "sharing", timeLim = 1000), 
                      solve_problem(costMaps$sharing, featureMaps$sharingSTY, targetsList$STY, "sharingSTY", timeLim = TLim))
    solZones$demand = TD/1e6
    return(solZones)
    
  })
  stopCluster(cl)
  
  scenariOptim = lapply(1:dim(scenariOptim)[2], function(i) { 
    data.table(do.call(cbind, scenariOptim[,i])) })
  
  scenariOptim = do.call(rbind, scenariOptim)
  
  scenariOptim$zone = as.numeric(scenariOptim$zone)
  scenariOptim$demand = as.numeric(scenariOptim$demand)
  scenariOptim$long = as.numeric(scenariOptim$long)
  scenariOptim$lat = as.numeric(scenariOptim$lat)
  
  save(scenariOptim, file="outputs/scenariOptim.Rdata")
  
} else {load("outputs/scenariOptim.Rdata")}

### pixel area
scenariOptim = merge(scenariOptim, data.table(coordinates(grd),area = grd$area*grd$pHarv, areaAR=grd$area*grd$pHarvAR), by=c("long","lat"))
scenariOptim$area[grep("sharing",scenariOptim$scenario)] <- scenariOptim$areaAR[grep("sharing",scenariOptim$scenario)]
scenariOptim = merge(scenariOptim, df_zones, by="zone")
scenariOptim = subset(scenariOptim, area>0)

### scenarios names
scenariOptim$scenario <- as.factor(scenariOptim$scenario)
levels(scenariOptim$scenario) <- c("Balanced","Biodiversity","Carbon","Current","Road building","STY + Road building","STY","Timber") 
scenariOptim$scenario <- factor(scenariOptim$scenario, levels = c("Timber","Carbon","Biodiversity","Balanced","Current","STY","Road building","STY + Road building") )

### ES costs 
scenariOptim = scenariOptim[,.(Cemi = raster::extract(cost_carbon[[zone]], cbind(long,lat)), 
                               timbLoss = raster::extract(cost_vrec[[zone]], cbind(long,lat)), 
                               biodLoss = raster::extract(cost_diversity[[zone]], cbind(long,lat)), 
                               vextReal = raster::extract(feat_prod[[zone]], cbind(long,lat)), 
                               long=long, lat=lat, demand = demand, scenario=scenario),.(zone)]
# get initial volume per pixel
scenariOptim = merge(scenariOptim, Mvcom0, by=c("long","lat"))

scenariOptim = merge(scenariOptim, df_zones, by="zone")
scenariOptim = merge(scenariOptim, 
                     data.table(coordinates(grd), area=grd$area, 
                                pHarv = grd$pHarv, pHarvAR = grd$pHarvAR, 
                                acs0=grd$acs, Rich0 = grd$mammals+grd$amphi), 
                     by=c("long","lat"))
scenariOptim$areaTot = scenariOptim$area*scenariOptim$pHarvAR  ## total area (accessible & not)
scenariOptim$areaLogging = scenariOptim$area*scenariOptim$pHarv  ## available area
scenariOptim$areaLogging[grep("build", scenariOptim$scenario)] = scenariOptim$areaTot[grep("build", scenariOptim$scenario)]
scenariOptim$areaLogging[scenariOptim$zname=="NL"] <- 0

scenCost = scenariOptim[demand == 35,.(timber=(1-sum(timbLoss*areaLogging)/sum(Vcom0*areaTot))*100,
                                       carbon=(1-sum(Cemi*areaLogging)/sum(acs0*areaTot))*100,
                                       biodiv=(1-sum(biodLoss*areaLogging)/sum(Rich0*areaTot))*100),
                        .(scenario)]
scenCost = melt(scenCost, id.vars = c("scenario"), variable.name = "ES", value.name = "loss")

demandFinal = scenariOptim[,.(areaTot = sum(areaLogging)*1e-4,
                              logIntens = weighted.mean(x = vextReal*trot, w = areaLogging),
                              cutCycle = weighted.mean(x = trot, w = areaLogging),
                              TtimbLoss=(1-sum(timbLoss*areaLogging)/sum(Vcom0*areaTot))*100,
                              TcarbLoss=(1-sum(Cemi*areaLogging)/sum(acs0*areaTot))*100,
                              TdivLoss=(1-sum(biodLoss*areaLogging)/sum(Rich0*areaTot))*100),
                           .(demand,scenario)]
demandFinal = melt(demandFinal, id.vars = c("demand","scenario") )

#### (3) Graphs ####

### maps ### 

palette <- function (name, indices = c(7,5,3)) {
  RColorBrewer::brewer.pal(9, name)[indices]
}
colours <- c(as.vector(sapply(c("Blues", "Purples","Reds"), palette)), "forestgreen")
colour_palette <- c("LS"= colours[1], "LM"=colours[2], "LL"=colours[3],
                    "MS"=colours[4], "MM"=colours[5], "ML"=colours[6],
                    "HS"=colours[7], "HM"=colours[8], "HL"=colours[9], "NL"=colours[10])

g1 <- ggplot(subset(scenariOptim, demand == 35)) +
  geom_point(aes(x=long,y=lat,colour=zname, size=area/1e3)) +
  theme_bw() + coord_fixed() + facet_wrap( ~ scenario, nrow = 2) +
  scale_colour_manual(name = "Zone", values = colour_palette) +
  labs(size = expression("Area available for logging ("*1000*" k"*m^2*")"), 
       x="",y="") +
  theme(legend.position = "top", text = element_text(size=40),
        axis.ticks=element_blank(), axis.text.x=element_blank(), 
        axis.text.y=element_blank(), 
        strip.text = element_text(colour = 'black'), 
        strip.background = element_blank()) +
  guides(colour=FALSE)

zone_legend <- ggplot(df_zones[-10],aes(x=vext,y=as.factor(trot),fill=zname,label=zname)) + geom_raster() + scale_fill_manual(values = colour_palette)+geom_text(size=15) + labs(x ="Extracted volume", y = "Cutting cycle (yrs)") + theme(legend.position="none", aspect.ratio=1, panel.background = element_blank(), text = element_text(size=30))

ggarrange(g1, zone_legend, ncol=2, widths=c(6,1))
ggsave("graphs/mapsScenarios.pdf", height=12, width=30)


### associated ES costs 

col_scenarios <- c("#6495ED","#458B00", "#E5C616", "#CD6600", "#E9967A","#483D8B", "#D33B44", "#8B008B")

levels(scenFinal$ES) = c("(a) Timber","(b) Carbon","(c) Biodiversity")
ggplot(scenCost, aes(x=scenario, fill=scenario, y=loss-100)) + 
  geom_histogram(stat="identity") + 
  facet_grid(~ES) + scale_fill_manual(values= col_scenarios) +
  labs(y="Variation (% initial value)", fill="Strategy", x="Strategy") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.background = element_rect(fill="white",colour="white"), 
        panel.background = element_rect(fill="white", colour = "black"),
        panel.grid = element_blank())
ggsave("graphs/costsScenario.pdf", height=4,width=7)


### ES = f(timber demand, scenario) ###

levels(demandFinal$variable) <- c("(a) Total area logged (Mha)",
                                  "(b) Mean logging intensity (m3/ha)", 
                                  "(c) Mean cutting cycle length (yr)",
                                  "(d) Timber retained (%)",
                                  "(e) Carbon retained (%)", 
                                  "(f) Biodiversity retained (%)" ) 

ggplot(demandFinal, aes(x=demand, y= value, colour=scenario))+ 
  geom_hline(data = data.frame(variable = levels(demandFinal$variable), h = c(rep(c(NA,100), each=3))),
             aes(yintercept = h), lty=2) + 
  geom_line() + scale_colour_brewer(palette = "Set1") +
  facet_wrap( ~ variable ,  nrow=3, scale="free_y", dir = "v") + 
  theme(panel.background = element_rect(fill="white", colour = "black"),
        panel.grid = element_blank(), strip.background = element_blank(), 
        legend.key = element_blank()) + 
  scale_color_manual(values= col_scenarios)+ 
  labs(x=expression("Timber production (M"*m^3*yr^{-1}*")"), y="",colour="Strategies") 
ggsave("graphs/increasingDemand.pdf", height=6, width=8)



##### Supplementary graphs #####

### damage - explanation diagram ###
pdf("LaTeX/graphs/schemaDam.pdf", height=2, width=3)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(c(0,1), c(0,1.1), col="white", bty="n",yaxt="n",xaxt="n",xlab="",ylab="")
rect(xleft = 0,ybottom = 0,xright = 1,ytop = 1,lwd=2)
rect(xleft = 0, ybottom = 0.8, xright = 1, ytop = 1, density = 5)
rect(xleft = 0,ybottom = 0.2,xright = 0.6,ytop = 0.8, col="grey")
segments(x0=0.1, y0=1,y1=1.08)
text(x=0.5,y=0.9, "Cext")
text(x=0.3,y=0.5, "Cdam")
text(x=0.1,y=1.1, "C0")
text(x=0.8,y=.4, "Cmin")
dev.off()

## damage - graph ##
ggplot(CIpred, aes(x=RatioExt)) + 
  geom_point(aes(y=RatioDam, colour=Forest, size=PlotSurface)) + 
  geom_line(aes(y=med), lwd=1.2) + 
  geom_ribbon(aes(ymin=inf_err, ymax=sup_err), alpha=0.1) +
  theme_bw() + theme(legend.position = "none") +
  labs(x="Proportion of carbon extracted", y="Proportion of carbon damaged")
ggsave("graphs/damModel.pdf", height=5, width=6)


## WDext 
pdf("graphs/map_WDext.pdf", height=4, width=5)
plot(raster(grd, "WDext"), col=rev(heat.colors(20))[-c(1,2)])
dev.off()

### ternary plots - cost analysis

limsPlots = range(costsTot[,grep("_st",colnames(costsTot)),with=F]) +c(-10,0)
paletteTern = colorRampPalette(c("blue","yellow", "orange", "red"))(100)

ggtern(data=costsTot,aes(x=alphaC,y=alphaB,z=alphaV, value=Ccost_st)) + 
  stat_interpolate_tern(geom="polygon", method=lm, colour="black", formula = value~x+y,
                        n=100,aes(fill=..level..),expand=1)  +
  # geom_point(cex=10) +
  theme(legend.position = "none") + 
  scale_fill_gradientn(colours = paletteTern, limits=limsPlots) +
  labs(x="Carbon\nweight", y="Biodiversity\nweight", z="Timber\nweight", fill="Carbon cost (%)") +
  geom_point(data = costsTot[which_balanced], size = 5)
ggsave("graphs/carbonLoss.pdf", height=4, width=4.5)

ggtern(data=costsTot,aes(x=alphaC,y=alphaB,z=alphaV, value=Bcost_st)) + 
  stat_interpolate_tern(geom="polygon", method=lm,n=100, aes(fill=..level..),expand=1, colour="black")  +
  # geom_point(cex=10) + 
  theme(legend.position = "none") + 
  scale_fill_gradientn(colours = paletteTern, limits=limsPlots) + 
  labs(x="Carbon\nweight", y="Biodiversity\nweight", z="Timber\nweight", fill="Diversity cost (%)") +
  geom_point(data = costsTot[which_balanced], size = 5)
ggsave("graphs/biodivLoss.pdf",  height=4, width=4.5)

g <- ggtern(data=costsTot,aes(x=alphaC,y=alphaB,z=alphaV, value=Vcost_st)) + 
  stat_interpolate_tern(geom="polygon", method=lm,n=100,
                        aes(fill=..level..), colour="black",expand=1)  +
  # geom_point(cex=10) + 
  scale_fill_gradientn(colours = paletteTern, limits=limsPlots) + 
  labs(x="Carbon\nweight", y="Biodiversity\nweight", z="Timber\nweight", fill="ES \nloss (%)") +
  geom_point(data = costsTot[which_balanced], size = 5)
g + theme(legend.position = "none")
ggsave("LaTeX/graphs/timberLoss.pdf",  height=4, width=4.5)

as_ggplot(get_legend(g + guides(colour = guide_colourbar(barheight = 15))))
ggsave("LaTeX/graphs/legendESloss.pdf", height=4, width=1.5)


## maps with changing demand 
scenariOptim$demand2 = paste(scenariOptim$demand, "Mm3/yr")
ggplot(subset(scenariOptim, demand != 35)) +
  geom_point(aes(x=long,y=lat,colour=zname, size=area/1e3)) +
  theme_bw() + coord_fixed() + facet_grid(demand2 ~ scenario) +
  scale_colour_manual(name = "Zone",values = colour_palette) +
  labs(size = "Area available for logging (1000 km2)", x="", y="") +
  theme(legend.position = "top", text = element_text(size=50),
        axis.ticks=element_blank(), axis.text.x=element_blank(), 
        axis.text.y=element_blank()) +
  guides(colour=FALSE)
ggsave("graphs/mapsChangeDemand.pdf", height=40, width=40)

