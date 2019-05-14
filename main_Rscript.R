library(prioritizr)
library(sp)
library(rgdal)
library(raster)
library(data.table)
library(rstan)
library(ggpubr)
library(parallel)
library(ggtern)

solveProblems <- TRUE
current_demand <- 35 ## in Mm3

#### study region & maps ####

load("data/maps.Rdata")
load("data/grd.Rdata")
source("codes/proportion_areas.R")

pdf("graphs/harv_areas.pdf", height=5, width=7)
par(mar=c(1,1,1,1))
plot(borders, col = "#F5F5F5", axes=FALSE)
plot(forest_cover_90, add=TRUE, legend = FALSE, col="#006400")
plot(protected_forest,legend=FALSE, col="#FF8C00", add=TRUE)
plot(area_avail,legend=FALSE, col="#7FFF00", add=TRUE)
plot(borders, col = NA, add=TRUE)
legend(x="bottomright", fill=c("#FF8C00","#006400","#7FFF00"),
       legend = paste(c("Protected forest","Inaccessible forest","Available forest"), " (",c(Pwdpa,Punacc,Pharv),"%)",sep=""), 
       bg = "white",box.col = "white")
dev.off()


#### ES costs & prioritisation ####

### (1) Timber ###

load("data/parametersTimber.Rdata")

df_zones <- data.table(zone = 1:10, 
                       zname = c("LS","LM","LL","MS","MM","ML","HS","HM","HL","NL"),
                       vext = c(rep(10*1:3, each=3),0), 
                       trot = c(rep(c(15,30,65), 3),1))

DTinput <- data.table(expand.grid(id = grd$id, zone = 1:9))
DTinput <- merge(DTinput, data.frame(grd), by="id")
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
                  .(trot,vext,id)]
DTinput = merge(DTinput, Carbon, by = c("trot","vext","id"))
DTinput = DTinput[,c("trot","vext","id","zone","zname","pAreaAvail","area","acs","long","lat", "Carbon","Crecov")]

## carbon cost rasters
DTinput$Carbon <- apply(cbind(0, DTinput$Carbon), 1, max)
cost_carbon <- dcast(DTinput, long + lat ~ zone, value.var = "Carbon")
colnames(cost_carbon) <- c("long","lat","LS","LM","LL","MS","MM","ML","HS","HM","HL")
cost_carbon$NL = 0
coordinates(cost_carbon) <- ~ long +lat
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


##### Solving optimisation problems #####

##### (1) Combination of ES costs ####


if (solveProblems){
  
  source("codes/solve_problem.R")
  
  coeffs_balanced = c(1,1,1)/3
  source("codes/createTargetsCostsFeatures.R")
  
  cost_comb = expand.grid(alphaC = seq(0,1,0.1), alphaB = seq(0,1,0.1))
  cost_comb = subset(cost_comb, alphaC + alphaB <= 1)
  cost_comb$alphaV = round(1 - cost_comb$alphaC - cost_comb$alphaB, digits=1)
  
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

changeCost <- merge(changeCost, data.frame(grd)[,c("long","lat","area","pAreaAvail")], by=c("long","lat"))
changeCost[, areaLogging := area * pAreaAvail]

source("codes/rel_ES_costs.R")

costs_analysis <- changeCost[,.(timber = rel_ES_costs(zone, long, lat, areaLogging, "timber"), 
                                carbon = rel_ES_costs(zone, long, lat, areaLogging, "carbon"), 
                                biodiversity = rel_ES_costs(zone, long, lat, areaLogging, "biodiversity")),
                             .(alphaV, alphaC, alphaB)]

costsTot = melt(costs_analysis, id.vars = c("alphaV", "alphaC", "alphaB"), variable.name = "ES", value.name = "loss")


### (2) Scenario comparision ###

if (solveProblems){
  
  library(parallel)
  
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores)
  
  clusterEvalQ(cl, library(prioritizr))
  clusterEvalQ(cl, library(data.table))
  
  clusterExport(cl, varlist = c("costMaps","featureMaps","areaIFL","targetsList","solve_problem","grd"))
  
  scenariOptim <- parSapply(cl, c(current_demand, seq(10,80,10))*1e6, function(TD){
    
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
scenariOptim = merge(scenariOptim, data.table(coordinates(grd), 
                                              areaLogging = grd$area*grd$pAreaAvail, 
                                              areaUnprotec = grd$area*grd$pAreaUnprotect), by=c("long","lat"))
scenariOptim$areaLogging[grep("sharing", scenariOptim$scenario)] <- scenariOptim$areaUnprotec[grep("sharing",scenariOptim$scenario)]
scenariOptim = merge(scenariOptim, df_zones, by="zone")
scenariOptim = subset(scenariOptim, areaLogging > 0)

### scenarios names
scenariOptim$scenario <- as.factor(scenariOptim$scenario)
levels(scenariOptim$scenario) <- c("Balanced","Biodiversity","Carbon","Current","Road building","STY + Road building","STY","Timber") 
scenariOptim$scenario <- factor(scenariOptim$scenario, levels = c("Timber","Carbon","Biodiversity","Balanced","Current","STY","Road building","STY + Road building") )

## get vextreal 
scenariOptim = merge(scenariOptim, scenariOptim[,.(vextReal = raster::extract(feat_prod[[zone]], cbind(long,lat)),
                                                   long=long, lat=lat, demand = demand, scenario=scenario), .(zone)], 
                     by = c("long","lat","demand","scenario","zone"))

## get costs
scenariOptim$areaLogged = scenariOptim$areaLogging
scenariOptim$areaLogged[scenariOptim$zname=="NL"] <- 0
demandFinal <- scenariOptim[,.(areaTot = sum(areaLogged)*1e-6,
                               logIntens = weighted.mean(x = vextReal*trot, w = areaLogged),
                               cutCycle = weighted.mean(x = trot, w = areaLogged),
                               timber = rel_ES_costs(zone, long, lat, areaLogged, "timber"), 
                               carbon = rel_ES_costs(zone, long, lat, areaLogged, "carbon"), 
                               biodiv = rel_ES_costs(zone, long, lat, areaLogged, "biodiversity")),
                            .(demand, scenario)]
demandFinal = melt(demandFinal, id.vars = c("demand","scenario"))


#### (3) Graphs ####

### maps ### 

palette <- function (name, indices = c(7,5,3)) {
  RColorBrewer::brewer.pal(9, name)[indices]
}
colours <- c(as.vector(sapply(c("Blues", "Purples","Reds"), palette)), "forestgreen")
colour_palette <- c("LS"= colours[1], "LM"=colours[2], "LL"=colours[3],
                    "MS"=colours[4], "MM"=colours[5], "ML"=colours[6],
                    "HS"=colours[7], "HM"=colours[8], "HL"=colours[9], "NL"=colours[10])

g1 <- ggplot(subset(scenariOptim, demand == current_demand)) +
  geom_point(aes(x=long,y=lat,colour=zname, size=areaLogging/1e6)) +
  theme_bw() + coord_fixed() + facet_wrap( ~ scenario, nrow = 3) +
  scale_colour_manual(name = "Zone", values = colour_palette) +
  labs(size = "Area available for logging (Mha)", 
       x="",y="") +
  theme(text = element_text(size = 30), 
        axis.ticks=element_blank(), axis.text.x=element_blank(), 
        axis.text.y=element_blank(),
        panel.background = element_rect(fill="white", colour = "black"), 
        panel.grid = element_blank(),
        strip.text = element_text(colour = 'black'), 
        strip.background = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 35)) +
  guides(colour=FALSE)

source("codes/splitFacet.R")
new_plots <- splitFacet(g1 + theme(legend.position = "none"))

legend_size <- as_ggplot(get_legend(g1 + theme(legend.position = "bottom")))

zone_legend <- ggplot(df_zones[-10],aes(x=vext,y=as.factor(trot),fill=zname,label=zname)) + 
  geom_raster() + scale_fill_manual(values = colour_palette)+geom_text(size=10) + 
  labs(x =expression("Extracted volume ("*m^3*ha^{-1}*")"), y = "Cutting cycle (yrs)") + 
  theme(legend.position="none", aspect.ratio=1, 
        panel.background = element_blank(), 
        plot.margin = margin(2, 2, 2, 2, "cm"),
        text = element_text(size=25), 
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) 

g_all <- ggarrange(new_plots[[1]], new_plots[[2]], new_plots[[3]],
                   new_plots[[4]], new_plots[[5]], new_plots[[6]],
                   new_plots[[7]], new_plots[[8]], zone_legend, 
                   ncol=3, nrow = 3)
ggarrange(g_all, legend_size, nrow = 2, heights = c(10,1))
ggsave("graphs/mapsScenarios.pdf", height=18, width=18)

## proportion of area per zone per strategy ##
dfAreaZone = scenariOptim[demand == current_demand, .(area = sum(areaLogging)), .(zname, scenario)]
dfAreaZone = subset(dfAreaZone, area > 1e4)  ## keep only zones that represent over 1% of the total area
dfAreaZone = merge(dfAreaZone, dfAreaZone[,.(areaTot = sum(area)), .(scenario)], by = "scenario")
dfAreaZone[, pArea := area/areaTot*100]

ggplot(dfAreaZone, aes(zname, weight = pArea, fill = zname)) + geom_bar() + facet_wrap( ~ scenario, ncol = 4) + 
  scale_fill_manual(name = "Zone", values = colour_palette) + theme(legend.position = "none") +
  geom_text(aes(label=round(pArea, digits = 1), y = pArea), vjust=-0.2)
ggsave("graphs/proportionAreaZone.pdf", height = 6, width = 12)


### associated ES costs 

col_scenarios <- c("#6495ED","#458B00", "#E5C616", "#CD6600", "#E9967A","#483D8B", "#D33B44", "#8B008B")

scenCost = subset(demandFinal, demand == current_demand & variable %in% c("timber", "carbon", "biodiv"))
scenCost$ES = factor(scenCost$variable)
levels(scenCost$ES) = c("(a) Timber","(b) Carbon","(c) Biodiversity")

g2 <- ggplot(scenCost, aes(x=scenario, fill=scenario, y=value)) + 
  geom_histogram(stat="identity") + 
  facet_grid(~ ES) + scale_fill_manual(values= col_scenarios) +
  labs(y="Variation (% initial value)", fill="Strategy", x="Strategy") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.background = element_rect(fill="white",colour="white"), 
        panel.background = element_rect(fill="white", colour = "black"),
        strip.background = element_blank(),
        panel.grid = element_blank())
g2
ggsave("graphs/costsScenario.pdf", height=4, width=7)

g2 + geom_text(aes(label = round(value, 1)))
ggsave("graphs/costsScenario_annotated.pdf", height=4, width=7)

### ES = f(timber demand, scenario) ###

levels(demandFinal$variable) <- c("Total area logged (Mha)",
                                  "Mean logging intensity (m3/ha)",
                                  "Mean cutting cycle length (yr)",
                                  "Timber variation (%)",
                                  "Carbon variation (%)", 
                                  "Biodiversity variation (%)") 
df_annotate = data.frame(variable = levels(demandFinal$variable), 
                         lttr = paste0("(", letters[1:6], ")"), 
                         side = c(rep(-Inf,3), rep(Inf,3)), 
                         hl = c(rep(-0.5,3), rep(1.5,3)), 
                         vl = c(1.5, 3, 3, 1.5, 2.5, 2.5))

legend_strategies <- as_ggplot(get_legend(g2)) 

g3 <- ggplot(demandFinal, aes(x=demand, y= value, colour=scenario)) + 
  geom_hline(data = data.frame(variable = levels(demandFinal$variable), h = c(rep(c(NA,0), each=3))),
             aes(yintercept = h), lty=2) + 
  geom_line(lwd=0.7) + #scale_colour_brewer(palette = "Set1") +
  facet_wrap( ~ variable, scale="free_y", nrow=3, dir = "v", strip.position = "left") + 
  theme(panel.background = element_rect(fill="white", colour = "black"),
        panel.grid = element_blank(), strip.background = element_blank(), 
        legend.position = "none", strip.placement = "outside") +
  geom_text(data = df_annotate, aes(x = side, y = Inf, label = lttr, hjust = hl, vjust = vl), colour = "black") +
  scale_colour_manual(values= col_scenarios) + 
  labs(x=expression("Timber production (M"*m^3*yr^{-1}*")"), y="",colour="Strategy") 
ggarrange(g3, legend_strategies, ncol = 2, widths =  c(4,1))
ggsave("graphs/increasingDemand.pdf", height=6, width=8)



##### Supplementary graphs #####

### damage - explanation diagram ###
pdf("graphs/schemaDam.pdf", height=2, width=3)
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

paletteTern = colorRampPalette(c("red","orange", "yellow", "blue"))(100)

costsTot$loss_rel
levels(costsTot$ES) <- c("(a) Timber","(b) Carbon","(c) Biodiversity")

ggtern(data=costsTot,aes(x=alphaV,y=alphaC,z=alphaB, value=loss)) + 
  stat_interpolate_tern(geom="polygon", method=lm, colour="black", formula = value~x+y,
                        n=100,aes(fill=..level..),expand=1)  +
  geom_point(size = 1) + 
  facet_wrap( ~ ES) + scale_fill_gradientn(colours = paletteTern) +
  labs(x=expression(alpha[T]), y=expression(alpha[C]), z=expression(alpha[B]), fill="ES\nvariation (%)") +
  theme(strip.background = element_blank(), legend.position = "bottom",
        strip.text = element_text(hjust = 0, size = 15))
ggsave("graphs/changingESweights.pdf", height=4, width=10)

## maps with changing demand 
scenariOptim$demand2 = paste(scenariOptim$demand, "Mm3/yr")
ggplot(subset(scenariOptim, demand != current_demand)) +
  geom_point(aes( x = long, y = lat, colour = zname, size = areaLogging/1e6)) +
  theme_bw() + coord_fixed() + facet_grid(demand2 ~ scenario) +
  scale_colour_manual(name = "Zone",values = colour_palette) +
  labs(size = "Area available for logging (^Mha)", x="", y="") +
  theme(legend.position = "top", text = element_text(size=50),
        axis.ticks=element_blank(), axis.text.x=element_blank(), 
        axis.text.y=element_blank()) +
  guides(colour=FALSE)
ggsave("graphs/mapsChangeDemand.pdf", height=40, width=40)

