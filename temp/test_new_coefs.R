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

#### study region & maps ####

load("data/maps.Rdata")
load("data/grd.Rdata")
source("codes/proportion_areas.R")


#### ES costs & prioritisation ####

coeffs_balanced = c(0.4,0.3,0.3)
new_name = "433"

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


source("codes/solve_problem.R")

source("codes/createTargetsCostsFeatures.R")

solZones <- rbind(solve_problem(costMaps$timber, featureMaps$base, targetsList$base, "volume", timeLim = 1000), 
                  solve_problem(costMaps$carbon, featureMaps$base, targetsList$base, "carbon", timeLim = 1000),
                  solve_problem(costMaps$biodiversity, featureMaps$base, targetsList$base, "biodiv", timeLim = 1000),
                  solve_problem(costMaps$balanced, featureMaps$base, targetsList$base, "balanced", timeLim = 1000),
                  solve_problem(costMaps$balanced, featureMaps$medium, targetsList$base, "medium", timeLim = 1000), 
                  solve_problem(costMaps$balanced, featureMaps$STY, targetsList$STY, "STY", timeLim = 1000),
                  solve_problem(costMaps$sharing, featureMaps$sharing, targetsList$base, "sharing", timeLim = 1000), 
                  solve_problem(costMaps$sharing, featureMaps$sharingSTY, targetsList$STY, "sharingSTY", timeLim = 1000))
solZones$demand = 35
scenariOptim$zone = as.numeric(scenariOptim$zone)
scenariOptim$demand = as.numeric(scenariOptim$demand)
scenariOptim$long = as.numeric(scenariOptim$long)
scenariOptim$lat = as.numeric(scenariOptim$lat)
save(solZones, file = paste0("interm_results/solZones_35_", new_name, ".Rdata"))

scenariOptim = solZones

source("codes/rel_ES_costs.R")

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

g1 <- ggplot(subset(scenariOptim, demand == 35)) +
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
ggsave(paste0("graphs/mapsScenarios_", new_name, ".pdf"), height=18, width=18)

## proportion of area per zone per strategy ##
dfAreaZone = scenariOptim[demand == 35, .(area = sum(areaLogging)), .(zname, scenario)]
dfAreaZone = subset(dfAreaZone, area > 1e4)  ## keep only zones that represent over 1% of the total area
dfAreaZone = merge(dfAreaZone, dfAreaZone[,.(areaTot = sum(area)), .(scenario)], by = "scenario")
dfAreaZone[, pArea := area/areaTot*100]

ggplot(dfAreaZone, aes(zname, weight = pArea, fill = zname)) + geom_bar() + facet_wrap( ~ scenario, ncol = 4) + 
  scale_fill_manual(name = "Zone", values = colour_palette) + theme(legend.position = "none") +
  geom_text(aes(label=round(pArea, digits = 1), y = pArea), vjust=-0.2)
ggsave(paste0("graphs/proportionAreaZone_", new_name, ".pdf"), height = 6, width = 12)


### associated ES costs 

col_scenarios <- c("#6495ED","#458B00", "#E5C616", "#CD6600", "#E9967A","#483D8B", "#D33B44", "#8B008B")

scenCost = subset(demandFinal, demand == 35 & variable %in% c("timber", "carbon", "biodiv"))
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
ggsave(paste0("graphs/costsScenario_", new_name, ".pdf"), height=4, width=7)
