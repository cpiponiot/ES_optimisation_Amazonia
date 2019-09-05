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
current_demand <- 30 ## in Mm3
IFLconservation <- 0.8 ## proportion of intact forest landscapes to preserve per ecoregion
biodivRecov <- 0 ## slope of biodiversity recovery (yr^-1)

#### study region & maps ####

load("data/maps.Rdata")
load("data/grd.Rdata")
source("codes/proportion_areas.R")

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


##### Solving the optimisation problem #####

if (solveProblems){
  
  source("codes/solve_problem.R")
  
  coeffs_balanced = c(1,1,1)/3
  source("codes/createTargetsCostsFeatures.R")
  
  library(parallel)
  
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores)
  
  clusterEvalQ(cl, library(prioritizr))
  clusterEvalQ(cl, library(data.table))
  
  clusterExport(cl, varlist = c("cost_carbon","cost_vrec","df_zones", "harv", "df_val0", "grd",
                                "featureMaps","areaIFL","targetsList", "solve_problem"))
  
  scenariOptim <- parSapply(cl, c(0, 0.5, 1, 2, 5)/100, function(bR){
    
    ## add recovery to biodiversity cost
    cost_diversity <- stack(sapply(1:nrow(df_zones), function(i) {
      grd$rich_mamm = grd$mammals * ( 0.0144 * df_zones$vext[i])  ## species loss, according to Burivalova et al (2014)
      grd$rich_mamm = grd$rich_mamm * (1 - bR * df_zones$trot[i] )  ## recovery, relative to initial value & loss
      grd$rich_mamm [ grd$rich_mamm < 0 ] <- 0
      grd$rich_amph = grd$amphi * ( 0.0153 * df_zones$vext[i])  ## species loss, according to Burivalova et al (2014)
      grd$rich_amph = grd$rich_amph * (1 - bR * df_zones$trot[i] )  ## recovery, relative to initial value & loss
      grd$rich_amph[grd$rich_amph < 0] <- 0
      return(crop(raster(grd, 'rich_mamm') + raster(grd, 'rich_amph'),extent(cost_carbon)))
    }))
    
    cost_balanced <- (cost_diversity / mean(df_val0$B0) + cost_carbon / mean(df_val0$C0) + cost_vrec / mean(df_val0$V0)) * harv
    
    ## optimization
    solZones <- rbind(solve_problem(cost_diversity * harv, featureMaps$base, targetsList$base, "biodiv", timeLim = 1000),
                      solve_problem(cost_balanced, featureMaps$base, targetsList$base, "balanced", timeLim = 1000))
    solZones$bioRec = bR
    return(solZones)
  })
  stopCluster(cl)
  scenariOptim = lapply(1:dim(scenariOptim)[2], function(i) { 
    data.table(do.call(cbind, scenariOptim[,i])) })
  
  scenariOptim = do.call(rbind, scenariOptim)
  
  scenariOptim$zone = as.numeric(scenariOptim$zone)
  scenariOptim$long = as.numeric(scenariOptim$long)
  scenariOptim$lat = as.numeric(scenariOptim$lat)
  
  save(scenariOptim, file="outputs/scenariOptim_biodivRecov.Rdata")
  
} else {load("outputs/scenariOptim_biodivRecov.Rdata")}

### pixel area
scenariOptim = merge(scenariOptim, data.table(coordinates(grd), 
                                              areaLogging = grd$area*grd$pAreaAvail, 
                                              areaUnprotec = grd$area*grd$pAreaUnprotect), by=c("long","lat"))
scenariOptim = merge(scenariOptim, df_zones, by="zone")
scenariOptim = subset(scenariOptim, areaLogging > 0)

### scenarios names
scenariOptim$scenario <- as.factor(scenariOptim$scenario)
levels(scenariOptim$scenario) <- c("Balanced","Biodiversity")
scenariOptim$scenario <- factor(scenariOptim$scenario, levels = c("Biodiversity","Balanced") )

## get costs
scenariOptim$areaLogged = scenariOptim$areaLogging
scenariOptim$areaLogged[scenariOptim$zname=="NL"] <- 0

source("codes/rel_ES_costs.R")

bRecFinal <- scenariOptim[,.(timber = rel_ES_costs(zone, long, lat, areaLogged, "timber"), 
                               carbon = rel_ES_costs(zone, long, lat, areaLogged, "carbon"), 
                               biodiv = rel_ES_costs(zone, long, lat, areaLogged, "biodiversity")),
                            .(bioRec, scenario)]
bRecFinal = melt(bRecFinal, id.vars = c("bioRec","scenario"))


### map results

palette <- function (name, indices = c(7,5,3)) {
  RColorBrewer::brewer.pal(9, name)[indices]
}
colours <- c(as.vector(sapply(c("Blues", "Purples","Reds"), palette)), "forestgreen")
colour_palette <- c("LS"= colours[1], "LM"=colours[2], "LL"=colours[3],
                    "MS"=colours[4], "MM"=colours[5], "ML"=colours[6],
                    "HS"=colours[7], "HM"=colours[8], "HL"=colours[9], "NL"=colours[10])

scenariOptim$scenario_label = scenariOptim$scenario
levels(scenariOptim$scenario_label) = paste(levels(scenariOptim$scenario), "strategy")
scenariOptim$bioRec_label = paste0("Recovery = ", as.numeric(scenariOptim$bioRec)*100,"%/yr")

g1 <- ggplot(scenariOptim) +
  geom_point(aes(x=long,y=lat,colour=zname, size=areaLogging/1e6)) +
  theme_bw() + coord_fixed() + facet_grid( scenario_label ~ bioRec_label) +
  scale_colour_manual(name = "Zone", values = colour_palette) +
  labs(size = "Area available for logging (Mha)", 
       x="",y="") +
  theme(legend.position = "bottom", 
        text = element_text(size = 30), 
        axis.ticks=element_blank(), axis.text.x=element_blank(), 
        axis.text.y=element_blank(),
        panel.background = element_rect(fill="white", colour = "black"), 
        panel.grid = element_blank(),
        strip.text = element_text(colour = 'black'), 
        strip.background = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 30)) +
  guides(colour=FALSE)

zone_legend <- ggplot(df_zones[-10],aes(x=vext,y=as.factor(trot),fill=zname,label=zname)) + 
  geom_raster() + scale_fill_manual(values = colour_palette)+geom_text(size=10) + 
  labs(x =expression("Extracted volume ("*m^3*ha^{-1}*")"), y = "Cutting cycle (yrs)") + 
  theme(legend.position="none", aspect.ratio=1, 
        panel.background = element_blank(), 
        plot.margin = margin(2, 2, 2, 2, "cm"),
        text = element_text(size=25), 
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) 

ggarrange(g1, zone_legend, ncol = 2, widths = c(10,2))
ggsave("graphs/mapsBiodRec.pdf", height=10, width=30)
