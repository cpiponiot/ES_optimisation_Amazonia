library(prioritizr)
library(data.table)
library(rstan)
# library(sp)
# library(rgdal)
# library(raster)


## get all maps

source("codes/getMaps.R")

## get parameters

# source("codes/getParams.R")
load("pars.Rdata")
rm(list = setdiff(ls(),c("grd","pars")))

## define zones 

grd$pu = 1:nrow(grd)
df_zones <- data.table(zone = 1:10, 
                       zname = c("LS","LM","LL","MS","MM","ML","HS","HM","HL","NL"),
                       vext = c(rep(10*1:3, each=3),0), 
                       trot = c(rep(c(15,30,65), 3),1))

# # # # # # # # # # # # # # # # # # # # # 

source("codes/solve_problem.R")

# combination of costs
cost_comb = expand.grid(alphaC = seq(0,1,0.1), alphaB = seq(0,1,0.1))
cost_comb = subset(cost_comb, alphaC + alphaB <= 1)
cost_comb$alphaV = round(1 - cost_comb$alphaC - cost_comb$alphaB, digits=1)

## predictions 
# op = list.files("outputs")
# i1 <- as.numeric(gsub("changeDemand|_|.Rdata", "", op[grep("changeDemand", op)]))
# i2 <- as.numeric(gsub("changeCost|_|.Rdata", "", op[grep("changeCost", op)]))
# i_left <- setdiff(1:nrow(pars), c(i1,i2))

for (i in 1:nrow(pars)){
  
  pred <- data.table(expand.grid(pu = grd$pu, zone = 1:10))
  pred <- cbind(pred, pars[i])
  pred = merge(pred, data.frame(grd), by="pu")
  pred = merge(pred, df_zones, by="zone")
  
  source("codes/mapsVolume.R")
  source("codes/mapsCarbon.R")
  # biodiversity 
  cost_diversity <- stack(sapply(df_zones$vext, function(x) {
    grd$rich_mamm = apply(cbind(grd$mammals*(pars[i]$slopeMam*x),0), 1, max)
    grd$rich_amph = apply(cbind(grd$amphi*(pars[i]$slopeAmp*x),0), 1, max)
    return(crop(raster(grd, 'rich_mamm') + raster(grd, 'rich_amph'),extent(cost_carbon)))
  }))
  
  source("codes/createTargetsCostsFeatures.R")
  
  ## costs ##
  df_cost = lapply(1:4, function(j){
    X = list(cost_carbon, cost_diversity, cost_vrec, stack(feat_prod))[[j]]
    names(X) = 1:10
    name = c("carbon","biodiversity","timber","vextReal")[j]
    dfC = subset(data.table(as.data.frame(X, xy=TRUE)), !is.na(X1))
    colnames(dfC) = c("long","lat",1:10)
    dfC = melt(dfC, id.vars = c("long", "lat"), variable.name = "zone", value.name = "cost")
    dfC$ES = name
    return(dfC)
  })
  df_cost = do.call(rbind, df_cost)
  df_cost = dcast(df_cost, long+lat+zone~ES, value.var = "cost")
  df_cost$zone = as.numeric(df_cost$zone)
  df_cost = merge(df_cost, Mvcom0, by=c("long","lat"))
  
  rm(list = setdiff(ls(),c("grd","pars","costMaps","featureMaps","areaIFL",
                           "targetsList","solve_problem","cost_comb","i","df_zones", "df_cost")))
  
  library(parallel)
  
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores)
  
  clusterEvalQ(cl, library(prioritizr))
  clusterEvalQ(cl, library(data.table))
  
  clusterExport(cl, varlist = c("costMaps","featureMaps","areaIFL","targetsList","solve_problem","grd","cost_comb"))
  
  changeDemand <- parSapply(cl, seq(10,80,10)*1e6, function(TD){
    
    for (y in 1:length(targetsList)) 
      targetsList[[y]]$target[1] <- TD
    
    ## increase time limit for sty-landsharing 
    TLim = c(rep(1000,5), 2000, 3000, 5000)[TD*1e-7]
    
    ## carbon optimization
    solZones <- rbind(solve_problem(costMaps$carbon, featureMaps$base, targetsList$base, "carbon", timeLim = 1000),
                      solve_problem(costMaps$timber, featureMaps$base, targetsList$base, "volume", timeLim = 1000),
                      solve_problem(costMaps$balanced, featureMaps$base, targetsList$base, "balanced", timeLim = 1000),
                      solve_problem(costMaps$balanced, featureMaps$STY, targetsList$STY, "STY", timeLim = 1000),
                      solve_problem(costMaps$sharing, featureMaps$sharing, targetsList$base, "sharing", timeLim = 1000), 
                      solve_problem(costMaps$sharing, featureMaps$sharingSTY, targetsList$STY, "sharingSTY", timeLim = TLim))
    solZones$demand = TD/1e6
    return(solZones)
    
  })
  
  changeDemand = lapply(1:dim(changeDemand)[2], function(j) data.table(do.call(cbind, changeDemand[,j])) )
  changeDemand = do.call(rbind, changeDemand)
  changeDemand$zone = as.numeric(changeDemand$zone)
  changeDemand$demand = as.numeric(changeDemand$demand)
  changeDemand$long = as.numeric(changeDemand$long)
  changeDemand$lat = as.numeric(changeDemand$lat)  
  changeDemand = merge(changeDemand,df_cost, by=c("long","lat","zone"))
  
  if (pars$i[i] == "maxL")  i <- "maxL"
  
  save(changeDemand, file=paste("outputs/changeDemand_", i, ".Rdata", sep=""))
  
  changeCost <- parSapply(cl, 1:nrow(cost_comb), function(x){
    
    costsWeighted <- cost_comb[x,1] * costMaps$carbon + cost_comb[x,2] * costMaps$biodiversity + cost_comb[x,3] * costMaps$timber 
    
    ##  optimization
    solZones <- solve_problem(costsWeighted, featureMaps$base, targetsList$base, x, timeLim = 600)
    solZones$zone[solZones$zone==0] <- 10
    solZones = data.table(solZones, alphaC = cost_comb$alphaC[x], 
                          alphaB = cost_comb$alphaB[x], alphaV = cost_comb$alphaV[x])
    return(solZones)
    
  })
  
  stopCluster(cl)
  
  changeCost = lapply(1:dim(changeCost)[2], function(j) data.table(do.call(cbind, changeCost[,j])) )
  changeCost = do.call(rbind, changeCost)
  changeCost = merge(changeCost, df_cost, by=c("long","lat","zone"))
  
  save(changeCost, file=paste("outputs/changeCost_", i, ".Rdata", sep=""))
  
}

# rm(list = setdiff(ls(),c("grd","pars")))