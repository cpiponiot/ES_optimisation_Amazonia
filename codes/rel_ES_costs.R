### long, lat : coordinates of all cells of zone "zone" 
### area: PPF area of each cell 

rel_ES_costs <- function(zone, long, lat, areaLogging, cost) {
  
  if (length(areaLogging) != length(long) | length(lat) != length(long) )
    stop("long, lat and areaLogging must be of the same length")
  
  if (length(cost) != 1 | !(cost %in% c("timber","carbon","biodiversity")))
    stop("Cost must be either timber, carbon or biodiversity.")
  
  df = data.table(long, lat, areaLogging, zone)
  
  if (cost == "timber") {
    maps_zones = cost_vrec
    initial_values = Mvcom0
    initial_values$val0 = initial_values$Vcom0
  } else if (cost == "carbon") {
    maps_zones = cost_carbon
    initial_values = data.frame(coordinates(grd), val0 = grd$acs)
  } else if (cost == "biodiversity") {
    maps_zones = cost_diversity
    initial_values = data.frame(coordinates(grd), val0 = grd$mammals + grd$amphi)
  }
  
  df = df[,.(cost_cell = raster::extract(maps_zones[[zone]], cbind(long,lat)), 
             long, lat, areaLogging), .(zone)]
  
  df$areaTot = raster::extract(raster(grd, "areaForest"), df[,c("long","lat")])
  
  # get initial volume per pixel
  df = merge(df, initial_values, by=c("long","lat"))
  
  return(df[,(-sum(cost_cell*areaLogging)/sum(val0*areaTot))*100])
  
}