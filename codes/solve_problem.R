solve_problem <- function(costs, zn, targ, name, timeLim=300){
  
  # zonePenalty = diag(10)
  # for (i in 1:9) for (j in c(1:9)[-i])  zonePenalty[i,j] = 0.8
  
  prob <- problem(costs, zn) %>%
    add_min_set_objective() %>%
    add_manual_targets(targ) %>%
    add_binary_decisions() 
    
  
  sol <- try(solve(prob %>% 
                     add_rsymphony_solver(time_limit=timeLim)), 
             silent = TRUE)
  
  if (!is(sol,"try-error")) {
    
    if (!is(try(raster::extract(category_layer(sol), grd), silent = TRUE),"try-error")){
      
      ## when no problem occured
      fin_zones <- raster::extract(category_layer(sol), grd)
      
    } else {
      
      # when the time limit was too short: try doubling it, and if it doesn't work return NA
      sol <- try(solve(prob %>% add_rsymphony_solver(time_limit = 2*timeLim)), 
                 silent = TRUE)
      
      if (!is(try(raster::extract(category_layer(sol), grd), silent = TRUE),"try-error")){
        
        ## doubling the time limit was enough
        fin_zones <- raster::extract(category_layer(sol), grd)
        
      } else {
        ## the time limit was still too short: return NAs
        fin_zones <- rep(NA, length(grd$pHarv))
      }
    }
  } else {
    ## when the problem is infeasible: solver returns error
    fin_zones <- rep(NA, length(grd$pHarv))
  }
  
  return(data.table(coordinates(grd), zone = fin_zones, scenario = name))
}