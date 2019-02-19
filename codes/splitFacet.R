splitFacet <- function(x){
  facet_vars <- names(x$facet$params$facets)         
  x$facet    <- ggplot2::ggplot()$facet              
  datasets   <- split(x$data, x$data[,facet_vars,with = F])    
  new_plots  <- lapply(1:length(datasets),function(i) { 
    x$data <- datasets[[i]]
    x$labels$title <- paste0("(",letters[i], ") ", unlist(datasets[[i]][1,facet_vars,with=FALSE]))
    x})
}    