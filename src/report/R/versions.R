get_synder_version <- function(){
  system2('synder', args=list('-v'), stdout=TRUE)
}

get_fagin_version <- function(config){
  # TODO: Make this portable
  system2('cat', args=list(file.path(config$d_home, 'VERSION')), stdout=TRUE) 
}
