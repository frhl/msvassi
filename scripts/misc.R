take_last <- function(x, split = '\\|'){
  vec <- unlist(strsplit(x, split = split))
  return(vec[length(vec)])
}