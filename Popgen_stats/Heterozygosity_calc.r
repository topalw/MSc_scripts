Het_calc <- function(start_row, fin_row, path_to_gds){
  require(SNPRelate)
  require(ggplot2)
  gen.tab <- snpgdsGetGeno(path_to_gds)
  start <- start_row
  fin <- fin_row #if pop is not serially distributed in file give intermediate value(s) and update column
  pre_total <- fin - start + 1
  He <- rep(NA, dim(gen.tab)[2])
  Ho <- rep(NA, dim(gen.tab)[2])
  for(i in 1:dim(gen.tab)[2]){
    column <- gen.tab[start : fin , i] #to be updated if  needed  c() + gen.tab[extrasample,i] or intervals 
    nas <- sum(is.na(column))
    total <- pre_total - nas
    if (total == 0){
      print('Found site with only NAs. Consider filtering for missing data!')
      next
    } else{
      sum0 <- length(which(column == 0))
      sum1 <- length(which(column == 1))
      sum2 <- length(which(column == 2))
      p <- (2*sum0 + sum1) / (2*total)
      q <- (2*sum2 + sum1) / (2*total)
      He[i] <- 2*p*q
      Ho[i] <- sum1 / total
    }
  }
  seHo <- sd(Ho)/sqrt(length(Ho))
  seHe <- sd(He)/sqrt(length(He))
  return(paste('Mean He:', round(mean(He),3), " ", 'Mean Ho:', round(mean(Ho),3), 
               'SE He:', round(seHe,3), ' ', 'SE Ho:', round(seHo,3)))
}
