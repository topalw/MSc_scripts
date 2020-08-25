private.allele.count <- function(start, fin, path_to_gds){
  require(SNPRelate)
  gen.tab <- snpgdsGetGeno(path_to_gds)
  count <- 0
  tmpfin <- fin +1
  tmpstart <- start - 1
  for (i in 1:dim(gen.tab)[2]){
    if (start == 1){
      rest_column <- gen.tab[tmpfin:dim(gen.tab)[1], i]
      rest_nas <- sum(is.na(rest_column))
      rest_total <- 2*length(rest_column) - 2*rest_nas
    } else {
      rest_column <- c(gen.tab[1:tmpstart,i], gen.tab[tmpfin:dim(gen.tab)[1],i])
      rest_nas <- sum(is.na(rest_column))
      rest_total <- 2*length(rest_column) - 2*rest_nas
    }
    if( sum(rest_column,na.rm=T)  == 0 | sum(rest_column,na.rm=T) == rest_total){
    count <- count + 1
    }
  }
  return(paste('Private alleles =', count, sep = ''))
}

