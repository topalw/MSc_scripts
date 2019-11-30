#From gds to PCA analysis (first 4 components + parallel coordinate plot) and distance PCoA 
  require(SNPRelate)
  require(MASS)
  require(ggplot2)
  #path parameters to two files
  setwd('')                   #files will be saved here
  path_to_gds <- ''           #where is the gds?
  ind_pop <- ''               #tab-delimited file with col1 = individual names, col2= population names
  
  #open files and get samples, populations
  genofile <- snpgdsOpen(path_to_gds)
  popinfo <- read.delim(ind_pop,sep='\t')
  samples <- popinfo$V1
  pops <- popinfo$V2
  
  #do the pca
  pca <- snpgdsPCA(genofile, autosome.only = F)
  
  #variance explained by each PC
  variance_explained <- head(round(pca$varprop*100,3))
  variance_explained
  
  #parcoord plot
  tmp <- pops[match(pca$sample.id,samples)]
  parcoord(pca$eigenvect[,1:10],col=tmp)
  title(main='')
  legend(legend= levels(tmp), col = 1:nlevels(tmp))
  
  #create dataframe with 4 first PCs as in SNPRelate tutorial
  tab <- data.frame(sample.id=pca$sample.id, 
    pop=factor(pops)[match(pca$sample.id, samples)],
    EV1=pca$eigenvect[,1],
    EV2=pca$eigenvect[,2],
    EV3=pca$eigenvect[,3],
    EV4=pca$eigenvect[,4],
    stringsAsFactors = FALSE)
  
#plotting 
  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black",size = 2),
                   panel.grid.major = element_blank(),
                   legend.title = element_blank(),
                   panel.grid.minor = element_line(size=0.4, colour='grey'))
  
#Loop plot plots PC1 - PC2:4
for (i in 2:4){
  tmp <- paste('EV', i, sep='') # column to take
  pca_plot <- ggplot(tab, aes(x=EV1, y=tmp, colour=pop)) + geom_point(size=6) + mytheme +
    geom_text(label=tab$sample.id, size=4, hjust=1.5,vjust=1.5)
  
  new_title <- paste('corrPCA_PC1-PC', i, '.jpeg',sep='') #file title
  #file
  jpeg(file = new_title, width = 809, height = 1024)
  #axis labels with variance explained
  xname <- paste('PC1-', variance_explained[1], sep = '')
  yname <- paste('PC', i, '-', variance_explained[i], sep = '')
  
  pca_plot + labs(y = yname, x = xname, colour='Populations')
  dev.off()
}
  
#Distance PCoA
  Inbreeding <- snpgdsIndivBeta(genofile, autosome.only = F)
  b <- snpgdsIndivBetaRel(Inbreeding, min(Inbreeding$beta))
  tmp <- b$beta[order(b$sample.id), order(b$sample.id)]     #ordered beta table
  #from beta to distances as in Hierfstat
  diag(tmp) <- NA
  minx <- min(tmp, na.rm=T)
  distmat <- 1- (tmp-minx)/(1-minx)
  diag(distmat) <- 0
  #do PCoA and make data frame with sampleid and pop code
  v1 <- cmdscale(distmat, k=12)
  pca.dist <- as.data.frame(v1)
  pca.dist$sample.id <- b$sample.id
  pca.dist$pop <- pops[order(pops)]
  
#plotting
  
for(i in 2:4){
  tmp <- paste(V,i,sep = '')
  dist_plot <- ggplot(pca.dist, aes(x=V1, y=tmp, colour=pca.dist$pop)) + geom_point(size=6)+ mytheme + 
    geom_text(label=pca.dist$sample.id, size=4, hjust=1.5,vjust=1.5)
  yname <- paste('PC', i, sep = '')
  file_name <- paste('distPCoA_PC1-PC', i, '.jpeg', sep='')
  
  jpeg(file=file_name, width = 809, height = 1024)
  dist_plot + labs(y=yname , x='PC1', colour='Populations') 
  dev.off
}
  
