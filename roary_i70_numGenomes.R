library(parallel)


setwd("/mnt/ubi/iferres/pewitEval/numGenomes")

dir.create('roary_i70_resu')

count <- function(x){
  x <- read.csv(x, stringsAsFactors = FALSE)
  rn <- as.character(x$Gene)
  x <- x[ , 15:ncol(x)]
  rownames(x) <- rn
  
  xx <- apply(x, 2, function(i){
    vapply(i, function(j){ length(strsplit(j, '\t')[[1]])}, 1L)
  })
  
  sum(xx)
}


sq <- seq(10,180,10)

gffs <- list.files(path = 'gffs_cfetus', pattern = 'gff$', full.names = TRUE)

fin <- mclapply(sq, function(d){
  
  
  df <- mclapply(1:10, function(i){
    set.seed(i)
    gfs <- sample(gffs, d)
    out <- paste0('roary_i70_resu/',sub('[./]','',dirs[d]),'_out_',i)
    
    roary <- paste0('roary -p 1 -cd 100 -i 70 -f ',out,' ',paste(gfs, collapse = ' '))
    stime <- system.time(system(roary))
    
    xx <- read.csv(paste0(out,'/gene_presence_absence.Rtab'),header = T, sep = '\t')
    
    
    o <- NULL
    
    o[1] <- list(gfs)
    o[2] <- out
    o[3] <- count(paste0(out,'/gene_presence_absence.csv')) 
    o[4] <- ncol(xx)-1 
    o[5] <- nrow(xx)
    o[6] <- length(which(apply(xx[,-1],1,function(x){all(x==1L)})))
    o[7] <- length(which(rowSums(xx[,-1])>=round((ncol(xx)-1)*0.95))) 
    o[8] <- length(which(rowSums(xx[,-1])==1)) 
    o[9] <- length(which(rowSums(xx[,-1])<round(ncol(xx[,-1])*0.95)))-o[[8]][1]
    o[10] <- stime[[3]]
    
    names(o) <-  c("Orgs", "OutDir", "Num_CDS", "Num_Orgs", 
                   "Num_Clusters", "Core", "SCore", "Singles", 
                   "Accs", "Sys_time")
    
    
    return(o)
    
    
  }, mc.cores = 5)
  
  resu <- do.call(rbind, df)
  
  return(resu)
  
  
}, mc.cores = 2)


saveRDS(fin, file = 'roary_i70_resu/resu_numGenomes.RDS')


