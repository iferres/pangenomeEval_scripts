setwd("/mnt/ubi/iferres/pewitEval/genomeSize")

dir.create('roary_i95_resu')

dirs <- list.dirs(recursive = FALSE)
dirs <- grep('_gff', dirs, value = TRUE)



fin <- mclapply(1:length(dirs), function(d){
  
  gffs <- list.files(path = dirs[d], pattern = 'gff$', full.names = TRUE)
  
  
  
  df <- mclapply(1:5, function(i){
    
    set.seed(i)
    gfs <- sample(gffs, 10)
    out <- paste0('roary_i95_resu',sub('[./]','',dirs[d]),'_out_',i)
    
    roary <- paste0('roary -p 1 -cd 100 -i 95 -f ',out,' ',paste(gfs, collapse = ' '))
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


saveRDS(fin, file = 'roary_i95_resu/resu_genomeSize.RDS')

