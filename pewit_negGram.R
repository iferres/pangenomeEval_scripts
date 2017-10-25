# PEWIT
library(parallel)
library(pewit)

pfam <- '/mnt/ubi/iferres/pfam'
hmmpfam <- paste0(pfam,'/Pfam-A.hmm')
datpfam <- paste0(pfam,'/Pfam-A.hmm.dat')

setwd("/mnt/ubi/iferres/pewitEval/negGram/")

dir.create('pewit_resu')

dirs <- list.dirs()
dirs <- grep('_gff', dirs, value = TRUE)


fin <- mclapply(1:length(dirs), function(d){
  
  gffs <- list.files(path = dirs[d], pattern = 'gff$', full.names = TRUE)
  
  
  
  df <- mclapply(1:5, function(i){
    set.seed(i)
    gfs <- sample(gffs, 25)
    out <- paste0('pewit_resu',sub('[./]','',dirs[d]),'_out_',i)
    stime <- system.time( p <- pangenome(gffs = gfs,
                                         hmmPfam = hmmpfam,
                                         datPfam = datpfam,
                                         n_threads = 1L,
                                         dir_out = out, 
                                         writeFfns = FALSE,
                                         writeFastas = FALSE,
                                         pmOutfileType = 'npara',
                                         alignCore = FALSE,
                                         accuAli = FALSE,
                                         coreLevel = 0.95))
    
    o <- NULL
    
    o[1] <- list(gfs)
    o[2] <- attr(p,'output')
    o[3] <- attr(p,'ncds') 
    o[4] <- attr(p,'norgs') 
    o[5] <- attr(p,'nclust') 
    o[6] <- attr(p, 'ncore100') 
    o[7] <- attr(p, 'ncore95') 
    o[8] <- attr(p, 'nsingles') 
    o[9] <- attr(p, 'naccs')
    o[10] <- stime[[3]]
    
    names(o) <-  c("Orgs", "OutDir", "Num_CDS", "Num_Orgs", 
                   "Num_Clusters", "Core", "SCore", "Singles", 
                   "Accs", "Sys_time")
      
    return(o)
    
    
  }, mc.cores = 5)
  
  resu <- do.call(rbind, df)
  
  return(resu)
  
  
}, mc.cores = 3)

saveRDS(fin, file='pewit_resu/pewit_posGram.RDS')
