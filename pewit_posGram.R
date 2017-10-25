library(parallel)
library(pewit)

pfam <- '/mnt/ubi/iferres/pfam'
hmmpfam <- paste0(pfam,'/Pfam-A.hmm')
datpfam <- paste0(pfam,'/Pfam-A.hmm.dat')
dir.create('pewit_resu')

dirs <- list.dirs()
dirs <- grep('_gff', dirs, value = TRUE)


fin <- mclapply(1:length(dirs), function(d){
  
  gffs <- list.files(path = dirs[d], pattern = 'gff$', full.names = TRUE)
  
  
  # df <- data.frame(Orgs=rep(NA, 5),
  #                  OutDir=rep(NA, 5),
  #                  Num_CDS=rep(NA, 5),
  #                  Num_Orgs=rep(NA, 5),
  #                  Num_Clusters=rep(NA, 5),
  #                  Core=rep(NA, 5),
  #                  SCore=rep(NA, 5),
  #                  Singles=rep(NA, 5),
  #                  Accs=rep(NA, 5),
  #                  Sys_time=rep(NA, 5))
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
    
    # df$Orgs[i] <- list(gfs)
    # df$OutDir[i] <- attr(p,'output')
    # df$Num_CDS[i] <- attr(p,'ncds') 
    # df$Num_Orgs[i] <- attr(p,'norgs') 
    # df$Num_Clusters[i] <- attr(p,'nclust') 
    # df$Core[i] <- attr(p, 'ncore100') 
    # df$SCore[i] <- attr(p, 'ncore95') 
    # df$Singles[i] <- attr(p, 'nsingles') 
    # df$Accs[i] <- attr(p, 'naccs')
    # df$Sys_time[i] <- stime[[3]]
    
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


