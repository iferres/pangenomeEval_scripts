#Micropan - blast
library(micropan)
library(parallel)

setwd("/mnt/ubi/iferres/pewitEval/genomeSize/")

dirs <- list.dirs(recursive = FALSE)
dirs <- grep('_gff', dirs, value = TRUE)

dir.create('micropan_blast_resu')




fin <- mclapply(1:length(dirs), function(d){
  
  gffs <- list.files(path = dirs[d], pattern = 'gff$', full.names = TRUE)
  gid <- paste0("%0",nchar(length(gffs)),"d")
  gid <- paste0('GID',sprintf(gid,1:length(gffs)))
  ref <- cbind(sapply(strsplit(gffs,'/'),function(x){rev(x)[1]}), gid)
  write.csv(ref, file = 'micropan_blast_resu/ref_gid.csv',quote = F)
  
  df <- mclapply(1:5, function(i){
    set.seed(i)
    gfs <- sample(gffs, 10)
    out <- paste0('micropan_blast_resu',sub('[./]','',dirs[d]),'_out_',i,'/')
    dir.create(out)
    faas <- sapply(gfs, function(x){
      pewit:::extractSeqsFromGff3(x, 
                                  in.path = out, 
                                  keep = 'none', 
                                  write.in.path = 'aa')
      paste0(out,sub('gff$','faa', rev(strsplit(x,'/')[[1]])[1]))
    })
    
    
    pprep <- sapply(faas, function(x){
      
      gd <- ref[which(ref[,1]==sub('faa$','gff',rev(strsplit(x,'/')[[1]])[1])),2]
      panPrep(x, GID.tag = gd, out.file = x,protein = TRUE)
      sfx <- paste0('_',gd,'.faa')
      sub('[.]faa$',sfx, x)
    })
    
    #out blast folder
    oblf <- paste0(out,'blast_out')
    dir.create(oblf)
    #SPECIFY CHANGED EVALUE FROM DEFAULT (1)!!
    runMicropanBlast <- function(d, i, pprep, oblf){
      njb <- paste0(d,i,sprintf('%02d',sample(1:20, 2)), collapse = '')
      blastAllAll(pprep, 
                  oblf, 
                  e.value = 1e-10,
                  threads = 1L,
                  job = as.integer(njb),
                  verbose = FALSE)
      bl <- list.files(path = oblf, full.names = TRUE)
      df <- bDist(bl, verbose = FALSE)
      bcl <- bClust(df)
      pm <- panMatrix(bcl)
      return(pm)
    }
    
    
    stime <- system.time(xx <- runMicropanBlast(d = d, 
                                                i = i, 
                                                pprep = pprep,
                                                oblf = oblf))
    
    write.table(xx, file = paste0(out,'panmatrix_micropan.tsv'), 
                quote = F, sep = '\t', row.names = T, col.names = T)
    
    o <- NULL
    
    o[1] <- list(gfs)
    o[2] <- out
    o[3] <- sum(xx) 
    o[4] <- dim(xx)[1] 
    o[5] <- dim(xx)[2]
    o[6] <- length(which(apply(xx,2,function(x){all(x==1L)})))
    o[7] <- length(which(colSums(xx)>=round((dim(xx)[1]-1)*0.95))) 
    o[8] <- length(which(colSums(xx)==1)) 
    o[9] <- length(which(colSums(xx)<round(dim(xx)[1]*0.95)))-o[[8]][1]
    o[10] <- stime[[3]]
    
    names(o) <-  c("Orgs", "OutDir", "Num_CDS", "Num_Orgs", 
                   "Num_Clusters", "Core", "SCore", "Singles", 
                   "Accs", "Sys_time")
    
    return(o)
    
    
  }, mc.cores = 5)
  
  resu <- do.call(rbind, df)
  
  return(resu)
  
  
}, mc.cores = 2)

saveRDS(fin, file = 'micropan_blast_resu/resu_genomeSize.RDS')