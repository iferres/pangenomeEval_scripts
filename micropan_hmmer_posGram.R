#Micropan - hmmer
library(micropan)
library(parallel)
library(ulimit)

#
hmmerScan <- function( in.files, db, out.folder, threads=0, verbose=TRUE ){
  basic <- paste( "hmmscan -o hmmer3.txt --cut_ga --noali --cpu", threads )
  db.name <- rev( unlist( strsplit( db, split="/" ) ) )[1]
  for( i in 1:length( in.files ) ){
    gi <- gregexpr( "GID[0-9]+", in.files[i], extract=T )
    rname <- paste( gi, "_vs_", db.name, ".txt", sep="" )
    res.files <- dir( out.folder )
    if( !(rname %in% res.files) ){
      if( verbose ) cat( "hmmerScan: Scanning", in.files[i], "...\n" )
      command <- paste( basic, "--domtblout", file.path( out.folder, rname ), db, in.files[i]  )
      print( command )
      system( command )
      file.remove( "hmmer3.txt" )
    }
  }
}

ulimit::memory_limit(20000)
setwd("/mnt/ubi/iferres/pewitEval/posGram/")

dirs <- list.dirs(recursive = FALSE)
dirs <- grep('_gff', dirs, value = TRUE)

dir.create('micropan_hmmer_resu')

db <- '/mnt/ubi/iferres/pfam/Pfam-A.hmm'


fin <- mclapply(1:length(dirs), function(d){
  
  gffs <- list.files(path = dirs[d], pattern = 'gff$', full.names = TRUE)
  gid <- paste0("%0",nchar(length(gffs)),"d")
  gid <- paste0('GID',sprintf(gid,1:length(gffs)))
  ref <- cbind(sapply(strsplit(gffs,'/'),function(x){rev(x)[1]}), gid)
  write.csv(ref, file = 'micropan_hmmer_resu/ref_gid.csv',quote = F)
  
  df <- mclapply(1:5, function(i){
    set.seed(i)
    gfs <- sample(gffs, 25)
    out <- paste0('micropan_hmmer_resu',sub('[./]','',dirs[d]),'_out_',i,'/')
    dir.create(out)
    faas <- sapply(gfs, function(x){
      pewit:::extractSeqsFromGff3(x, 
                                  in.path = out, 
                                  keep = 'none', 
                                  write.in.path = 'aa')
      paste0(out,sub('gff$','faa', rev(strsplit(x,'/')[[1]])[1]))
    })
    
    faas <- gsub('#','_', faas)
    
    pprep <- sapply(faas, function(x){
      
      gd <- ref[which(gsub('#','_',ref[,1])==sub('faa$','gff',rev(strsplit(x,'/')[[1]])[1])),2]
      panPrep(x, GID.tag = gd, out.file = x,protein = TRUE)
      sfx <- paste0('_',gd,'.faa')
      sub('[.]faa$',sfx, x)
    })
    
    
    #out blast folder
    oblf <- paste0(out,'hmmer_out')
    dir.create(oblf)
    
    runMicropanHmmer <- function(faas, db, out.folder){
      hmmerScan(faas, db, out.folder = out.folder, threads = 0L, verbose = FALSE)
      pfam.files <- list.files(path = out.folder, full.names = TRUE)
      pfam.table <- lapply(pfam.files, function(i){
        tab <- readHmmer(i)
        tab <- hmmerCleanOverlap(tab)
        tab
      })
      pfam.table <- do.call(rbind, pfam.table)
      
      cluster.pfam <- dClust(pfam.table)
      pm.pfam <- panMatrix(cluster.pfam)
      
      return(pm.pfam)
    }
    
    
    stime <- system.time(xx <- runMicropanHmmer(pprep, db, out.folder = oblf))
    
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


saveRDS(fin, file = 'micropan_hmmer_resu/resu_posGram.RDS')
