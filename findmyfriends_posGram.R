# FindMyFriends
library(parallel)
library(FindMyFriends)

setwd("/mnt/ubi/iferres/pewitEval/posGram/")

dir.create('findmyfriends_resu')

dirs <- list.dirs()
dirs <- grep('_gff$', dirs, value = TRUE)

run_fmf <- function(faas){
  
  pan <- pangenome(faas, 
                   translated = TRUE, 
                   geneLocation = 'prodigal', 
                   lowMem = F)
  cdh <- cdhitGrouping(pan)
  nspl <- neighborhoodSplit(cdh)
  pm <- as(nspl, 'matrix')
  pm[which(pm>1)] <- 1L
  pm
}


fin <- mclapply(1:length(dirs), function(d){
  
  gffs <- list.files(path = dirs[d], pattern = 'gff$', full.names = TRUE)
  
  
  
  df <- mclapply(1:5, function(i){
    set.seed(i)
    gfs <- sample(gffs, 25)
    out <- paste0('findmyfriends_resu/',sub('[./]','',dirs[d]),'_out_',i,'/')
    dir.create(out)
    
    faas <- sapply(gfs, function(x){
      pewit:::extractSeqsFromGff3(x, 
                                  in.path = out, 
                                  keep = 'none', 
                                  write.in.path = 'aa')
      paste0(out,sub('gff$','faa', rev(strsplit(x,'/')[[1]])[1]))
    })
    
    #Reformat faas headers
    sapply(gfs, function(x){
      rl <- readLines(x)
      tb <- pewit:::extractGffTable(rl)
      tb <- tb[which(tb$Type=='CDS'), ]
      nhe <- apply(tb, 1, function(y){
        
        paste0('>', 
               y[1], 
               '_', 
               y[2], 
               ' # ', 
               y[7], 
               ' # ', 
               y[8],
               ' # ', 
               ifelse(y[9]=='+', 1, -1))
        
      })
      ff <- grep(paste0(sub('gff$', 'faa', basename(x)), '$'), faas, value = T)
      rl <- readLines(ff)
      rl[grep('^>', rl)] <- nhe
      writeLines(rl, con = ff)
      
    })
    
    stime <- system.time( p <- run_fmf(faas))
    
    o <- NULL
    
    o[1] <- list(gfs)
    o[2] <- out
    o[3] <- sum(p) 
    o[4] <- dim(p)[2] 
    o[5] <- dim(p)[1]
    o[6] <- length(which(apply(p, 1, function(x){all(x>=1L)})))
    o[7] <- length(which(rowSums(p)>=round((ncol(p))*0.95))) 
    o[8] <- length(which(rowSums(p)==1)) 
    o[9] <- length(which(rowSums(p)<round(ncol(p)*0.95)))-o[[8]][1]
    o[10] <- stime[[3]]
    
    names(o) <-  c("Orgs", "OutDir", "Num_CDS", "Num_Orgs", 
                   "Num_Clusters", "Core", "SCore", "Singles", 
                   "Accs", "Sys_time")
    
    return(o)
    
    
  }, mc.cores = 5)
  
  resu <- do.call(rbind, df)
  
  return(resu)
  
  
}, mc.cores = 2)

saveRDS(fin, file='findmyfriends_resu/findmyfriends_posGram.RDS')
