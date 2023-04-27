  mrna.count.DESeqnormalized <- readRDS('data/yeast/mrna.count.DESeqnormalized.rds')
  dim(mrna.count.DESeqnormalized)
  
  library(pbmcapply)
  
  sizefactors <- readRDS('data/yeast/mrna.count.DESeq.sizeFactors.rds')
  head(sizefactors)
  
  s_hat <- psych::geometric.mean(sizefactors)
  
  row <- 2
  me.var.log2 <- pbmclapply(X = 1:nrow(mrna.count.DESeqnormalized[ , ]), FUN = function(row) {
    dummy.bootstraps <- 300  
    dummy.samples <- ncol(mrna.count.DESeqnormalized)
    dummy.mu <- mean(as.numeric(mrna.count.DESeqnormalized[row, ]))
    
    me.bootstraps.me.var.log2 <- sapply(1:dummy.bootstraps, FUN = function(i) {
      dummy.count <- rpois(n = dummy.samples, lambda = s_hat*dummy.mu)
      dummy.count.log2 <- log2(dummy.count + 0.5)
      dummy.var.log2 <- var(dummy.count.log2)  
      return(dummy.var.log2)
    })
    
    ggplot() + geom_histogram(aes(me.bootstraps.me.var.log2))
    
    me.bootstraps.me.var.log2.mu <- mean(me.bootstraps.me.var.log2)
    result <- list('gene.id' = rownames(mrna.count.DESeqnormalized)[row] , 'total.var.log2' = var(log2(as.numeric(mrna.count.DESeqnormalized[row, ] + 0.5))), 
                   'me.var.log2.estimated' = me.bootstraps.me.var.log2.mu)
    
    cat(unlist(result), '\n')
    
    return(result)
  }, mc.cores = 20, ignore.interactive = F)  
  
  me.var.log2 <- do.call(rbind, me.var.log2)  
  me.var.log2 <- as.data.frame(me.var.log2)
  me.var.log2$gene.id <- as.character(me.var.log2$gene.id)    
  me.var.log2$total.var.log2 <- as.numeric(me.var.log2$total.var.log2)
  
  me.var.log2$me.var.log2.estimated <- as.numeric(me.var.log2$me.var.log2.estimated)
  me.var.log2$me.to.total.var <- me.var.log2$me.var.log2.estimated / me.var.log2$total.var.log2 * 100
  
  
  saveRDS(me.var.log2, file = 'output/estimated.me.var.log2.DESeq.using.dummy.rds')
  
  


