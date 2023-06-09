#Generating the sample sets for yeast analysis
sub.sample.deseq.v1 <- function() { 
  rm(list = ls())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  library(dplyr)
  
  ver <- 1.0
  
  subset.size <- c(600.3)
  runs <- 4
  
  mc.cores <- 50
  
  # this cutoff for filtering trios on noise to total var
  cutoff <- 0.4
  
  path.dir <- paste0(as.character(subset.size), '/dna.and.expression/')
  print(path.dir)
  dir.create(path.dir, showWarnings = FALSE)
  
  trios.all <- readRDS(paste0(path.dir, 'findr.trios.rds'))
  print(str(trios.all)); print(colnames(trios.all)); print(dim(trios.all));
  exp.log2.deseq.residuals <- readRDS(paste0(path.dir,'findr.exp.log2.deseq.residual.rds') )
  print(colnames(exp.log2.deseq.residuals)[1:7]); print(dim(exp.log2.deseq.residuals));
  genotype <-readRDS(paste0(path.dir,'findr.genotypes.binary.strongest.eqtl.rds'))
  print(colnames(genotype)[1:5]); print(dim(genotype));
  
  all(trios.all$L.id %in% colnames(genotype))
  all(trios.all$A.id %in% colnames(exp.log2.deseq.residuals))
  all(trios.all$B.id %in% colnames(exp.log2.deseq.residuals))
  
  colnames(exp.log2.deseq.residuals) <- sub(colnames(exp.log2.deseq.residuals), pattern = '[.]', replacement = '-')
  all(trios.all$A.id %in% colnames(exp.log2.deseq.residuals))
  all(trios.all$B.id %in% colnames(exp.log2.deseq.residuals))
  
  estimated.me.var.log2.DESeq <- readRDS('estimated.me.var.log2.DESeq.using.dummy.rds')
  estimated.me.var.log2.DESeq$gene.id <- sub('[.]', '-', estimated.me.var.log2.DESeq$gene.id)
  
  str(estimated.me.var.log2.DESeq); dim(estimated.me.var.log2.DESeq)  
  trios.all$A.id.me.var <- estimated.me.var.log2.DESeq$me.var.log2.estimated[match(trios.all$A.id, 
                                                                               estimated.me.var.log2.DESeq$gene.id)]
  trios.all$B.id.me.var <- estimated.me.var.log2.DESeq$me.var.log2.estimated[match(trios.all$B.id,
                                                                               estimated.me.var.log2.DESeq$gene.id)]
  
  str(trios.all); sapply(trios.all, FUN = function(x) {
    return(any(is.na(x)))
  })
  
  getcorrelation <- function(L, A, B) {
    return(c('p.LA' = cor.test(L, A)$p.value, 'p.LB' = cor.test(L, B)$p.value, 'p.AB' = cor.test(A, B)$p.value))
  }
  
  library(pbmcapply)
  library(cit)
  source('../../scripts/modelstage2/adj_cit1_4.yeast.R')
  
  # trios.all <- trios.all[1:300, ]
  
  trios.all <- trios.all[trios.all$A.id != trios.all$B.id, ]
  
  runs.i <- 1
  for (runs.i in 1:runs) {
    print(runs.i)
    indices <- sample(1:nrow(exp.log2.deseq.residuals), size = subset.size, replace = F)
    exp.log2.deseq.residuals.subset <- exp.log2.deseq.residuals[indices, ]
    genotype.subset <- genotype[indices, ]
    
    all(rownames(genotype.subset) == rownames(exp.log2.deseq.residuals.subset))
    
    exp.log2.deseq.residuals.subset.var <- apply(exp.log2.deseq.residuals.subset, 2, var)
    estimated.me.var.log2.DESeq$total.var.log2.residual <- exp.log2.deseq.residuals.subset.var[estimated.me.var.log2.DESeq$gene.id]
    estimated.me.var.log2.DESeq$me.to.total.var.residual <- estimated.me.var.log2.DESeq$me.var.log2.estimated / estimated.me.var.log2.DESeq$total.var.log2.residual
    
    trios.all$A.id.me.var.to.totalvar <- estimated.me.var.log2.DESeq$me.to.total.var.residual[
      match(trios.all$A.id, estimated.me.var.log2.DESeq$gene.id)]
    
    trios.all$B.id.me.var.to.totalvar <- estimated.me.var.log2.DESeq$me.to.total.var.residual[
      match(trios.all$B.id, estimated.me.var.log2.DESeq$gene.id)
    ]
    
    print(any(sapply(trios.all, function(x) return(any(is.na(x))))))
    
    
    trios <- trios.all %>% filter(A.id.me.var.to.totalvar >= cutoff | B.id.me.var.to.totalvar >= cutoff)
    
    # trios <- trios[1:300, ]
    
    i <- 1
    result.cit.ecit <- pbmclapply(X = 1:nrow((trios[ , ])), FUN = function(i) {
      if(i %% 500 == 0) {
        print(trios[i, ])
      }
      L <- as.integer(genotype[, trios$L.id[i]])
      A <- as.numeric(exp.log2.deseq.residuals[, trios$A.id[i]])
      B <- as.numeric(exp.log2.deseq.residuals[, trios$B.id[i]])
      v_eA <- as.numeric(trios$A.id.me.var[i])
      v_eB <- as.numeric(trios$B.id.me.var[i])
      
      cit.AB.result <- cit.cp(L = L, G = A, T = B)
      # cit.BA.result <- cit.cp(L = L, G = B, T = A)
    
      ecit.AB.result <- get_adj_cit_pvals(L = L, Gp = A, Tp = B, v_eG = v_eA, v_eT = v_eB, 
                                       bootstrap = 1000, resampl = 100, rseed = 1)
      
      cor.AB.result <- getcorrelation(L, A, B)
      
      # ecit.BA.result <- get_adj_cit_pvals(L = L, Gp = B, Tp = A, v_eG = v_eB, v_eT = v_eA,
      #                                     bootstrap = 1000, resampl = 100, rseed = 1)
       
      #return(c(cit.AB = cit.AB.result, cit.BA = cit.BA.result, ecit.AB = ecit.AB.result, ecit.BA = ecit.BA.result))
      
      return(c(cit.AB = cit.AB.result, ecit.AB = ecit.AB.result, cor.AB = cor.AB.result))
    }, mc.cores = mc.cores, ignore.interactive = F)  
    
    # print(warnings())
    
    trios.run.i <- cbind(trios, t(as.data.frame(result.cit.ecit)))
    rownames(trios.run.i) <- NULL
    # print(str(trios.run.i))
    path.output.cit.ecit <- paste0(path.dir, paste0( 'findr.result.cit.ecit.runs.i', runs.i,'.','ver', ver,'.rds'))
    print(path.output.cit.ecit)
    path.output.samples <- paste0(path.dir, 'samples.', runs.i,'.','ver', ver,'.rds')
    print(path.output.samples)
    saveRDS(rownames(genotype)[indices], file = path.output.samples)
    saveRDS(trios.run.i, file = path.output.cit.ecit)
    
  }
  
}

