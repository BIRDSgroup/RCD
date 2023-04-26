library(RColorBrewer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
plot.save.manuscript <- function(plt,path, filename) {
  for (exten in c('.pdf', '.jpeg')) {
    units = "in"
    height = 4.5
    width = 5.7
    ggsave(plot = plt,
           units = units,
           height = height,
           width = width,
           filename = paste(path, filename, exten, sep = ''))
  }
  
}

# manuscript 3(B)
pAUPR.noisy <- function() {
  rm(list = ls())
  library(dplyr)

  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  runs <- c('data/yeast/600.3/dna.and.expression/3.findr.result.cit.ecit.runs.i1.rds',
            'data/yeast/600.3/dna.and.expression/3.findr.result.cit.ecit.runs.i2.rds',
            'data/yeast/600.3/dna.and.expression/3.findr.result.cit.ecit.runs.i3.rds',
            'data/yeast/600.3/dna.and.expression/3.findr.result.cit.ecit.runs.i4.rds')
  
  path.yeastract.relabeled <- 'data/yeast/RegulationMatrix_Documented_2021223_1119_1568177990.80cross3394dnaplusexpression.csv.rds'
  ground.truth <- readRDS(path.yeastract.relabeled)
  rownames(ground.truth) <- ground.truth$TF.TG
  
  aupr.k.percentage <- function(result.cit.ecit.deseq, top.k.percentage = c(0.05, 0.10, 0.30, 0.50, 1, 2, 5, 10, 20, 40, 80, 100)) {
    result.cit.ecit.deseq.cit.sorted <- result.cit.ecit.deseq %>% arrange(cit.AB.p_cit.fdr)
    result.cit.ecit.deseq.ecit.sorted <- result.cit.ecit.deseq %>% arrange(ecit.AB.adj_p_cit.fdr)
    result.cit.ecit.deseq.random.sorted <- result.cit.ecit.deseq[sample(1:nrow(result.cit.ecit.deseq), replace = F), ]
    
    
    library(PRROC)
    
    # top.k <- seq(0.05, 1.0, 0.05)
    top.k <- top.k.percentage/100
    top.k.i <- 1
    aupr.at.k <- list()
    for(top.k.i in 1:length(top.k)) {
      
      result.cit.ecit.deseq.cit.sorted.top.k <- result.cit.ecit.deseq.cit.sorted[
        1:(nrow(result.cit.ecit.deseq.cit.sorted)*top.k[top.k.i]), ]
      result.cit.ecit.deseq.ecit.sorted.top.k <- result.cit.ecit.deseq.ecit.sorted[
        1:(nrow(result.cit.ecit.deseq.ecit.sorted)*top.k[top.k.i]), ]
      
      result.cit.ecit.deseq.random.sorted.top.k <- result.cit.ecit.deseq.random.sorted[
        1:(nrow(result.cit.ecit.deseq.cit.sorted) * top.k[top.k.i]), ]
      
      pr.cit <- pr.curve(
        scores.class0 = 1 - result.cit.ecit.deseq.cit.sorted.top.k$cit.AB.p_cit.fdr[result.cit.ecit.deseq.cit.sorted.top.k$groundtruth == 1],
        scores.class1 = 1 - result.cit.ecit.deseq.cit.sorted.top.k$cit.AB.p_cit.fdr[result.cit.ecit.deseq.cit.sorted.top.k$groundtruth == 0],
        rand.compute = T)
      
      pr.ecit <- pr.curve(
        scores.class0 = 1 - result.cit.ecit.deseq.ecit.sorted.top.k$ecit.AB.adj_p_cit.fdr[result.cit.ecit.deseq.ecit.sorted.top.k$groundtruth == 1],
        scores.class1 = 1 - result.cit.ecit.deseq.ecit.sorted.top.k$ecit.AB.adj_p_cit.fdr[result.cit.ecit.deseq.ecit.sorted.top.k$groundtruth == 0],
        rand.compute = T
      )
      
      rand.aupr = sum(result.cit.ecit.deseq.random.sorted.top.k$groundtruth == 1)/nrow(result.cit.ecit.deseq.cit.sorted.top.k)
      
      aupr.at.k[[length(aupr.at.k) + 1]] <- data.frame(k = top.k.percentage[top.k.i],
                                                       trios = nrow(result.cit.ecit.deseq.cit.sorted.top.k), 
                                                       # cit.causal.pairs = sum(result.cit.ecit.deseq.cit.sorted.top.k$groundtruth == 1),
                                                       # ecit.causal.pairs = sum(result.cit.ecit.deseq.ecit.sorted.top.k$groundtruth == 1),
                                                       cit.aupr = pr.cit$auc.integral,
                                                       ecit.aupr = pr.ecit$auc.integral,
                                                       rand.aupr = rand.aupr, 
                                                       cit.aupr.fld = pr.cit$auc.integral/rand.aupr, 
                                                       ecit.aupr.fld = pr.ecit$auc.integral/rand.aupr
                                                       # cit.tp = sum(result.cit.ecit.deseq.cit.sorted.top.k$groundtruth == 1),
                                                       # ecit.tp = sum(result.cit.ecit.deseq.ecit.sorted.top.k$groundtruth == 1), 
                                                       # cit.fp = sum(result.cit.ecit.deseq.cit.sorted.top.k$groundtruth == 0),
                                                       # ecit.fp = sum(result.cit.ecit.deseq.ecit.sorted.top.k$groundtruth == 0)
      )
      
    }
    
    aupr.at.k <- do.call(rbind, aupr.at.k)
    return(aupr.at.k)

  }
  
  # top.k.percentage <- c(0.05, 0.10, 0.30, 0.50, 1, 2, 5, 10, 20, 40, 80, 100)
  top.k.percentage <- seq(5, 100, 5)
  aupr.at.k.runs <- list()
  
  error.ratio.cutoff <- 0.4
  runs.i <- 1
  for (runs.i in 1:length(runs)) {
    path.input <- runs[runs.i]
    print(path.input)
    result.cit.ecit.deseq <- readRDS(path.input)
    result.cit.ecit.deseq <- result.cit.ecit.deseq %>% filter(A.id.me.var.to.totalvar >= error.ratio.cutoff |
                                                                B.id.me.var.to.totalvar >= error.ratio.cutoff)
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$A.id != result.cit.ecit.deseq$B.id, ]
    
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$A.id %in% rownames(ground.truth), ]
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$B.id %in% colnames(ground.truth), ]
    
    result.cit.ecit.deseq <- result.cit.ecit.deseq[!is.na(result.cit.ecit.deseq$ecit.AB.adj_p_cit), ]
    result.cit.ecit.deseq <- result.cit.ecit.deseq[!is.na(result.cit.ecit.deseq$cit.AB.p_cit), ]
    
    # sapply(result.cit.ecit.deseq, function(x) return(any(is.na(x))))
    
    result.cit.ecit.deseq$ecit.AB.adj_p_cit.fdr <- p.adjust(result.cit.ecit.deseq$ecit.AB.adj_p_cit, method = 'BH')
    result.cit.ecit.deseq$cit.AB.p_cit.fdr <- p.adjust(result.cit.ecit.deseq$cit.AB.p_cit, method = 'BH')
    
    library(pbmcapply)
    
    i <- 10
    result.lapply <- pbmclapply(X = 1:nrow(result.cit.ecit.deseq[,]), FUN = function(i) {
      cat(i, paste(result.cit.ecit.deseq[i, c(1, 2, 3)]), '\n')
      result.AB <- ground.truth[result.cit.ecit.deseq$A.id[i], result.cit.ecit.deseq$B.id[i]]
      # result.BA <- ground.truth[result.cit.ecit.deseq$B.id[i], result.cit.ecit.deseq$A.id[i]]
      
      model <- result.AB
      
      return(as.integer(model))
    }, mc.cores = 30, ignore.interactive = F)
    
    result.cit.ecit.deseq$groundtruth <- do.call(c, result.lapply)
    result.cit.ecit.deseq$groundtruth <- as.character(result.cit.ecit.deseq$groundtruth)
    
    aupr.at.k <- aupr.k.percentage(result.cit.ecit.deseq, top.k.percentage)
    aupr.at.k.runs[[1 + length(aupr.at.k.runs)]] <- aupr.at.k 
    
  }
  
  aupr.at.k.runs <- do.call(rbind, aupr.at.k.runs)
  
  aupr.at.k.runs.summary <- aupr.at.k.runs %>% group_by(k) %>% 
    dplyr::summarise(cit.aupr.mu = mean(cit.aupr),
              cit.apur.sd = sd(cit.aupr),
              ecit.aupr.mu = mean(ecit.aupr), 
              ecit.aupr.sd = sd(ecit.aupr),
              rand.aupr.mu = mean(rand.aupr))
  
  aupr.at.k.runs.summary.mu <- melt(aupr.at.k.runs.summary %>% select(k, CIT = cit.aupr.mu, ECIT = ecit.aupr.mu),
        id.vars = 'k',
        variable.name = c('method'),
        value.name = 'aupr.mu')
  
  
  aupr.at.k.runs.summary.sd <- melt(aupr.at.k.runs.summary %>% select(k, CIT = cit.apur.sd, ECIT = ecit.aupr.sd),
                                    id.vars = 'k',
                                    variable.name = c('method'),
                                    value.name = 'aupr.sd')
  
  aupr.at.k.runs.summary.melt <- merge(aupr.at.k.runs.summary.mu, aupr.at.k.runs.summary.sd, 
                                  by.x = c('k', 'method'), by.y = c('k', 'method') )
  
  library(ggplot2)
  library(reshape2)
  
  
  pAUPR.noisy <- aupr.at.k.runs.summary.melt %>% ggplot(aes(x = k, y = aupr.mu)) + 
    geom_line(aes(color = method)) +
    geom_point(aes(color = method)) +
    geom_errorbar(aes(ymin = aupr.mu - aupr.sd, ymax = aupr.mu + aupr.sd, color = method)) +
    geom_abline(slope = 0, intercept = mean(aupr.at.k.runs$rand.aupr), linetype = "dashed") +
    xlab('top % predictions') +
    ylab('Avg AUPR') + 
    ggtitle(paste0('A.me.var.to.totalvar >= error.ratio.cutoff \n |B.me.var.to.totalvar >= error.ratio.cutoff \n',error.ratio.cutoff)) + 
    scale_x_continuous(breaks = seq(5,100, 10), labels = as.character(seq(5, 100, 10)), guide = guide_axis(check.overlap = T))
  
  pAUPR.noisy <- pAUPR.noisy + 
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")

  print(pAUPR.noisy)
  
  exten <- '.pdf'
  for (exten in c('.pdf', '.jpeg')) {
    ggsave(plot = pAUPR.noisy,
           units = "in",
           height = 4.5,
           width = 5.7,
           filename = paste("output/yeast.pAUPR.noisy.", 
                            error.ratio.cutoff, exten, sep = ''))
  }
  
  
  # library(ggplot2)
  # library(gridExtra)
  # 
  # grid.arrange(pAUPR.all, pAUPR.noisy, nrow = 2, ncol = 1, legend)
  
  
  
  temp <- 1

}
pAUPR.noisy()
#manuscript 3(A)
pAUPR.all <- function() {
  rm(list = ls())
  library(dplyr)

  runs <- c('data/yeast/600.4/dna.and.expression/4.findr.result.cit.ecit.runs.i1.rds',
            'data/yeast/600.4/dna.and.expression/4.findr.result.cit.ecit.runs.i2.rds',
            'data/yeast/600.4/dna.and.expression/4.findr.result.cit.ecit.runs.i3.rds',
            'data/yeast/600.4/dna.and.expression/4.findr.result.cit.ecit.runs.i4.rds')
  
  path.yeastract.relabeled <- 'data/yeast/RegulationMatrix_Documented_2021223_1119_1568177990.80cross3394dnaplusexpression.csv.rds'
  ground.truth <- readRDS(path.yeastract.relabeled)
  rownames(ground.truth) <- ground.truth$TF.TG
  
  aupr.k.percentage <- function(result.cit.ecit.deseq) {
    result.cit.ecit.deseq.cit.sorted <- result.cit.ecit.deseq %>% arrange(cit.AB.p_cit.fdr)
    result.cit.ecit.deseq.ecit.sorted <- result.cit.ecit.deseq %>% arrange(ecit.AB.adj_p_cit.fdr)
    result.cit.ecit.deseq.random.sorted <- result.cit.ecit.deseq[sample(1:nrow(result.cit.ecit.deseq), replace = F), ]
    
    
    library(PRROC)
    
    top.k <- seq(0.05, 1.0, 0.05)
    top.k.i <- 1
    aupr.at.k <- list()
    for(top.k.i in 1:length(top.k)) {
      
      result.cit.ecit.deseq.cit.sorted.top.k <- result.cit.ecit.deseq.cit.sorted[
        1:(nrow(result.cit.ecit.deseq.cit.sorted)*top.k[top.k.i]), ]
      result.cit.ecit.deseq.ecit.sorted.top.k <- result.cit.ecit.deseq.ecit.sorted[
        1:(nrow(result.cit.ecit.deseq.ecit.sorted)*top.k[top.k.i]), ]
      
      result.cit.ecit.deseq.random.sorted.top.k <- result.cit.ecit.deseq.random.sorted[
        1:(nrow(result.cit.ecit.deseq.cit.sorted) * top.k[top.k.i]), ]
      
      pr.cit <- pr.curve(
        scores.class0 = 1 - result.cit.ecit.deseq.cit.sorted.top.k$cit.AB.p_cit.fdr[result.cit.ecit.deseq.cit.sorted.top.k$groundtruth == 1],
        scores.class1 = 1 - result.cit.ecit.deseq.cit.sorted.top.k$cit.AB.p_cit.fdr[result.cit.ecit.deseq.cit.sorted.top.k$groundtruth == 0],
        rand.compute = T)
      
      pr.ecit <- pr.curve(
        scores.class0 = 1 - result.cit.ecit.deseq.ecit.sorted.top.k$ecit.AB.adj_p_cit.fdr[result.cit.ecit.deseq.ecit.sorted.top.k$groundtruth == 1],
        scores.class1 = 1 - result.cit.ecit.deseq.ecit.sorted.top.k$ecit.AB.adj_p_cit.fdr[result.cit.ecit.deseq.ecit.sorted.top.k$groundtruth == 0],
        rand.compute = T
      )
      
      rand.aupr = sum(result.cit.ecit.deseq.random.sorted.top.k$groundtruth == 1)/nrow(result.cit.ecit.deseq.cit.sorted.top.k)
      
      aupr.at.k[[length(aupr.at.k) + 1]] <- data.frame(k = top.k[top.k.i],
                                                       trios = nrow(result.cit.ecit.deseq.cit.sorted.top.k), 
                                                       # cit.causal.pairs = sum(result.cit.ecit.deseq.cit.sorted.top.k$groundtruth == 1),
                                                       # ecit.causal.pairs = sum(result.cit.ecit.deseq.ecit.sorted.top.k$groundtruth == 1),
                                                       cit.aupr = pr.cit$auc.integral,
                                                       ecit.aupr = pr.ecit$auc.integral,
                                                       rand.aupr = rand.aupr, 
                                                       cit.aupr.fld = pr.cit$auc.integral/rand.aupr, 
                                                       ecit.aupr.fld = pr.ecit$auc.integral/rand.aupr
                                                       # cit.tp = sum(result.cit.ecit.deseq.cit.sorted.top.k$groundtruth == 1),
                                                       # ecit.tp = sum(result.cit.ecit.deseq.ecit.sorted.top.k$groundtruth == 1), 
                                                       # cit.fp = sum(result.cit.ecit.deseq.cit.sorted.top.k$groundtruth == 0),
                                                       # ecit.fp = sum(result.cit.ecit.deseq.ecit.sorted.top.k$groundtruth == 0)
      )
      
    }
    
    aupr.at.k <- do.call(rbind, aupr.at.k)
    return(aupr.at.k)
    
  }
  
  
  aupr.at.k.runs <- list()
  
  runs.i <- 1
  error.ratio <- 0.0
  for (runs.i in 1:length(runs)) {
    path.input <- runs[runs.i]
    print(path.input)
    result.cit.ecit.deseq <- readRDS(path.input)
    result.cit.ecit.deseq <- result.cit.ecit.deseq %>% filter(A.id.me.var.to.totalvar >= error.ratio |
                                                                B.id.me.var.to.totalvar >= error.ratio)
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$A.id != result.cit.ecit.deseq$B.id, ]
    
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$A.id %in% rownames(ground.truth), ]
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$B.id %in% colnames(ground.truth), ]
    
    result.cit.ecit.deseq <- result.cit.ecit.deseq[!is.na(result.cit.ecit.deseq$ecit.AB.adj_p_cit), ]
    result.cit.ecit.deseq <- result.cit.ecit.deseq[!is.na(result.cit.ecit.deseq$cit.AB.p_cit), ]
    
    # sapply(result.cit.ecit.deseq, function(x) return(any(is.na(x))))
    
    result.cit.ecit.deseq$ecit.AB.adj_p_cit.fdr <- p.adjust(result.cit.ecit.deseq$ecit.AB.adj_p_cit, method = 'BH')
    result.cit.ecit.deseq$cit.AB.p_cit.fdr <- p.adjust(result.cit.ecit.deseq$cit.AB.p_cit, method = 'BH')
    
    library(pbmcapply)
    
    i <- 10
    result.lapply <- pbmclapply(X = 1:nrow(result.cit.ecit.deseq[,]), FUN = function(i) {
      cat(i, paste(result.cit.ecit.deseq[i, c(1, 2, 3)]), '\n')
      result.AB <- ground.truth[result.cit.ecit.deseq$A.id[i], result.cit.ecit.deseq$B.id[i]]
      # result.BA <- ground.truth[result.cit.ecit.deseq$B.id[i], result.cit.ecit.deseq$A.id[i]]
      
      model <- result.AB
      
      return(as.integer(model))
    }, mc.cores = 30, ignore.interactive = F)
    
    result.cit.ecit.deseq$groundtruth <- do.call(c, result.lapply)
    result.cit.ecit.deseq$groundtruth <- as.character(result.cit.ecit.deseq$groundtruth)
    
    aupr.at.k <- aupr.k.percentage(result.cit.ecit.deseq)
    aupr.at.k.runs[[1 + length(aupr.at.k.runs)]] <- aupr.at.k 
    
  }
  
  aupr.at.k.runs <- do.call(rbind, aupr.at.k.runs)
  
  library(dplyr)
  aupr.at.k.runs.summary <- aupr.at.k.runs %>% group_by(k) %>% 
    dplyr::summarise(cit.aupr.mu = mean(cit.aupr),
              cit.apur.sd = sd(cit.aupr),
              ecit.aupr.mu = mean(ecit.aupr), 
              ecit.aupr.sd = sd(ecit.aupr),
              rand.aupr.mu = mean(rand.aupr))
  
  aupr.at.k.runs.summary.mu <- melt(aupr.at.k.runs.summary %>% select(k, CIT = cit.aupr.mu, ECIT = ecit.aupr.mu),
                                    id.vars = 'k',
                                    variable.name = c('method'),
                                    value.name = 'aupr.mu')
  
  
  aupr.at.k.runs.summary.sd <- melt(aupr.at.k.runs.summary %>% select(k, CIT = cit.apur.sd, ECIT = ecit.aupr.sd),
                                    id.vars = 'k',
                                    variable.name = c('method'),
                                    value.name = 'aupr.sd')
  
  aupr.at.k.runs.summary.melt <- merge(aupr.at.k.runs.summary.mu, aupr.at.k.runs.summary.sd, 
                                       by.x = c('k', 'method'), by.y = c('k', 'method') )
  
  aupr.at.k.runs.summary.melt$k <- aupr.at.k.runs.summary.melt$k * 100
  
  library(ggplot2)
  library(reshape2)
  
  min.y <- min(mean(aupr.at.k.runs$rand.aupr), 
               min(aupr.at.k.runs.summary.melt$aupr.mu - aupr.at.k.runs.summary.melt$aupr.sd))
  max.y <- max(aupr.at.k.runs.summary.melt$aupr.mu + aupr.at.k.runs.summary.melt$aupr.sd)
  
  
  pAUPR.all <- aupr.at.k.runs.summary.melt %>% ggplot(aes(x = k, y = aupr.mu)) + 
    geom_line(aes(color = method)) +
    geom_point(aes(color = method)) +
    geom_errorbar(aes(ymin = aupr.mu - aupr.sd, ymax = aupr.mu + aupr.sd, color = method)) +
    geom_abline(slope = 0, intercept = mean(aupr.at.k.runs$rand.aupr), linetype = "dashed") +
    ylim(min.y, max.y) + 
    xlab('top k% predictions') +
    ylab('Avg AUPR') + 
    ggtitle(paste0('error_ratio >= ', error.ratio)) + 
    scale_x_continuous(breaks = seq(5,100, 10), labels = as.character(seq(5, 100, 10)), guide = guide_axis(check.overlap = T))
  
  pAUPR.all <- pAUPR.all + 
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  print(pAUPR.all)
  
  ggsave(plot = pAUPR.all,
         units = "in",
         height = 4.5,
         width = 5.7,
         filename = paste("output/pAUPR.error.ratio",error.ratio ,".pdf", sep=''))
  
  
  
}



# manuscirpt fig 3.(D, E)
precision.recall.noisy.fdr <- function() {
  rm(list = ls())
  library(dplyr)
  runs <- c('data/yeast/600.3/dna.and.expression/3.findr.result.cit.ecit.runs.i1.rds',
            'data/yeast/600.3/dna.and.expression/3.findr.result.cit.ecit.runs.i2.rds',
            'data/yeast/600.3/dna.and.expression/3.findr.result.cit.ecit.runs.i3.rds',
            'data/yeast/600.3/dna.and.expression/3.findr.result.cit.ecit.runs.i4.rds')
  
  path.yeastract.relabeled <- 'data/yeast/RegulationMatrix_Documented_2021223_1119_1568177990.80cross3394dnaplusexpression.csv.rds'
  ground.truth <- readRDS(path.yeastract.relabeled)
  rownames(ground.truth) <- ground.truth$TF.TG
  
  precision.recall.fdr <- function(result.cit.ecit.deseq, cutoff = c(0.001,0.01,0.05,0.1,0.2,0.3)) {
    a <- result.cit.ecit.deseq
    print(addmargins(table(a$groundtruth)))
    print(sum(a$groundtruth==1)/nrow(a))
    a$cit.padj = p.adjust(a$cit.AB.p_cit, method="BH")
    a$ecit.padj = p.adjust(a$ecit.AB.adj_p_cit, method="BH")
    
    
    df = data.frame(cutoff= cutoff, 
                    cit.truecalls=NA, cit.falsecalls=NA, 
                    ecit.truecalls=NA, ecit.falsecalls=NA, 
                    stringsAsFactors = F)
    for (i in seq_len(nrow(df))) 
    {  
      cf = df[i,'cutoff']
      df[i,2:5] = c(sum(a$cit.padj <= cf & a$groundtruth==1), 
                    sum(a$cit.padj <= cf & a$groundtruth==0), 
                    sum(a$ecit.padj <= cf & a$groundtruth==1), 
                    sum(a$ecit.padj <= cf & a$groundtruth==0))
    }
    df$cit.calls = df$cit.truecalls + df$cit.falsecalls
    df$ecit.calls = df$ecit.truecalls + df$ecit.falsecalls
    df$cit.prec = df$cit.truecalls / df$cit.calls
    df$ecit.prec = df$ecit.truecalls / df$ecit.calls
    print(df)
    
    return(df)
    
  }
  
  PR.runs <- list()
  
  runs.i <- 1
  for (runs.i in 1:length(runs)) {
    path.input <- runs[runs.i]
    print(path.input)
    result.cit.ecit.deseq <- readRDS(path.input)
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$A.id != result.cit.ecit.deseq$B.id, ]
    
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$A.id %in% rownames(ground.truth), ]
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$B.id %in% colnames(ground.truth), ]
    
    result.cit.ecit.deseq <- result.cit.ecit.deseq[!is.na(result.cit.ecit.deseq$ecit.AB.adj_p_cit), ]
    result.cit.ecit.deseq <- result.cit.ecit.deseq[!is.na(result.cit.ecit.deseq$cit.AB.p_cit), ]
    
    # result.cit.ecit.deseq$ecit.AB.adj_p_cit.fdr <- p.adjust(result.cit.ecit.deseq$ecit.AB.adj_p_cit, method = 'BH')
    # result.cit.ecit.deseq$cit.AB.p_cit.fdr <- p.adjust(result.cit.ecit.deseq$cit.AB.p_cit, method = 'BH')
    
    library(pbmcapply)
    
    i <- 10
    result.lapply <- pbmclapply(X = 1:nrow(result.cit.ecit.deseq[,]), FUN = function(i) {
      cat(i, paste(result.cit.ecit.deseq[i, c(1, 2, 3)]), '\n')
      result.AB <- ground.truth[result.cit.ecit.deseq$A.id[i], result.cit.ecit.deseq$B.id[i]]
      
      model <- result.AB
      
      return(as.integer(model))
    }, mc.cores = 30, ignore.interactive = F)
    
    result.cit.ecit.deseq$groundtruth <- do.call(c, result.lapply)
    result.cit.ecit.deseq$groundtruth <- as.character(result.cit.ecit.deseq$groundtruth)
    
    PR.run.i <- precision.recall.fdr(result.cit.ecit.deseq)
    PR.runs[[1 + length(PR.runs)]] <- PR.run.i
    
  }
  
  PR.runs <- do.call(rbind, PR.runs)
  
  PR.runs.call.melt <- melt(PR.runs %>% select(cutoff, CIT = cit.calls, ECIT = ecit.calls), 
       id.vars = 'cutoff',
       variable.name = 'method',
       value.name = 'calls') 
  
  PR.runs.call.melt.summary <- PR.runs.call.melt %>% group_by(cutoff, method) %>% 
    summarise(calls.mu = mean(calls), calls.sd = sd(calls))
  
  PR.runs.precision.melt <- melt(PR.runs %>% select(cutoff, CIT = cit.prec, ECIT = ecit.prec),
                                 id.vars = 'cutoff',
                                 variable.name = 'method',
                                 value.name = 'precision')

  
  PR.runs.precision.melt.summary <- PR.runs.precision.melt %>% group_by(cutoff, method) %>%
    summarise(precision.mu = mean(precision), precision.sd = sd(precision))
  
  
  PR.runs.TP.melt <- melt(PR.runs %>% select(cutoff, CIT = cit.truecalls, ECIT = ecit.truecalls),
                          id.vars = 'cutoff',
                          variable.name = 'method',
                          value.name = 'TP')
  
  PR.runs.TP.melt.summary <- PR.runs.TP.melt %>% group_by(cutoff, method) %>%
    summarise(TP.mu = mean(TP), TP.sd = sd(TP))
  
  
  library(ggplot2)
  library(reshape2)
  
  
  
  pcalls <- PR.runs.call.melt.summary %>% ggplot(aes(x = as.factor(cutoff), y = calls.mu, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = calls.mu - calls.sd, ymax = calls.mu + calls.sd), width = .2, position = position_dodge(1)) +
    labs(title="Performance on \n yeast dataset (noisy genes)", x="FDR Cutoff", y="Avg # Causal Calls", fill="")
    
  pcalls <- pcalls + 
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  
  #manuscript 3.(D) 
  pcalls
  
ggsave(plot = pcalls,
         units = "in",
         height = 4.5,
         width = 5.7,
         filename = 'output/yeast.causalcall.pdf')

  
  pprecision <- PR.runs.precision.melt.summary %>% ggplot(aes(x = as.factor(cutoff), y = precision.mu, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = precision.mu - precision.sd, ymax = precision.mu + precision.sd), width = .2, position = position_dodge(1)) +
    labs(title="Performance on \n yeast dataset (noisy genes)", x="FDR Cutoff", y="Avg Precision", fill="") 
  
  pprecision <- pprecision +
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  #manuscript 3.(E)
  pprecision
  ggsave(plot = pprecision,
         units = "in",
         height = 4.5,
         width = 5.7,
         filename = 'output/yeast.precision.noisy.pdf')
  
  
  pTP <- PR.runs.TP.melt.summary %>% ggplot(aes(x = as.factor(cutoff), y = TP.mu, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = TP.mu - TP.sd, ymax = TP.mu + TP.sd), width = .2, position = position_dodge(1)) +
    labs(title="Performance on \n yeast dataset (noisy genes)", x="FDR Cutoff", y="Avg True Positive", fill="")
  
  pTP <- pTP + 
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  
  
  pTP
  
  ggsave(plot = pTP,
         units = "in",
         height = 4.5,
         width = 5.7,
         filename = 'output/yeast.TP.pdf')
  
  
}


precision.recall.all.fdr <- function() {
  rm(list = ls())
  library(dplyr)
  
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  runs <- c('data/yeast/600.4/dna.and.expression/4.findr.result.cit.ecit.runs.i1.rds',
            'data/yeast/600.4/dna.and.expression/4.findr.result.cit.ecit.runs.i2.rds',
            'data/yeast/600.4/dna.and.expression/4.findr.result.cit.ecit.runs.i3.rds',
            'data/yeast/600.4/dna.and.expression/4.findr.result.cit.ecit.runs.i4.rds')
  
  path.yeastract.relabeled <- 'data/yeast/RegulationMatrix_Documented_2021223_1119_1568177990.80cross3394dnaplusexpression.csv.rds'
  ground.truth <- readRDS(path.yeastract.relabeled)
  rownames(ground.truth) <- ground.truth$TF.TG
  
  precision.recall.fdr <- function(result.cit.ecit.deseq, cutoff = c(0.001,0.01,0.05,0.1,0.2,0.3)) {
    a <- result.cit.ecit.deseq
    print(addmargins(table(a$groundtruth)))
    print(sum(a$groundtruth==1)/nrow(a))
    a$cit.padj = p.adjust(a$cit.AB.p_cit, method="BH")
    a$ecit.padj = p.adjust(a$ecit.AB.adj_p_cit, method="BH")
    
    
    df = data.frame(cutoff= cutoff, 
                    cit.truecalls=NA, cit.falsecalls=NA, 
                    ecit.truecalls=NA, ecit.falsecalls=NA, 
                    stringsAsFactors = F)
    for (i in seq_len(nrow(df))) 
    {  
      cf = df[i,'cutoff']
      df[i,2:5] = c(sum(a$cit.padj <= cf & a$groundtruth==1), 
                    sum(a$cit.padj <= cf & a$groundtruth==0), 
                    sum(a$ecit.padj <= cf & a$groundtruth==1), 
                    sum(a$ecit.padj <= cf & a$groundtruth==0))
    }
    df$cit.calls = df$cit.truecalls + df$cit.falsecalls
    df$ecit.calls = df$ecit.truecalls + df$ecit.falsecalls
    df$cit.prec = df$cit.truecalls / df$cit.calls
    df$ecit.prec = df$ecit.truecalls / df$ecit.calls
    print(df)
    
    return(df)
    
  }
  
  PR.runs <- list()
  
  runs.i <- 1
  for (runs.i in 1:length(runs)) {
    path.input <- runs[runs.i]
    print(path.input)
    result.cit.ecit.deseq <- readRDS(path.input)
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$A.id != result.cit.ecit.deseq$B.id, ]
    
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$A.id %in% rownames(ground.truth), ]
    result.cit.ecit.deseq <- result.cit.ecit.deseq[result.cit.ecit.deseq$B.id %in% colnames(ground.truth), ]
    
    result.cit.ecit.deseq <- result.cit.ecit.deseq[!is.na(result.cit.ecit.deseq$ecit.AB.adj_p_cit), ]
    result.cit.ecit.deseq <- result.cit.ecit.deseq[!is.na(result.cit.ecit.deseq$cit.AB.p_cit), ]
    
    # result.cit.ecit.deseq$ecit.AB.adj_p_cit.fdr <- p.adjust(result.cit.ecit.deseq$ecit.AB.adj_p_cit, method = 'BH')
    # result.cit.ecit.deseq$cit.AB.p_cit.fdr <- p.adjust(result.cit.ecit.deseq$cit.AB.p_cit, method = 'BH')
    
    library(pbmcapply)
    
    i <- 10
    result.lapply <- pbmclapply(X = 1:nrow(result.cit.ecit.deseq[,]), FUN = function(i) {
      cat(i, paste(result.cit.ecit.deseq[i, c(1, 2, 3)]), '\n')
      result.AB <- ground.truth[result.cit.ecit.deseq$A.id[i], result.cit.ecit.deseq$B.id[i]]
      
      model <- result.AB
      
      return(as.integer(model))
    }, mc.cores = 30, ignore.interactive = F)
    
    result.cit.ecit.deseq$groundtruth <- do.call(c, result.lapply)
    result.cit.ecit.deseq$groundtruth <- as.character(result.cit.ecit.deseq$groundtruth)
    
    PR.run.i <- precision.recall.fdr(result.cit.ecit.deseq)
    PR.runs[[1 + length(PR.runs)]] <- PR.run.i
    
  }
  
  PR.runs <- do.call(rbind, PR.runs)
  
  PR.runs.call.melt <- melt(PR.runs %>% select(cutoff, CIT = cit.calls, ECIT = ecit.calls), 
                            id.vars = 'cutoff',
                            variable.name = 'method',
                            value.name = 'calls') 
  
  PR.runs.call.melt.summary <- PR.runs.call.melt %>% group_by(cutoff, method) %>% 
    summarise(calls.mu = mean(calls), calls.sd = sd(calls))
  
  PR.runs.precision.melt <- melt(PR.runs %>% select(cutoff, CIT = cit.prec, ECIT = ecit.prec),
                                 id.vars = 'cutoff',
                                 variable.name = 'method',
                                 value.name = 'precision')
  
  
  PR.runs.precision.melt.summary <- PR.runs.precision.melt %>% group_by(cutoff, method) %>%
    summarise(precision.mu = mean(precision), precision.sd = sd(precision))
  
  
  PR.runs.TP.melt <- melt(PR.runs %>% select(cutoff, CIT = cit.truecalls, ECIT = ecit.truecalls),
                          id.vars = 'cutoff',
                          variable.name = 'method',
                          value.name = 'TP')
  
  PR.runs.TP.melt.summary <- PR.runs.TP.melt %>% group_by(cutoff, method) %>%
    summarise(TP.mu = mean(TP), TP.sd = sd(TP))
  
  
  library(ggplot2)
  library(reshape2)
  
  
  
  pcalls <- PR.runs.call.melt.summary %>% ggplot(aes(x = as.factor(cutoff), y = calls.mu, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = calls.mu - calls.sd, ymax = calls.mu + calls.sd), width = .2, position = position_dodge(1)) +
    labs(title="Performance on \n yeast dataset (all genes)", x="FDR Cutoff", y="Avg # Causal Calls", fill="")
  
  pcalls <- pcalls + 
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  
  
  pcalls
  
  ggsave(plot = pcalls,
         units = "in",
         height = 4.5,
         width = 5.7,
         filename = 'output/yeast.causalcall.all.pdf')
  
  pprecision <- PR.runs.precision.melt.summary %>% ggplot(aes(x = as.factor(cutoff), y = precision.mu, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = precision.mu - precision.sd, ymax = precision.mu + precision.sd), width = .2, position = position_dodge(1)) +
    labs(title="Performance on \n yeast dataset (all genes)", x="FDR Cutoff", y="Avg Precision", fill="") 
  
  pprecision <- pprecision +
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  
  pprecision
  ggsave(plot = pprecision,
         units = "in",
         height = 4.5,
         width = 5.7,
         filename = 'output/yeast.precision.all.pdf')
  
  
  pTP <- PR.runs.TP.melt.summary %>% ggplot(aes(x = as.factor(cutoff), y = TP.mu, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = TP.mu - TP.sd, ymax = TP.mu + TP.sd), width = .2, position = position_dodge(1)) +
    labs(title="Performance on \n yeast dataset (all genes)", x="FDR Cutoff", y="Avg True Positive", fill="")
  
  pTP <- pTP + 
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  
  
  pTP
  
  ggsave(plot = pTP,
         units = "in",
         height = 4.5,
         width = 5.7,
         filename = 'output/yeast.TP.all.pdf')
  
  


}
precision.recall.all.fdr()
ground.truth.generation.yeastract.plus <- function() {
  # this chunk is from yeast1k_analysis.R
  # Directly downlaod from yeastract+. 
  rm(list = ls())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  path.ground.truth <- 'data/yeast/RegulationMatrix_Documented_2021519_720_541750229.80cross3394.Dnabinding.plus.expression.csv'
  ground.truth <- read.csv2(path.ground.truth, header = T, stringsAsFactors = F, strip.white = T, check.names = F)
  colnames(ground.truth) <- c('TF.TG', colnames(ground.truth[, -c(1)]))
  
  ensembl <- read.table(file = 'data/FindrCausalNetworkInferenceOnYeast/data/output/cleaned_genes_ensembl83_yeast_step3.txt')
  ensembl <- ensembl[, c(1, 4, 10, 11)]
  colnames(ensembl) <- c('chr', 'start', 'gene.id', 'gene.name')
  ensembl$gene.name <- sub('Name=', '', ensembl$gene.name)  
  head(ensembl)
  
  
  grepl('p$', ground.truth$TF.TG) %>% sum()
  ground.truth$TF.TG[grepl('p$', ground.truth$TF.TG) == F ]
  ground.truth$TF.TG <- sub(pattern = 'p$', '', ground.truth$TF.TG)
  ground.truth$TF.TG <- stringr::str_to_upper(ground.truth$TF.TG)
  
  colnames(ground.truth)[1:10]
  colnames(ground.truth) <- stringr::str_to_upper(colnames(ground.truth))
  
  sum(ground.truth$TF.TG %in% ensembl$gene.name)
  indices <- match(ground.truth$TF.TG, ensembl$gene.name)
  ground.truth$TF.TG[!is.na(indices)] <- as.character(ensembl$gene.id[ indices[!is.na(indices)] ])
  
  sum(colnames(ground.truth) %in% ensembl$gene.name)
  
  colnames(ground.truth)[(colnames(ground.truth) %in% ensembl$gene.name) == F ]
  indices <- match(colnames(ground.truth), ensembl$gene.name)
  colnames(ground.truth)[!is.na(indices)] <- as.character(ensembl$gene.id[indices[!is.na(indices)]])
  
  sum(colnames(ground.truth) %in% ensembl$gene.id)
  
  mrna <- readRDS('data/findr.exp.log2.deseq.residual.rds')
  sum(ground.truth$TF.TG %in% colnames(mrna))  
  
  ground.truth$TF.TG[(ground.truth$TF.TG %in% colnames(mrna)) == F]
  ground.truth <- ground.truth[ground.truth$TF.TG %in% colnames(mrna), ]
  
  
  present <- colnames(ground.truth) %in% colnames(mrna)
  table(present)
  colnames(ground.truth)[present == F]
  
  biomart <- read.table('data/mart_export.txt', header = T, check.names = F, stringsAsFactors = F, strip.white = T, fill = T)
  present <- (colnames(ground.truth) %in% biomart$Gene.name)
  colnames(ground.truth)[present]
  
  indices <- match(colnames(ground.truth), biomart$Gene.name )
  colnames(ground.truth)[!is.na(indices)] <- as.character(biomart$Gene.stable.ID[indices[!is.na(indices)]])
  present <- colnames(ground.truth) %in% colnames(mrna)
  table(present)
  
  present[1] <- T
  ground.truth <- ground.truth[, present]
  
  sum(rowSums(ground.truth[, -c(1)]))
  
  saveRDS(ground.truth, 'data/yeast/RegulationMatrix_Documented_2021519_720_541750229.80cross3394.Dnabinding.plus.expression.csv.rds')
  
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
DESeq_alternative_using_package <- function() {
  rm(list = ls())
  # https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
  

  mrna.count <- readRDS('data/yeast/mrna.count.rds')
  mrna.count[, 1:ncol(mrna.count)] <- sapply(mrna.count[, 1:ncol(mrna.count)], as.integer)
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = mrna.count, colData = as.data.frame(colnames(mrna.count)), design = ~ 1)
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds) %>% View
  mrna.count.normalized <- counts(dds, normalized=TRUE)
  
  saveRDS(sizeFactors(dds), file = 'data/yeast/mrna.count.DESeq.sizeFactors.rds')
  saveRDS(mrna.count.normalized, file = "data/yeast/mrna.count.DESeqnormalized.rds")
  
  mrna.count.normalized.mu <- apply(mrna.count.normalized[, ], 1, mean)
  mrna.count.normalized.var <- apply(mrna.count.normalized[, ], 1, var)
  
  ggplot() + geom_point(aes(x = mrna.count.normalized.mu, y = mrna.count.normalized.var)) + scale_x_log10() + scale_y_log10() + geom_abline(intercept = 0, slope = 1, color ='red')
  
  # dummy gene simulation
  s_hat <- psych::geometric.mean(dds$sizeFactor)
  dummy.rows <- 2000  
  dummy.cols <- 1012
  dummy.mrna.count <- matrix(nrow = dummy.rows, ncol = dummy.cols)  
  rownames(dummy.mrna.count) <- paste0("dummy", seq(1, dummy.rows))
  
  for (i in seq(1:nrow(dummy.mrna.count))) {
    lambda <- s_hat * sample(mrna.count.normalized.mu, size = 1)
    dummy.mrna.count[i, ] <- rpois(n = ncol(dummy.mrna.count), lambda = lambda)
  }
  
  dummy.mrna.count.mu <- apply(dummy.mrna.count, 1, mean)
  dummy.mrna.count.var <- apply(dummy.mrna.count, 1, var)
  
  mrna.count.normalized.and.dummy <- rbind(mrna.count.normalized, dummy.mrna.count)
  mrna.count.normalized.and.dummy.mu <- apply(mrna.count.normalized.and.dummy, 1, mean)
  mrna.count.normalized.and.dummy.var <- apply(mrna.count.normalized.and.dummy, 1, var)
  plot.data <- data.frame(mu = mrna.count.normalized.and.dummy.mu, var = mrna.count.normalized.and.dummy.var)
  plot.data$is_dummy <- grepl("dummy*", rownames(mrna.count.normalized.and.dummy))  
  
  
  result.cit.ecti <- readRDS('data/yeast/600.3/dna.and.expression/3.findr.result.cit.ecit.runs.i1.rds')
  cutoff <- 0.4
  noisy.genes.ids <- result.cit.ecti$A.id[result.cit.ecti$A.id.me.var.to.totalvar >= cutoff]
  noisy.genes.ids2 <- result.cit.ecti$B.id[result.cit.ecti$B.id.me.var.to.totalvar >= cutoff]
  noisy.genes.ids <- c(noisy.genes.ids, noisy.genes.ids2)
  noisy.genes.ids <- unique(noisy.genes.ids)
    
  noisy.genes.ids <- sub(pattern = '-', replacement = '.',  noisy.genes.ids)
    
  plot.data.noisy <- plot.data[noisy.genes.ids, ]
    
  temp <- 1
  
  
  pcount <- ggplot(data = plot.data , aes(x = log2(mu), y = log2(var))) + 
    geom_point(aes(color = as.factor(is_dummy)), alpha=1)   + 
    geom_abline(intercept = 0, slope = 1, color = 'red') + 
    ggtitle('DESeq Normarlized counts(Yeast) \n black highlights noisy genes') +
    labs(x = 'log2 mean', y = 'log2 variance', colour = 'Poisson Sampled' ) +
    geom_point(data = plot.data.noisy, aes(x = log2(mu), y = log2(var)), color='black', shape=1, alpha=0.7)
  
  pcount
  
  pcount <- pcount + 
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  pcount
  
  plot.save.manuscript(pcount, path='output/', filename = 'yeast.deseqnormalised.count')
  
  mrna.count.normalized.and.dummy.log2.transformed <- log2(mrna.count.normalized.and.dummy + 0.5)
  plot.data.normalized <- data.frame(mu = apply(mrna.count.normalized.and.dummy.log2.transformed, 1, mean), var = apply(mrna.count.normalized.and.dummy.log2.transformed, 1, var))        
  
  plot.data.normalized$is_dummy <- grepl("dummy*", rownames(mrna.count.normalized.and.dummy.log2.transformed))
  
  plot.data.normalized.noisy <- plot.data.normalized[noisy.genes.ids, ]
  
  ptranformed <- ggplot(data = plot.data.normalized, aes(x = mu,y = var)) + 
    geom_point(aes(color = as.factor(is_dummy)), alpha=1) + 
    ggtitle('log2 transformed DESeq \n Normarlized counts Yeast \n black highlights noisy genes') +
    labs(x = 'mean', y = 'variance', colour = 'Poisson Sampled') +
    geom_point(data = plot.data.normalized.noisy, aes(x = mu, y = var), shape=1, alpha=0.7)
  
  ptranformed
  
  ptranformed <- ptranformed + 
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  ptranformed
  plot.save.manuscript(ptranformed, 'output/','yeast.log2.tranformed.Deseq')
  
  temp <- 1  
}



measurement_error_estimation <- function() {
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
  
  ggplot() + geom_point(aes(x = me.var.log2$total.var.log2, y = me.var.log2$me.var.log2.estimated))
}







