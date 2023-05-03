setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('data/gtex/')


#'gtex gct have genenXsample.tissue count data.
#'get genesXsamples for specific tissue.
#'filter only protein coding genes

meqtl.cit.ecit()

meqtl.cit.ecit <- function() {
  rm(list = ls())
  library(tidyr)
  library(cit)
  
  library(parallel)
  
  
  
  path.tis2.tpm <- 'GTEx_Analysis_v7_eQTL_expression_matrices/Muscle_Skeletal.v7.normalized_expression.bed.gz'
  path.genotype <- 'L.Asub.Testis.NerTib.Thy.Pan.Msk.GT.FORMAT'
  path.cov <- 'GTEx_Analysis_v7_eQTL_covariates/Muscle_Skeletal.v7.covariates.txt'
  path.meqtl.expr <- 'Muscle-Skeletal/Muscle Skeletal/matrixeqtl/expression.txt'
  path.meqtl.geno <- 'Muscle-Skeletal/Muscle Skeletal/matrixeqtl/genotype.txt'
  path.meqtl.cov <- 'Muscle-Skeletal/Muscle Skeletal/matrixeqtl/covariates.txt'
  path.meqtl.genloc <- 'Muscle-Skeletal/Muscle Skeletal/matrixeqtl/geneloc.txt'
  path.meqtl.snploc <- 'Muscle-Skeletal/Muscle Skeletal/matrixeqtl/snpsloc.txt'
  path.gencode <- 'gencode.v19.gene.id.names.map.txt'
  path.egenes <- 'GTEx_Analysis_v7_eQTL/Muscle_Skeletal.v7.egenes.txt.gz'
  path.cis.egenes <- 'Muscle-Skeletal/Muscle Skeletal/qval_cis_egenes.txt'
  
  path.tis2.tpm
  expression <- read.csv2(gzfile(path.tis2.tpm), sep = '\t', 
                          header = T, check.names = F, strip.white = T, stringsAsFactors = F)
  expression[1:4, 1:4]
  
  expression$gene_id <- tstrsplit(expression$gene_id, '[.]')[[1]]
  
  
  gencode <- read.csv2(path.gencode, header = T, sep = '\t', check.names = F,
                       strip.white = T, stringsAsFactors = F)
  gencode$gene_id <- tstrsplit(gencode$gene_id, '[.]')[[1]]
  rownames(gencode) <- gencode$gene_id
  
  cis.egenes <- read.table(gzfile(path.egenes), 
                           header = T, check.names = F, fill = T)
  
  cis.egenes <- cis.egenes %>% group_by(gene_id) %>% dplyr::summarise(variant_id=variant_id[which.min(qval)], 
                                                                      qval.cis=qval[which.min(qval)]) 
  
  cis.egenes <- cis.egenes %>% filter(qval.cis < 0.05)
  
  path.cis.egenes
  write.table(cis.egenes, path.cis.egenes, row.names = F, quote = F)
  
  gencode <- gencode[expression$gene_id, ]
  gencode$chr <- tstrsplit(gencode$Chromosome, ':')[[1]]
  gencode$start <- tstrsplit(tstrsplit(gencode$Chromosome, ':')[[2]] , '-')[[1]]
  gencode$end <- tstrsplit(tstrsplit(gencode$Chromosome, ':')[[2]] , '-')[[2]]
  
  geneloc <- gencode %>% select(geneid= gene_id, chr=chr, s1=start, s2=end)
  rownames(geneloc) <- NULL
  
  expression <- expression %>% select(-c(`#chr` , start, end))
  
  path.genotype
  genotype <- read.table(path.genotype, header = T, check.names = F, strip.white = T,
                         stringsAsFactors = F)
  genotype[1:4, 1:4]
  
  
  # keep genotype only for significant cis-snp
  cis.snps <- unique(cis.egenes$variant_id)
  cis.snps[1:4]
  str.tokens <- tstrsplit(cis.snps, '_')
  cis.snps <- paste(str.tokens[[1]], str.tokens[[2]], sep = '_')
  
  common.snps <- intersect(cis.snps, paste(genotype$CHROM, genotype$POS, sep = '_'))
  
  genotype <- genotype[match(common.snps, paste(genotype$CHROM, genotype$POS, sep = '_')), ]
  
  snpsloc <- cbind(snp = paste(genotype$CHROM, genotype$POS, sep = '_'), genotype %>% select(CHROM, POS))
  
  genotype <- cbind(id = paste(genotype$CHROM, genotype$POS, sep = '_'), genotype %>% select(-c(CHROM, POS)))
  
  snpsloc[1:4, ]
  
  genotype[1:4, 1:4]
  
  
  path.cov
  covariates <- read.table(path.cov, header = T,
                           strip.white = T, stringsAsFactors = F, check.names = F)
  # keep rows and columns in same align
  indv.exp <- colnames(expression)[-c(1)]
  indv.geno <- colnames(genotype)[-c(1)]
  indv.cov <- colnames(covariates)[-c(1)]
  indv.comm <- intersect(indv.exp, intersect(indv.geno, indv.cov))
  length(indv.comm)
  
  expression <- expression %>% select(gene_id, all_of(indv.comm))
  genotype <- genotype %>% select(id, all_of(indv.comm))
  covariates <- covariates %>% select(ID, all_of(indv.comm))
  
  
  # all colunms of genotype, expression, covariates same aligned
  all((colnames(genotype) == colnames(expression))[-c(1)])
  all((colnames(genotype) == colnames(covariates))[-c(1)])
  
  all(geneloc$gene_id == expression$gene_id)
  all(snpsloc$snp == genotype$id)
  path.gencode <- 'gencode.v19.gene.id.names.map.txt'
  gencode <- read.csv2(path.gencode, header = T, sep = '\t', check.names = F,
                       strip.white = T, stringsAsFactors = F)
  gencode$gene_id <- tstrsplit(gencode$gene_id, '[.]')[[1]]
  rownames(gencode) <- gencode$gene_id
  
  
  gencode$chr <- tstrsplit(gencode$Chromosome, ':')[[1]]
  gencode$start <- tstrsplit(tstrsplit(gencode$Chromosome, ':')[[2]] , '-')[[1]]
  gencode$end <- tstrsplit(tstrsplit(gencode$Chromosome, ':')[[2]] , '-')[[2]]
  

  
  library(dplyr)
  
  path.tis1.log2.expr.resi <- 'Muscle-Skeletal/tissue.mrna.log2.DESeqnormalized.rediual.rds'
  path.tis2.log2.expr.resi <- path.tis1.log2.expr.resi
  path.tis1.geno <- 'Muscle-Skeletal/Muscle Skeletal/matrixeqtl/genotype.txt'
  path.trios <- 'Muscle-Skeletal/Muscle Skeletal/trios.rds'
  path.tiss1.me <- 'Muscle-Skeletal/me.var.log2.rds'
  path.tiss2.me <- path.tiss1.me
  path.cit.ecit <- 'Muscle-Skeletal/Muscle Skeletal/trios.cit.ecit.result.rds'
  
  tis1.log2.expr.resi <- readRDS(path.tis1.log2.expr.resi)
  tis2.log2.expr.resi <- readRDS(path.tis2.log2.expr.resi)
  
  path.tis1.geno
  
  genotype <- read.table(path.tis1.geno, header = T, strip.white = T,
                         stringsAsFactors = F, check.names = F)
  rownames(genotype) <- genotype$id
  genotype <- genotype[, -c(1)]
  
  table(genotype == './.')
  genotype[(genotype == './.')] <- NA
  table(is.na(genotype))
  


  all.trios <- readRDS(path.cit.ecit)
  
  trios <- all.trios
  
  
    biomart <- read.table('biomart.GRCh37.p13.txt', header = F, strip.white = T, stringsAsFactors = F)
    biomart <- biomart %>% select(gene.id = V1, gene.symbol = V3, gene.type = V4)
    head(biomart)
    
    trios <- merge(trios, biomart, by.x = 'gene_id.x', by.y = 'gene.id', all.x = T)
    trios <- merge(trios, biomart, by.x = 'gene_id.y', by.y = 'gene.id', all.x = T)
    
    
    temp <- 1
  
  
  # go through this code to add mappability columns
    gene.mappbility <- read.table(gzfile('hg19_gencode19_75merExon_36merUTR_2mismatch_gene_mappability.txt.gz'),
                                  header = F, strip.white = T, stringsAsFactors = F)
    colnames(gene.mappbility) <- c('gene.id', 'score')
    gene.mappbility$gene.id <- tstrsplit(gene.mappbility$gene.id, '[.]')[[1]]
    head(gene.mappbility)
    
    match.indices <- match(trios$gene_id.x, gene.mappbility$gene.id)
    table(is.na(match.indices))
    
    trios$gene_id.x.mappability[!is.na(match.indices)] <- gene.mappbility$score[match.indices[!is.na(match.indices)]]
    
    match.indices <- match(trios$gene_id.y, gene.mappbility$gene.id)
    table(is.na(match.indices))
    
    trios$gene_id.y.mappability[!is.na(match.indices)] <- gene.mappbility$score[match.indices[!is.na(match.indices)]]
    
    mappability.thres <- 0.70
    table(trios$gene_id.x.mappability > mappability.thres & trios$gene_id.y.mappability > mappability.thres)
    
    trios <- trios[trios$gene_id.x.mappability > mappability.thres & trios$gene_id.y.mappability > mappability.thres, ]
    
  
  
  trios$cit.AB.p_cit.BH <- p.adjust(trios$cit.AB.p_cit, method = 'BH')
  trios$ecit.AB.adj_p_cit.BH <- p.adjust(trios$ecit.AB.adj_p_cit, method = 'BH')
  trios$cit.BA.p_cit.BH <- p.adjust(trios$cit.BA.p_cit, method = 'BH')
  trios$ecit.BA.adj_p_cit.BH <- p.adjust(trios$ecit.BA.adj_p_cit, method = 'BH')
  
  
  
  thres <- 0.05
  trios$cit.direction <- get_cit_direction(trios$cit.AB.p_cit, trios$cit.BA.p_cit, thres)
  
  trios$ecit.direction <- get_cit_direction(trios$ecit.AB.adj_p_cit, trios$ecit.BA.adj_p_cit, thres)
  
  table(trios$cit.direction)
  table(trios$ecit.direction)
  
  library(caret)
  cfmtx <- confusionMatrix(data = as.factor(trios$ecit.direction), reference = as.factor(trios$cit.direction))
  cfmtx$table
  
  View(trios)
  
  
  trios %>% select(variant_id, gene.symbol.x, gene.symbol.y, 
                   gene_id.x.mappability, gene_id.y.mappability,
                   #gene.type.x, gene.type.y,
                   ecit.AB.adj_p_cit, ecit.BA.adj_p_cit,
                   ecit.AB.adj_p_cit.BH, ecit.BA.adj_p_cit.BH,
                   ecit.direction) %>% View()
  
  cytoscape <- function() {
    
    causal <- trios %>% select(source.id = gene_id.x, target.id = gene_id.y, 
                               source.sym = gene.symbol.x, target.sym = gene.symbol.y, 
                               source.map = gene_id.x.mappability, target.map = gene_id.y.mappability, 
                               rho.AB,
                               ecit.direction) %>% filter(ecit.direction == 1)
    
    reactive <- trios %>% select(source.id = gene_id.y, target.id = gene_id.x, 
                                 source.sym = gene.symbol.y, target.sym = gene.symbol.x, 
                                 source.map = gene_id.y.mappability, target.map = gene_id.x.mappability, 
                                 rho.AB,
                                 ecit.direction) %>% filter(ecit.direction == 2)
    
    
    cytoscape.result <- rbind(causal, reactive)
    cytoscape.result <- cytoscape.result %>% select(-c(ecit.direction))
    cytoscape.result$edgetype <- 'up'
    cytoscape.result$edgetype[cytoscape.result$rho.AB < 0] <- 'down'
    
    cytoscape.result$source.sym[is.na(cytoscape.result$source.sym)] <- cytoscape.result$source.id[is.na(cytoscape.result$source.sym)]
    cytoscape.result$target.sym[is.na(cytoscape.result$target.sym)] <- cytoscape.result$target.id[is.na(cytoscape.result$target.sym)]
    
    table(cytoscape.result$ecit.direction)
    write.table(cytoscape.result, 'Muscle-Skeletal/Muscle Skeletal/cytoscape.muscle.txt', quote = F, row.names = F, sep = '\t')
    temp <- 1  
  }
  # causal example
  L.id <- '13_74110412'
  A.id <- 'ENSG00000102554'
  B.id <- 'ENSG00000198324'
  filename<-"cid.causal"
  boxplot.analysis()
  # independent example 16_81080366
  L.id <- '11_9242168'
  A.id <- 'ENSG00000184014'
  B.id <- 'ENSG00000106089'
  filename<-"cid.indep"
  boxplot.analysis()
  
  boxplot.analysis <- function() {

    A.symbol <- biomart$gene.symbol[match(A.id, biomart$gene.id)]
    B.symbol <- biomart$gene.symbol[match(B.id, biomart$gene.id)]
    
    L <- as.integer(genotype[L.id, ])
    
    indv.not.na <- !is.na(L)
    
    L <- L[indv.not.na]
    
    
    L1 <- rep(1, length(L))
    L2 <- rep(1, length(L))  
    
    L1[L == 0] <- 0
    L2[L == 0] <- 0
    
    L1[L == 1] <- 0
    L2[L == 1] <- 1
    
    A <- as.numeric(tis1.log2.expr.resi[A.id, indv.not.na])
    A <- scale(A, scale = F)
    B <- as.numeric(tis2.log2.expr.resi[B.id, indv.not.na])
    B <- scale(B, scale = F)
    
    table(L)
    
    table(L1, L2)
    
    me.var.A <- tis1.me.var$me.var.log2.estimated[match(A.id, tis1.me.var$gene.id)]
    me.var.B <- tis2.me.var$me.var.log2.estimated[match(B.id, tis2.me.var$gene.id)]
    
    cit.cp(cbind(L1, L2), A, B, n.resampl = 100)
    cit.cp(cbind(L1, L2), B, A, n.resampl = 100)
    
    
    AB <- qplot(A, B) + geom_smooth() + xlab(A.symbol) + ylab(B.symbol) +
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
    
    AB
    LB <- ggplot()  +geom_point(aes(as.factor(L), B))+geom_boxplot(aes(as.factor(L), B), notch = T) + xlab(L.id) + ylab(B.symbol) +
      
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")# L ~ B
    
    
    LB
    LA_B <- ggplot() + geom_point(aes(as.factor(L), A))+geom_boxplot(aes(as.factor(L), A), notch = T) + xlab(L.id) + 
      ylab(paste(A.symbol, '|', B.symbol, sep = ' ')) +  
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")#L ~ A|B
    
    AB_L <- qplot(A, residuals(lm(B~as.factor(L)))) + geom_smooth() +
      xlab(A.symbol) + ylab(paste(B.symbol, '|', L.id, sep = ' ')) +
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")# A ~ B|L
    
    
    L.ind.B_A <- ggplot() + geom_point(aes(as.factor(L), residuals(lm(B~A))))+geom_point(mapping = aes(x=as.factor(L), y=B))+geom_boxplot(aes(as.factor(L), residuals(lm(B~A))), notch = T) + 
      
      xlab(L.id) + ylab(paste(B.symbol, '|', A.symbol, sep = ' ')) +
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")# L indp B|A
    
    library(cowplot)
    
    ymin <- -2
    ymax <- 2
    plot.cid <- plot_grid(
      LA_B + ylab('Relative Expression') + ggtitle(A.symbol) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylim(ymin, ymax), 
      LB + ggtitle(B.symbol) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) +
        ylim(ymin, ymax),
      
      L.ind.B_A + ggtitle(paste(B.symbol, '|', A.symbol, sep = '')) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        ylim(ymin, ymax),
      align = 'hv',
      nrow = 1, ncol = 3,
      labels = 'auto'
    )
    
    print(plot.cid)
    
    
    save.manuscript <- function(plot, filename) {
      for (exten in c('.pdf', '.jpeg', '.svg')) {
        ggsave(plot = plot,
               units = "in",
               height = 4.5,
               width = 5.7,
               filename = paste("../../output/", filename, exten, sep = ''))
      }
    }
    
    save.manuscript(plot.cid, filename='cid.indep')
    
    
    
    temp <- 1
    
  }
}
