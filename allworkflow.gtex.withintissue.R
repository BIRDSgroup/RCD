setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('data/gtex/')

common.samples.list <- function() {
  #'extract samples and genes available for both tissues from recount
  #'https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
  
  rm(list = ls())
  loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  
  library(SummarizedExperiment)
  tissue1.rse <- loadRData('input/rse_gene_muscle.Rdata')  
  tissue1.counts <- data.frame(assays(tissue1.rse)$counts)
  
  tissue2.rse <- loadRData('input/rse_gene_muscle.Rdata')  
  tissue2.counts <- data.frame(assays(tissue1.rse)$counts)
  
  
  
  
  temp <- 1  
}

#'gtex gct have genenXsample.tissue count data.
#'get genesXsamples for specific tissue.
#'filter only protein coding genes

deseq.count.deseq.resi.deseq.me.var.all <- function() {
  rm(list = ls())
    # gtex.gct <- read.table('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct', skip = 2,
  #                        header = T, strip.white = T, stringsAsFactors = F, check.names = F)
  # 
  
  # saveRDS(gtex.gct, 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.rds')
  
  path.gtex.gct <- 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.rds'
  path.tissue.mrna.log2.DESeqnormalized <- 'Muscle-Skeletal/tissue.mrna.log2.DESeqnormalized.rds'
  path.me.var.log2.rds<- 'Muscle-Skeletal/me.var.log2.rds'
  path.cov <- 'GTEx_Analysis_v7_eQTL_covariates/Muscle_Skeletal.v7.covariates.txt'
  path.tissue.mrna.log2.DESeqnormalized.rediual <- 'Muscle-Skeletal/tissue.mrna.log2.DESeqnormalized.all.cov.rediual.rds'
  
  gtex.gct <- readRDS(path.gtex.gct)
  
  gtex.gct$Name <- tstrsplit(gtex.gct$Name, '[.]')[[1]]
  
  
  gtex.attributes <- read.csv2('input/GTEx_v7_Annotations_SampleAttributesDS.txt', header = T,
                               strip.white = T, stringsAsFactors = F, sep = '\t')
  gtex.attributes <- as.data.frame(gtex.attributes)
  gtex.attributes %>% select(SAMPID, SMTS, SMTSD ) %>% View()
  
  tissue.name <- 'Muscle - Skeletal'  #'Thyroid' #'Adipose - Subcutaneous' #'Muscle - Skeletal'  #'Pancreas'
  samples.ids.tissue <- gtex.attributes$SAMPID[grepl(pattern = tissue.name, gtex.attributes$SMTSD)]
  
  length(samples.ids.tissue)
  
  colnames(gtex.gct) <- gsub(pattern = '[.]', '-', colnames(gtex.gct))
  
  
  
  table(samples.ids.tissue %in% colnames(gtex.gct))
  
  samples.ids.tissue <- samples.ids.tissue[samples.ids.tissue %in% colnames(gtex.gct)]
  
  tissue.mrna.count <- gtex.gct[, samples.ids.tissue]
  rownames(tissue.mrna.count) <- gtex.gct$Name
  
  tissue.mrna.count[1:4, 1:4]
  
  # keep <- apply(tissue.mrna.count, 1, function(x) {
  #   a <- sum(x >= 6)
  #   return(a >= 0.1*length(x))
  # })
  # 
  # table(keep)
  # 
  # tissue.mrna.count <- tissue.mrna.count[keep, ]
  
  
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = tissue.mrna.count, colData = as.data.frame(colnames(tissue.mrna.count)), design = ~ 1)
  dds <- estimateSizeFactors(dds)
  tissue.mrna.count.DESeqnormalized <- counts(dds, normalized=TRUE)
  
  
  dim(tissue.mrna.count.DESeqnormalized)
  
  library(pbmcapply)
  
  sizefactors <- sizeFactors(dds)
  head(sizefactors)
  
  s_hat <- psych::geometric.mean(sizefactors)
  
  row <- 2
  me.var.log2 <- pbmclapply(X = 1:nrow(tissue.mrna.count.DESeqnormalized[ , ]), FUN = function(row) {
    dummy.bootstraps <- 500  
    dummy.samples <- ncol(tissue.mrna.count.DESeqnormalized)
    dummy.mu <- mean(as.numeric(tissue.mrna.count.DESeqnormalized[row, ]))
    
    me.bootstraps.me.var.log2 <- sapply(1:dummy.bootstraps, FUN = function(i) {
      dummy.count <- rpois(n = dummy.samples, lambda = s_hat*dummy.mu)
      dummy.count.log2 <- log2(dummy.count + 0.5)
      dummy.var.log2 <- var(dummy.count.log2)  
      return(dummy.var.log2)
    })
    
    ggplot() + geom_histogram(aes(me.bootstraps.me.var.log2))
    
    me.bootstraps.me.var.log2.mu <- mean(me.bootstraps.me.var.log2)
    result <- list('gene.id' = rownames(tissue.mrna.count.DESeqnormalized)[row] , 'total.var.log2' = var(log2(as.numeric(tissue.mrna.count.DESeqnormalized[row, ] + 0.5))), 
                   'me.var.log2.estimated' = me.bootstraps.me.var.log2.mu, 'mu.log2' = mean(log2(as.numeric(tissue.mrna.count.DESeqnormalized[row, ] + 0.5))))
    
    cat(unlist(result), '\n')
    
    return(result)
  }, mc.cores = 50, ignore.interactive = F)  
  
  me.var.log2 <- do.call(rbind, me.var.log2)  
  me.var.log2 <- as.data.frame(me.var.log2)
  me.var.log2$gene.id <- as.character(me.var.log2$gene.id) 
  
  # here total.var.log2 is before cov adjustement.
  me.var.log2$total.var.log2 <- as.numeric(me.var.log2$total.var.log2)
  me.var.log2$me.var.log2.estimated <- as.numeric(me.var.log2$me.var.log2.estimated)
  me.var.log2$mu.log2 <- as.numeric(me.var.log2$mu.log2)
  # me.var.log2$me.to.total.var <- me.var.log2$me.var.log2.estimated / me.var.log2$total.var.log2 * 100
  
  ggplot() + geom_point(aes(x = me.var.log2$total.var.log2, y = me.var.log2$me.var.log2.estimated)) +
    scale_x_log10() + scale_y_log10()
  
  ggplot() + geom_point(aes(x = me.var.log2$mu.log2, y = me.var.log2$me.var.log2.estimated)) +
    ggtitle(tissue.name)
  
  library(reshape2)
  data.melt <- melt(select(me.var.log2, mu.log2, total.var.log2, me.var.log2.estimated), 
                    id = c('mu.log2'), 
                    variable.name = 'variance',
                    value.name = 'log2.var')
  
  data.melt %>% ggplot(aes(x = mu.log2, y = log2.var, color = as.factor(variance))) +
    geom_point() + ggtitle(tissue.name)
  
  
  
  table(me.var.log2$me.var.log2.estimated < me.var.log2$total.var.log2)
  
  
  ggplot() + geom_point(aes(log2(tissue.mrna.count.DESeqnormalized[,1]), log2(tissue.mrna.count.DESeqnormalized[,2]))) +
    ggtitle(tissue.name)
  
  # path.output.dir <- paste0(tissue.name)
  # print(path.output.dir)
  # dir.create(path.output.dir, showWarnings = T)
  # 
  tissue.mrna.log2.DESeqnormalized <- log2(tissue.mrna.count.DESeqnormalized + 0.5)
  colnames(tissue.mrna.log2.DESeqnormalized) <- paste(tstrsplit(colnames(tissue.mrna.log2.DESeqnormalized), '-')[[1]], 
                                                      tstrsplit(colnames(tissue.mrna.log2.DESeqnormalized), '-')[[2]], sep = '-')
  tissue.mrna.log2.DESeqnormalized[1:3, 1:4]
  
  
  
  path.tissue.mrna.log2.DESeqnormalized
  saveRDS(tissue.mrna.log2.DESeqnormalized, path.tissue.mrna.log2.DESeqnormalized)
  
  
  path.me.var.log2.rds
  saveRDS(me.var.log2, path.me.var.log2.rds)
  
  # adjust for covariated
  path.tissue.mrna.log2.DESeqnormalized
  tissue.mrna.log2.DESeqnormalized <- readRDS(path.tissue.mrna.log2.DESeqnormalized)
  
  path.cov
  cov <- read.table(path.cov, 
                    header = T, strip.white = T, stringsAsFactors = F, check.names = F)
  
  rownames(cov) <- cov$ID
  cov <- cov[ , -c(1)]
  
  cov[1:4, 1:4]
  
  samples.com <- intersect(colnames(tissue.mrna.log2.DESeqnormalized), 
                           colnames(cov))
  length(samples.com)
  
  tissue.mrna.log2.DESeqnormalized <- tissue.mrna.log2.DESeqnormalized[, samples.com]
  cov <- cov[, samples.com]
  
  library(pbmcapply)
  
  row.i <- 1
  
  tissue.mrna.log2.DESeqnormalized.rediual <- pbmclapply(1:nrow(tissue.mrna.log2.DESeqnormalized), FUN = function(row.i) {
    A <- (tissue.mrna.log2.DESeqnormalized[row.i, ])
    
    model <- lm(A ~ t(cov))
    
    A.residual <- residuals(model) + model$coefficients['(Intercept)']
    
    return(A.residual)    
  }, mc.cores = 50, ignore.interactive = F)
  
  tissue.mrna.log2.DESeqnormalized.rediual <- as.data.frame(do.call(rbind, tissue.mrna.log2.DESeqnormalized.rediual))
  rownames(tissue.mrna.log2.DESeqnormalized.rediual) <- rownames(tissue.mrna.log2.DESeqnormalized)
  
  
  print(path.tissue.mrna.log2.DESeqnormalized.rediual)
  saveRDS(tissue.mrna.log2.DESeqnormalized.rediual, path.tissue.mrna.log2.DESeqnormalized.rediual)
  
  temp <- 1
}

deseq.mean.var.plots.generate <- function() {
  rm(list = ls())
  
  
  # gtex.gct <- read.table('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct', skip = 2,
  #                        header = T, strip.white = T, stringsAsFactors = F, check.names = F)
  # 
  
  # saveRDS(gtex.gct, 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.rds')
  
  path.gtex.gct <- 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.rds'
  path.tissue.mrna.log2.DESeqnormalized <- 'Muscle-Skeletal/tissue.mrna.log2.DESeqnormalized.rds'
  path.me.var.log2.rds<- 'Muscle-Skeletal/me.var.log2.rds'
  path.cov <- 'GTEx_Analysis_v7_eQTL_covariates/Muscle_Skeletal.v7.covariates.txt'
  path.tissue.mrna.log2.DESeqnormalized.rediual <- 'Muscle-Skeletal/tissue.mrna.log2.DESeqnormalized.all.cov.rediual.rds'
  
  gtex.gct <- readRDS(path.gtex.gct)
  
  gtex.gct$Name <- tstrsplit(gtex.gct$Name, '[.]')[[1]]
  
  
  gtex.attributes <- read.csv2('input/GTEx_v7_Annotations_SampleAttributesDS.txt', header = T,
                               strip.white = T, stringsAsFactors = F, sep = '\t')
  gtex.attributes <- as.data.frame(gtex.attributes)
  gtex.attributes %>% select(SAMPID, SMTS, SMTSD ) %>% View()
  
  tissue.name <- 'Muscle - Skeletal'  #'Thyroid' #'Adipose - Subcutaneous' #'Muscle - Skeletal'  #'Pancreas'
  samples.ids.tissue <- gtex.attributes$SAMPID[grepl(pattern = tissue.name, gtex.attributes$SMTSD)]
  
  length(samples.ids.tissue)
  
  
  colnames(gtex.gct) <- gsub(pattern = '[.]', '-', colnames(gtex.gct))
  
  table(samples.ids.tissue %in% colnames(gtex.gct))
  
  samples.ids.tissue <- samples.ids.tissue[samples.ids.tissue %in% colnames(gtex.gct)]
  
  tissue.mrna.count <- gtex.gct[, samples.ids.tissue]
  rownames(tissue.mrna.count) <- gtex.gct$Name
  
  tissue.mrna.count[1:4, 1:4]
  
  # keep <- apply(tissue.mrna.count, 1, function(x) {
  #   a <- sum(x >= 6)
  #   return(a >= 0.1*length(x))
  # })
  # 
  # table(keep)
  # 
  # tissue.mrna.count <- tissue.mrna.count[keep, ]
  
  
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = tissue.mrna.count, colData = as.data.frame(colnames(tissue.mrna.count)), design = ~ 1)
  dds <- estimateSizeFactors(dds)
  tissue.mrna.count.DESeqnormalized <- counts(dds, normalized=TRUE)
  
  
  dim(tissue.mrna.count.DESeqnormalized)
  
  library(pbmcapply)
  
  sizefactors <- sizeFactors(dds)
  head(sizefactors)
  
  mrna.count.normalized.mu <- apply(tissue.mrna.count.DESeqnormalized[, ], 1, mean)
  mrna.count.normalized.var <- apply(tissue.mrna.count.DESeqnormalized[, ], 1, var)
  
  length(mrna.count.normalized.mu)
  
  s_hat <- psych::geometric.mean(sizefactors)
  dummy.rows <- 2000  
  dummy.cols <- ncol(tissue.mrna.count.DESeqnormalized)
  dummy.mrna.count <- matrix(nrow = dummy.rows, ncol = dummy.cols)  
  rownames(dummy.mrna.count) <- paste0("dummy", seq(1, dummy.rows))
  
  for (i in seq(1:nrow(dummy.mrna.count))) {
    lambda <- s_hat * sample(mrna.count.normalized.mu, size = 1)
    dummy.mrna.count[i, ] <- rpois(n = ncol(dummy.mrna.count), lambda = lambda)
  }
  
  dummy.mrna.count.mu <- apply(dummy.mrna.count, 1, mean)
  dummy.mrna.count.var <- apply(dummy.mrna.count, 1, var)
  
  mrna.count.normalized.and.dummy <- rbind(tissue.mrna.count.DESeqnormalized, dummy.mrna.count)
  mrna.count.normalized.and.dummy.mu <- apply(mrna.count.normalized.and.dummy, 1, mean)
  mrna.count.normalized.and.dummy.var <- apply(mrna.count.normalized.and.dummy, 1, var)
  plot.data <- data.frame(mu = mrna.count.normalized.and.dummy.mu, var = mrna.count.normalized.and.dummy.var)
  plot.data$is_dummy <- grepl("dummy*", rownames(mrna.count.normalized.and.dummy))    
  
  pcount <- ggplot(data = plot.data , aes(x = log2(mu), y = log2(var))) + 
    geom_point(aes(color = as.factor(is_dummy)))   + 
    geom_abline(intercept = 0, slope = 1, color = 'red') + 
    ggtitle('DESeq Normarlized counts \n (Muscle Skeletal)') +
    labs(x = 'log2 mean', y = 'log2 variance', colour = 'Poisson Sampled' )
  
  
  pcount <- pcount + 
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  pcount
  
  source('../../scripts/savefigures.R')
  save.manuscript(plot = pcount, filename = 'muscleskeletal.deseqnormalised.count')
  
  mrna.count.normalized.and.dummy.log2.transformed <- log2(mrna.count.normalized.and.dummy + 0.5)
  plot.data.normalized <- data.frame(mu = apply(mrna.count.normalized.and.dummy.log2.transformed, 1, mean), 
                                     var = apply(mrna.count.normalized.and.dummy.log2.transformed, 1, var))        
  
  plot.data.normalized$is_dummy <- grepl("dummy*", rownames(mrna.count.normalized.and.dummy.log2.transformed))
  
  ptranformed <- ggplot(data = plot.data.normalized, aes(x = mu,y = var)) + 
    geom_point(aes(color = as.factor(is_dummy))) + 
    ggtitle('log2 transformed DESeq \n Normarlized counts muscle skeletal') +
    labs(x = 'mean', y = 'variance', colour = 'Poisson Sampled')
  
  ptranformed <- ptranformed + 
    theme_bw() + 
    theme(text = element_text(size=20), panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
          axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
  
  ptranformed
  save.manuscript(plot = ptranformed, filename = 'muscleskeletal.log2.tranformed.Deseq')
  
  
}


meqtl.cit.ecit <- function() {
  rm(list = ls())
  library(tidyr)
  library(cit)
  
  library(parallel)
  
   
  
  # EXTRACT TRIOS -----------------------------------------------------------
  # step1 -------------------------------------------------------------------
  # to get sig cis_egenes and L
  # tar -tf ../GTEx_Analysis_v7_eQTL.tar.gz | grep Adipose
  # tar -xf ../GTEx_Analysis_v7_eQTL.tar.gz GTEx_Analysis_v7_eQTL/Whole_Blood.v7.egenes.txt.gz
  
  
  
  
  # step2 -------------------------------------------------------------------
  #get genotype on these L
  # vcftools --gzvcf ../../GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz \
  # --out Muscle --extract-FORMAT-info GT --snps L.txt
  #time takes 40min
  # sed -ic  -e 's/0\/0/0/g' -e 's/0\/1/1/g' -e 's/1\/0/1/g' -e's/1\/1/2/g' Pancreas.GT.FORMAT 
  
  
  
  # step3 -------------------------------------------------------------------
  # make genotype, expression, covariates, snps loca, geneloc formatted as matrixeqtl
  # mkdir matrixeqtl
  # zcat ../GTEx_Analysis_v7_eQTL_expression_matrices/Whole_Blood.v7.normalized_expression.bed.gz >> matrixeqtl/expression.txt
  # 4.1- awk '{print $4 "\t" $1 "\t" $2 "\t" $3}' matrixeqtl/expression.txt > matrixeqtl/geneloc.txt then nano geneloc.txt
  # rename expression expression2 matrixeqtl/expression.txt
  # cut -f 4- matrixeqtl/expression2.txt > matrixeqtl/expression.txt
  # 4.2 sed  's/\t/_/1' Muscle.GT.FORMAT > matrixeqtl/genotype.txt
  # 4.3 cat ../GTEx_Analysis_v7_eQTL_covariates/Muscle_Skeletal.v7.covariates.txt > matrixeqtl/covariates.txt
  # 4.4 snpsloc awk '{print $1 "_" $2 "\t" $1 "\t" $2}' Muscle.GT.FORMAT > matrixeqtl/snpsloc.txt
  
  
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
  
  print(path.meqtl.expr)
  write.table(x = expression, file = path.meqtl.expr,
              sep = '\t', quote = F, row.names = F)
  
  path.meqtl.geno
  write.table(x = genotype, file = path.meqtl.geno,
              sep = '\t', quote = F, row.names = F)
  
  path.meqtl.cov
  write.table(x = covariates, file = path.meqtl.cov,
              sep = '\t', quote = F, row.names = F)
  
  path.meqtl.genloc
  write.table(x = geneloc, file = path.meqtl.genloc,
              sep = '\t', quote = F, row.names = F)
  
  path.meqtl.snploc
  write.table(x = snpsloc, file = path.meqtl.snploc,
              sep = '\t', quote = F, row.names = F)
  
  
  rm(list = ls())
  # Matrix eQTL by Andrey A. Shabalin
  # http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
  #
  # Be sure to use an up to date version of R and Matrix eQTL.
  
  # source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
  library(MatrixEQTL)
  
  ## Location of the package with the data files.
  # base.dir = "/data/private-data/GTEX/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/rahul2/Skin_Sun_Exposed/matrixeqtl/"
  base.dir = 'Muscle-Skeletal/Muscle Skeletal/matrixeqtl/'
  path.meqtl.result <- 'Muscle-Skeletal/Muscle Skeletal/matrixeqtl/matrix.cis.trans.rds'
  # setwd(base.dir)
  # base.dir = '.';
  
  ## Settings
  
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Genotype file name
  SNP_file_name = paste(base.dir, "genotype.txt", sep="");
  snps_location_file_name = paste(base.dir, "snpsloc.txt", sep="");
  
  # Gene expression file name
  expression_file_name = paste(base.dir, "expression.txt", sep="");
  gene_location_file_name = paste(base.dir, "geneloc.txt", sep="");
  
  # Covariates file name
  # Set to character() for no covariates
  covariates_file_name = paste(base.dir, "covariates.txt", sep="");
  
  # Output file name
  output_file_name_cis = tempfile();
  output_file_name_tra = tempfile();
  
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis = 2e-2;
  pvOutputThreshold_tra = 1e-5;
  
  # Error covariance matrix
  # Set to numeric() for identity.
  errorCovariance = numeric();
  # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
  
  # Distance for local gene-SNP pairs
  cisDist = 1e6;
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "./."; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  
  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  colnames(snpspos) <- c('snp', 'chr', 'pos')

  snpspos$chr <- paste('chr', snpspos$chr, sep='')
  
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  colnames(genepos) <- c('geneid',	'chr',	'left',	'right')  
  
  head(snpspos)
  head(genepos)
  
  all(snpspos$chr %in% genepos$chr)
  
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE);
  
  unlink(output_file_name_tra);
  unlink(output_file_name_cis);
  
  ## Results:
  
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
  cat('Detected local eQTLs:', '\n');
  show(me$cis$eqtls)
  cat('Detected distant eQTLs:', '\n');
  show(me$trans$eqtls)
  
  ## Plot the Q-Q plot of local and distant p-values
  
  plot(me)
  
  path.meqtl.result
  saveRDS(me, file = path.meqtl.result)
  
  # check some weak and strong eqtl-gene pairs
  expression <- as.matrix(gene)
  genotype <- as.matrix(snps)
  cov <- as.matrix(cvrt)
  
  head(me$trans$eqtls %>% arrange(FDR))
  tail(me$trans$eqtls %>% arrange(FDR))
  L <- as.factor(genotype['12_2590414', ])
  A <- as.numeric(expression['ENSG00000227081', ])
  
  model <- lm(A ~ t(cov))
  A.residue <- residuals(model) + model$coefficients[1]
  
  ggplot() + geom_boxplot(aes(x = as.factor(L), y = A.residue))
  
  anova(lm(A.residue ~ L))
  
  #Part2: Down analysis for best snp for tran gene
  #rm(list = ls())
  library(dplyr)
  # cis <- me$cis$eqtls %>% group_by(gene) %>% summarise(snps=snps[which.min(pvalue)],
  #                                                          statistic=statistic[which.min(pvalue)], pvalue=min(pvalue),
  #                                                          FDR=FDR[which.min(pvalue)], beta=beta[which.min(pvalue)])
  path.gencode <- 'gencode.v19.gene.id.names.map.txt'
  gencode <- read.csv2(path.gencode, header = T, sep = '\t', check.names = F,
                       strip.white = T, stringsAsFactors = F)
  gencode$gene_id <- tstrsplit(gencode$gene_id, '[.]')[[1]]
  rownames(gencode) <- gencode$gene_id
  
  
  gencode$chr <- tstrsplit(gencode$Chromosome, ':')[[1]]
  gencode$start <- tstrsplit(tstrsplit(gencode$Chromosome, ':')[[2]] , '-')[[1]]
  gencode$end <- tstrsplit(tstrsplit(gencode$Chromosome, ':')[[2]] , '-')[[2]]
  
  path.cis.egenes <- "Muscle-Skeletal/Muscle Skeletal/qval_cis_egenes.txt"
  path.meqtl.result <- 'Muscle-Skeletal/Muscle Skeletal/matrixeqtl/matrix.cis.trans.rds'
  path.trios <- 'Muscle-Skeletal/Muscle Skeletal/trios.rds'
  
  
  path.cis.egenes
  cis <- read.table(path.cis.egenes, header = T, check.names = F,
                    strip.white = T, stringsAsFactors = F)
  head(cis)
  
  cis <- cis %>% select(variant_id, gene_id, qval.cis)
  cis$variant_id <- paste(tstrsplit(cis$variant_id, '_')[[1]], tstrsplit(cis$variant_id, '_')[[2]], sep = '_')
  
  
  library(data.table)
  cis$gene_id <- tstrsplit(cis$gene_id, '[.]')[[1]]
  
  head(cis)
  cis <- merge(cis, gencode %>% select(gene_id, chr, start, end), by.x= 'gene_id', by.y = 'gene_id')
  cis$variant_id.chr <- tstrsplit(cis$variant_id, '_')[[1]]
  cis$variant_id.pos <- tstrsplit(cis$variant_id, '_')[[2]]
  cis$variant_id.pos <- as.numeric(cis$variant_id.pos)
  cis$start <- as.numeric(cis$start)
  
  cis$end <- as.numeric(cis$end)
  cis$tss <- pmin(abs(cis$variant_id.pos - cis$start), abs(cis$variant_id.pos -  cis$end))
  table(cis$tss < 1e6)

  
  path.meqtl.result
  me <- readRDS(path.meqtl.result)
  
  trans <- me$trans$eqtls
  trans <- trans %>% filter(FDR < 0.05)
  
  trans <- trans %>% select(variant_id = snps, gene_id = gene, FDR.y = FDR)
  
  head(cis)
  head(trans)
  trans <- merge(trans, gencode %>% select(gene_id, chr, start, end), by.x = 'gene_id', by.y = 'gene_id')
  trans$variant_id.chr <- tstrsplit(trans$variant_id, '_')[[1]]
  trans$variant_id.pos <- tstrsplit(trans$variant_id, '_')[[2]]
  trans$variant_id.pos <- as.numeric(trans$variant_id.pos)
  trans$start <- as.numeric(trans$start)
  trans$end <- as.numeric(trans$end)
  
  trans$tss <- pmin(abs(trans$variant_id.pos - trans$start), abs(trans$variant_id.pos -  trans$end))
  table(trans$tss < 1e6)
  
  trios <- merge(cis, trans, by.x = 'variant_id', by.y = 'variant_id')
  
  head(trios) %>% View()
  table(duplicated(trios %>% select(variant_id, gene_id.x, gene_id.y)))
  
  thres <- 1e-5
  table(trios$qval.cis < thres & trios$FDR.y < thres)
  
  trios <- trios %>% select(variant_id, gene_id.x, gene_id.y, qval.cis, FDR.y)
  
  table(trios$gene_id.x == trios$gene_id.y)
  
  path.trios
  saveRDS(trios, path.trios)
  dim(trios)
  
  
  # Part3 ECIT analysis ---------------------analysis results kept in the manuscript--------------------------------
  rm(list = ls())
  
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
  
  
  path.trios
  trios <- readRDS(path.trios)
  
  path.tiss1.me; path.tiss2.me
  tis1.me.var <- readRDS(path.tiss1.me)
  rownames(tis1.me.var) <- tis1.me.var$gene.id
  tis2.me.var <- readRDS(path.tiss2.me)
  rownames(tis2.me.var) <- tis2.me.var$gene.id
  
  head(tis1.log2.expr.resi[, 1:4])  
  head(tis2.log2.expr.resi[, 1:4])
  head(trios)
  head(tis1.me.var)
  head(tis2.me.var)
  head(genotype[, 1:4])
  
  # align all samples 
  indv.tis1 <- colnames(tis1.log2.expr.resi)
  indv.tis2 <- colnames(tis2.log2.expr.resi)
  indv.geno <- colnames(genotype)
  indv.com <- intersect(indv.tis1, intersect(indv.tis2, indv.geno))
  
  
  
  length(indv.tis1)
  
  
  length(indv.tis2)
  length(indv.geno)
  length(indv.com)
  
  tis1.log2.expr.resi <- tis1.log2.expr.resi[, indv.com]
  
  
  tis2.log2.expr.resi <- tis2.log2.expr.resi[, indv.com]
  genotype <- genotype[, indv.com]
  
  dim(tis1.log2.expr.resi)
  dim(tis2.log2.expr.resi)
  dim(genotype)
  
  table(trios$gene_id.x %in% rownames(tis1.log2.expr.resi))
  table(trios$gene_id.y %in% rownames(tis2.log2.expr.resi))
  table(trios$variant_id %in% rownames(genotype))
  
  source('../../scripts/modelstage2/adj_cit1_4.R')
  library(cit)
  getcorrelation <- function(L, A, B) {
    return(c('rho.LA' = cor(L, A, method = 'spearman'), 'rho.LB' = cor(L, B, method = 'spearman'), 'rho.AB' = cor(A, B, method = 'spearman')))
  }
  
  row.i <- 641
  
  rownames(trios) <- 1:nrow(trios)
  dim(trios)
  
  library(pbmcapply)
  
  result.trios <- pbmclapply(1:nrow(trios[ , ]), FUN = function(row.i) {
    if(row.i %% 100 == 0) {
      print(trios[row.i, ])
      
    }
    
    L <- as.integer(genotype[trios$variant_id[row.i], ])
    
    indv.not.na <- !is.na(L)
    
    L <- L[indv.not.na]
    
    
    L1 <- rep(1, length(L))
    L2 <- rep(1, length(L))  
    
    L1[L == 0] <- 0
    L2[L == 0] <- 0
    
    L1[L == 1] <- 0
    L2[L == 1] <- 1
    
    A <- as.numeric(tis1.log2.expr.resi[trios$gene_id.x[row.i], indv.not.na])
    B <- as.numeric(tis2.log2.expr.resi[trios$gene_id.y[row.i], indv.not.na])
    
    me.var.A <- tis1.me.var$me.var.log2.estimated[match(trios$gene_id.x[row.i], rownames(tis1.me.var))]  
    me.var.B <- tis2.me.var$me.var.log2.estimated[match(trios$gene_id.y[row.i], rownames(tis2.me.var))]
    
    result.cit.AB <- tryCatch({
      cit.cp(L = cbind(L1, L2), G = A, T = B, n.resampl = 100)
    },
    error = function(cond) {
      p_vals <- rep(NA, 5)
      names(p_vals) <- c("p_cit", "p_TassocL", "p_TassocGgvnL",  "p_GassocLgvnT",  "p_LindTgvnG" )
      return(p_vals)
    })
    
    result.cit.BA <- tryCatch({
      cit.cp(L = cbind(L1, L2), G = B, T = A, n.resampl = 100)
    },
    error = function(cond) {
      p_vals <- rep(NA, 5)
      names(p_vals) <- c("p_cit", "p_TassocL", "p_TassocGgvnL",  "p_GassocLgvnT",  "p_LindTgvnG" )
      return(p_vals)
    })
    
    result.ecit.AB <- get_adj_cit_pvals(L1 = L1, L2 = L2, Gp = A, Tp = B,
                                     v_eG = me.var.A, v_eT = me.var.B, 
                                     bootstrap = 1000, resampl = 100)

    result.ecit.BA <- get_adj_cit_pvals(L1 = L1, L2 = L2, Gp = B, Tp = A,
                                        v_eG = me.var.B, v_eT = me.var.A, 
                                        bootstrap = 1000, resampl = 100)
    
    cor.AB.result <- getcorrelation(L, A, B)
    
    result <- c(me.var.A = me.var.A, me.var.B = me.var.B, cit.AB = result.cit.AB, ecit.AB = result.ecit.AB, 
                cit.BA = result.cit.BA, ecit.BA = result.ecit.BA, cor.AB.result)
    
    return(result)
    # temp <- 1
  }, mc.cores = 80, ignore.interactive = F)
  
  result.trios <- as.data.frame(do.call(rbind, result.trios))
  
  trios <- cbind(trios, result.trios)
  
  path.cit.ecit
  saveRDS(trios, path.cit.ecit)
  
  # Result Analysis
  rm(list = ls())
  source('../../scripts/modelstage2/adj_cit1_4.R')
  path.cit.ecit <- 'Muscle-Skeletal/Muscle Skeletal/trios.cit.ecit.result.rds'
  all.trios <- readRDS(path.cit.ecit)
  
  trios <- all.trios
  
  gene.names <- function() {
    biomart <- read.table('biomart.GRCh37.p13.txt', header = F, strip.white = T, stringsAsFactors = F)
    biomart <- biomart %>% select(gene.id = V1, gene.symbol = V3, gene.type = V4)
    head(biomart)
    
    trios <- merge(trios, biomart, by.x = 'gene_id.x', by.y = 'gene.id', all.x = T)
    trios <- merge(trios, biomart, by.x = 'gene_id.y', by.y = 'gene.id', all.x = T)
    
    
    temp <- 1
  }
  
  mappability.analysis <- function() {
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
    
  }
  
  trios$cit.AB.p_cit.BH <- p.adjust(trios$cit.AB.p_cit, method = 'BH')
  trios$ecit.AB.adj_p_cit.BH <- p.adjust(trios$ecit.AB.adj_p_cit, method = 'BH')
  trios$cit.BA.p_cit.BH <- p.adjust(trios$cit.BA.p_cit, method = 'BH')
  trios$ecit.BA.adj_p_cit.BH <- p.adjust(trios$ecit.BA.adj_p_cit, method = 'BH')
  
  trios.backup <- trios
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
  
  cytoscape()
  boxplot.analysis <- function() {
    # cuasal example
    L.id <- '13_74110412'
    A.id <- 'ENSG00000102554'
    B.id <- 'ENSG00000198324'
    
    # independent example 16_81080366
    L.id <- '11_9242168'
    A.id <- 'ENSG00000184014'
    B.id <- 'ENSG00000106089'
    
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
    
    {
    # only for indep example remove outlier
    outlier <- A > 10 | A < 8.5
    table(outlier)
    L <- L[!outlier]
    L1 <- L1[!outlier]
    L2 <- L2[!outlier]  
    A <- A[!outlier]
    B <- B[!outlier]
    }
    
    table(L)
    table(L1, L2)
    
    me.var.A <- tis1.me.var$me.var.log2.estimated[match(A.id, tis1.me.var$gene.id)]
    me.var.B <- tis2.me.var$me.var.log2.estimated[match(B.id, tis2.me.var$gene.id)]
    
    cit.cp(cbind(L1, L2), A, B, n.resampl = 100)
    cit.cp(cbind(L1, L2), B, A, n.resampl = 100)
    
    get_adj_cit_pvals(L1 = L1, L2 = L2, Gp = A, Tp = B,
                      v_eG = me.var.A, v_eT = me.var.B, 
                      bootstrap = 1000, resampl = 100)
    
    get_adj_cit_pvals(L1 = L1, L2 = L2, Gp = B, Tp = A, 
                      v_eG = me.var.B, v_eT = me.var.A,
                      bootstrap = 1000, resampl = 100)
    
    AB <- qplot(A, B) + geom_smooth() + xlab(A.symbol) + ylab(B.symbol) +
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
    
    LB <- ggplot()  + geom_boxplot(aes(as.factor(L), B), notch = T) + xlab(L.id) + ylab(B.symbol) +
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")# L ~ B
    
    LA_B <- ggplot() + geom_boxplot(aes(as.factor(L), A), notch = T) + xlab(L.id) + 
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
    
    
    L.ind.B_A <- ggplot() + geom_boxplot(aes(as.factor(L), residuals(lm(B~A))), notch = T) + 
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
  
  GSEA.analysis <- function() {
    
    row.i <- 1
    
    gsea.result <- pbmclapply(X = 1:nrow(trios), FUN = function(row.i) {
      trios[row.i, ]
      AB.pval <- as.numeric(trios$ecit.AB.adj_p_cit[row.i])
      BA.pval <- as.numeric(trios$ecit.BA.adj_p_cit[row.i])
      
      result <- NULL
      if(AB.pval <= BA.pval) {
        result <- c('gene.id' = trios$gene_id.x[row.i], 'gene.symbol' = trios$gene.symbol.x[row.i], 
                    'pval' = AB.pval, ecit.direction = trios$ecit.direction[row.i])
      }else {
        result <- c('gene_id' = trios$gene_id.y[row.i], 'gene.symbol' = trios$gene.symbol.y[row.i], 
                    'pval' = BA.pval, ecit.direction = trios$ecit.direction[row.i])
      }
      
      return(result)
    }, mc.cores = 1)
    
    gsea.result <- as.data.frame(do.call(rbind, gsea.result))
    gsea.result$pval <- as.numeric(gsea.result$pval)
    gsea.result$score <- -log10(gsea.result$pval)
    
    temp <- gsea.result
    gsea.result <- gsea.result %>% group_by(gene_id) %>% dplyr::slice_max(order_by = score, n = 1)
    
    gsea.result <- as.data.frame(gsea.result)
    head(gsea.result)
    gsea.result <- gsea.result %>% select(gene.symbol, score) %>% arrange(desc(score))
    
    gsea.result <- na.omit(gsea.result)
    
    write.table(gsea.result, file = 'Muscle - Skeletal/Muscle Skeletal/gsea.regulators.txt', col.names = F,
                row.names = F, quote = F, sep = '\t')
    
    temp <- 1
  }
  

  temp <- 1
  
}

common.samples.list()
deseq.count.deseq.resi.deseq.me.var.all()
deseq.mean.var.plots.generate()
meqtl.cit.ecit()
