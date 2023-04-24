rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("scripts/modelstage2/adj_cit1_4.yeast.R")
# source("scripts/simulation/data_generate_causal.R")
# path.output <- "scripts/simulation/result_yeast.data_generate_causal3.rds"
source('scripts/simulation/data_generate_confounding.R')
path.output <- "scripts/simulation/result_yeast.data_generate_confounding3.rds"

# source('scripts/simulation/data_generate_causal_with_confounding.R')
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


parameters <- expand.grid(
  #hemani parameter
  n = c(300, 500, 1000),
  p = 0.5,
  r_ab = sqrt(seq(0, 1, by=0.2)),
  r_za = c(sqrt(0.1)),
  noisea = sqrt(seq(0, 1, by=0.2)),
  noiseb = sqrt(seq(0, 1, by=0.2)),
  nsim = 1:100
)

# parameters <- expand.grid(
#   n = c(300, 500, 1000),
#   p = 0.5,
#   r_za = sqrt(c(0.16, 0.49)),
#   r_ab = sqrt(c(0.00, 0.49)),
#   noisea = sqrt(seq(0, 1, by=0.2)),
#   # noisea = sqrt(c(0.6)),
#   noiseb = sqrt(seq(0, 1, by=0.2)),
#   nsim = 1:100
# )

i <- 10

run <- function(i)
{
  resampl <- 100
  message(i, " of ", nrow(parameters))
  dat <- with(parameters[i,], make_system_yeast(n, p, r_ab, r_za, noisea, noiseb))
  
  pval_AB <- cit.cp(L = dat$L, G = dat$Ap, T = dat$Bp, n.resampl = resampl)
  pval_BA <- cit.cp(L = dat$L, G = dat$Bp, T = dat$Ap, n.resampl = resampl)
  
  cit_res <- get_cit_direction(pval_AB[1], pval_BA[1],  thresh = 0.05)
  
  pvaladj_AB <- get_adj_cit_pvals(L = dat$L, Gp = dat$Ap, Tp = dat$Bp,
                                  v_eG = parameters[i,]$noisea^2, v_eT = parameters[i,]$noiseb^2, bootstrap = 1000 ,resampl = resampl)
  pvaladj_BA <- get_adj_cit_pvals(L = dat$L, Gp = dat$Bp, Tp = dat$Ap,
                                  v_eG = parameters[i,]$noiseb^2, v_eT = parameters[i,]$noisea^2, bootstrap = 1000 ,resampl = resampl)
  
  adj_cit_res <- get_cit_direction(pvaladj_AB[1], pvaladj_BA[1], thresh = 0.05)
  
  return(c(AB=pval_AB, BA=pval_BA, AB=pvaladj_AB, BA=pvaladj_BA, cit_res=cit_res, adjcit_res=adj_cit_res, var_G= var(dat$Ap), var_T=var(dat$Bp)))
  
}

library(pbmcapply)

result <- pbmclapply(FUN = run, X=1:nrow(parameters), mc.cores = 80, ignore.interactive = T)

# result <- mcmapply(run, i=1:nrow(parameters), mc.cores = 100)
parameters <- cbind(parameters, t(as.data.frame(result)))

#eibcccundbnvkvgcdktdcdhkrkneikignlbvuueibict
if (file.exists(path.output)) {
  message('change file name')
}

path.output
write(paste0('\n\n', path.output), file = 'scripts/simulation/remarks.txt', append = T)

for (col in colnames(parameters)[1:7]) {
  print(col)
  write(paste0(col, unique(parameters[, col]), collapse = '', sep = '; '), file = 'scripts/simulation/remarks.txt', append = T)
  
}
path.output
saveRDS(parameters, file = path.output)

# fig 2(A,B,C)
causal.model.plot.analysis <- function() {
  # # my way of plotting tpr i'e power and fpr ------------------------------
  # fig 2(A)
  plotting.tpr.fpr.v2 <- function() {
    rm(list = ls())
    
    library(RColorBrewer)
    source("scripts/modelstage2/adj_cit1_4.yeast.R")
    result.causal <- readRDS('scripts/simulation/result_yeast.data_generate_causal3.rds') #readRDS('scripts/simulation/result_yeast.data_generate_causal1.rds')
    str(result.causal)
    result.causal$model <- 1
    thres <- 0.05
    result.causal$predict.cit <- get_cit_direction(result.causal$AB.p_cit, result.causal$BA.p_cit, thres)
    result.causal$predict.ecit <- get_cit_direction(result.causal$AB.adj_p_cit, result.causal$BA.adj_p_cit, thres)
    rownames(result.causal) <- seq(1:nrow(result.causal))
    
    tp <- result.causal %>% group_by(n, noisea, noiseb, r_ab, r_za) %>% 
      dplyr::summarise(CIT = sum(as.logical(model == predict.cit)),
                       ECIT = sum(as.logical(model == predict.ecit)))
    
    library(reshape2)
    tp.melt <- melt(data = tp, id.vars = c('n', 'noisea', 'noiseb', 'r_ab', 'r_za'),
                    variable.name = 'method',
                    value.name = 'tp')
    
    tp.melt$lab <- 'Power'
    tp.melt$lab[abs(tp.melt$r_ab - 0.0) < 0.001] <- 'Type 1 error'
    tp.melt$proportion.causality.predicted <- tp.melt$tp / 100
    tp.melt$proportion.causality.predicted <- round(tp.melt$proportion.causality.predicted, digits = 2)
    
    tp.melt$n <- as.factor(tp.melt$n)
    levels(tp.melt$n)
    levels(tp.melt$n) <- paste("n=", levels(tp.melt$n), sep = '')
    levels(tp.melt$n)
    
    unique(result.causal$r_ab)^2
    plot.r_ab <- sqrt(0.4)
    
    unique(result.causal$noiseb)^2
    plot.noiseb <- sqrt(0.4)
    
    p.tpr.fpr <- tp.melt %>% filter(abs(r_ab - 0.0) < 0.001 | abs(r_ab - sqrt(0.4)) < 0.001, abs(noiseb - plot.noiseb) < 0.001) %>% 
      ggplot(aes(x = noisea^2, y = round(proportion.causality.predicted, 2), fill = method)) + 
      geom_bar(stat = 'identity', position = 'dodge') +
      facet_grid(lab ~ n) + 
      labs(y="Discovery rate", x = expression(sigma[eg]^2))
    
    p.tpr.fpr <- p.tpr.fpr +
      theme_bw() + 
      theme(text = element_text(size=15), panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
    
    # fig 2(A)
    print(p.tpr.fpr)
    
    for (exten in c('.pdf', '.jpeg')) {
      units = "in"
      height = 4.5
      width = 5.7
      ggsave(plot = p.tpr.fpr,
             units = units,
             height = height,
             width = width,
             filename = paste("output/p.tpr.fpr.simulation", exten, sep = ''))
      
      ggsave(plot = p.tpr.fpr,
             units = units,
             height = height,
             width = width,
             filename = paste("output/p.tpr.fpr.simulation", exten, sep = ''))
      
    }
    
    
  }
  plotting.tpr.fpr.v2()
  # fig 2(B)
  delta.tpr.matrix <- function() {
    rm(list = ls())
    library(RColorBrewer)
    source("scripts/modelstage2/adj_cit1_4.yeast.R")
    result.causal <- readRDS('scripts/simulation/result_yeast.data_generate_causal3.rds') #readRDS('scripts/simulation/result_yeast.data_generate_causal1.rds')
    str(result.causal)
    result.causal$model <- 1
    thres <- 0.05
    result.causal$predict.cit <- get_cit_direction(result.causal$AB.p_cit, result.causal$BA.p_cit, thres)
    result.causal$predict.ecit <- get_cit_direction(result.causal$AB.adj_p_cit, result.causal$BA.adj_p_cit, thres)
    rownames(result.causal) <- seq(1:nrow(result.causal))
    
    tp <- result.causal %>% group_by(n, noisea, noiseb, r_ab, r_za) %>% 
      dplyr::summarise(CIT = sum(as.logical(model == predict.cit)),
                       ECIT = sum(as.logical(model == predict.ecit)))
    tpr.confusion <- tp %>% filter(n == 500 & abs(r_ab - sqrt(0.4)) < 0.001)
    
    tpr.confusion$delta <- tpr.confusion$ECIT - tpr.confusion$CIT
    
    # tpr.confusion.matrix <- matrix(nrow = length(unique(tpr.confusion$noisea)), ncol = length(unique(tpr.confusion$noiseb)))
    # 
    # rownames(tpr.confusion.matrix) <- unique(tpr.confusion$noiseb^2)
    # 
    # colnames(tpr.confusion.matrix) <- unique(tpr.confusion$noisea^2)
    # 
    # tpr.confusion.matrix
    # 
    # for (row in rownames(tpr.confusion.matrix)) {
    #   for (col in colnames(tpr.confusion.matrix)) {
    #     
    #     # cat(row,',' ,col, ' ')
    #     tpr.confusion.matrix[row, col] <- tpr.confusion$delta[abs(tpr.confusion$noiseb^2 - as.numeric(row)) < 0.001 &
    #                                                             abs(tpr.confusion$noisea^2 - as.numeric(col)) < 0.001]
    #     cat(temp, ',')
    #     
    #   }
    #   cat('\n')
    # }
    # 
    # plot(tpr.confusion.matrix)
    
    tpr.confusion$delta <- tpr.confusion$delta/100
    
    tpr.confusion.plot <- tpr.confusion %>% ggplot(aes(noisea^2, noiseb^2, fill = delta)) + 
      geom_tile(color='gray') +
      geom_text(aes(label = delta)) +
      scale_fill_gradient(low = "white", high = "red") +
      labs(y=expression(sigma[et]^2), x = expression(sigma[eg]^2), fill = expression(Delta[TPR]))
    
    tpr.confusion.plot <- tpr.confusion.plot + 
      theme_bw() + 
      theme(text = element_text(size=15), panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
    
    # manuscript fig2(B)
    print(tpr.confusion.plot)
    
    filename <- 'tpr.matrix.simulation'
    
    for (path in c("output/")) {
      plot.save.manuscript(tpr.confusion.plot, path,
                           filename = filename)
      
    }
    
    # for (exten in c('.pdf', '.jpeg')) {
    #   units = "in"
    #   height = 4.5
    #   width = 5.7
    #   ggsave(plot = tpr.confusion.plot,
    #          units = units,
    #          height = height,
    #          width = width,
    #          filename = paste("/data/users/cs18s008/projects/ecit/manuscript/output/tpr.matrix.simulation", exten, sep = ''))
    #   
    #   ggsave(plot = tpr.confusion.plot,
    #          units = units,
    #          height = height,
    #          width = width,
    #          filename = paste("/data/users/cs18s008/projects/ecit/ECIT-ManuscriptFirst-Draft/figures/tpr.matrix.simulation", exten, sep = ''))
    #   
    # }
    
    temp <- 1  
  }
  delta.tpr.matrix()
  
  
  # fig 2(C)
  delta.type1.error.matrix <- function() {
    rm(list = ls())
    library(RColorBrewer)
    source("scripts/modelstage2/adj_cit1_4.yeast.R")
    result.causal <- readRDS('scripts/simulation/result_yeast.data_generate_causal3.rds') #readRDS('scripts/simulation/result_yeast.data_generate_causal1.rds')
    str(result.causal)
    result.causal$model <- 1
    thres <- 0.05
    result.causal$predict.cit <- get_cit_direction(result.causal$AB.p_cit, result.causal$BA.p_cit, thres)
    result.causal$predict.ecit <- get_cit_direction(result.causal$AB.adj_p_cit, result.causal$BA.adj_p_cit, thres)
    rownames(result.causal) <- seq(1:nrow(result.causal))
    
    tp <- result.causal %>% group_by(n, noisea, noiseb, r_ab, r_za) %>% 
      dplyr::summarise(CIT = sum(as.logical(model == predict.cit)),
                       ECIT = sum(as.logical(model == predict.ecit)))
    
    type1error <- tp %>% filter(abs(r_ab - 0) < 0.001 & n == 500)
  
    type1error$delta <- (type1error$ECIT - type1error$CIT)/100
    
    type1error.plot <- type1error %>% ggplot(aes(x = noisea^2, y = noiseb^2, fill = delta)) +
      geom_tile(color='gray') +
      geom_text(aes(label = delta)) +
      scale_fill_gradient(low = "white", high = "red") +
      labs(y=expression(sigma[et]^2), x = expression(sigma[eg]^2), fill = expression(Delta[Type1error]))
    
    type1error.plot <- type1error.plot + 
      theme_bw() + 
      theme(text = element_text(size=15), panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
    
    # manuscript 2.(C)
    print(type1error.plot)
    for (path in c("output/")) {
      plot.save.manuscript(type1error.plot, path,
                          filename = 'type1error.matrix')
      
    }
    
    
    temp <- 1
    
  }
  
  delta.type1.error.matrix()
    
}
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
causal.model.plot.analysis()

 independent.model.plot.analysis <- function() {
  # # my way of plotting tpr i'e power and fpr ------------------------------
  plotting.tpr.fpr.v2 <- function() {
    rm(list = ls())
    library(RColorBrewer)
    source("scripts/modelstage2/adj_cit1_4.yeast.R")
    result.independt <- readRDS('scripts/simulation/result_yeast.data_generate_confounding3.rds')
    str(result.independt)
    result.independt$model <- 3
    thres <- 0.05
    result.independt$predict.cit <- get_cit_direction(result.independt$AB.p_cit, result.independt$BA.p_cit, thres)
    result.independt$predict.ecit <- get_cit_direction(result.independt$AB.adj_p_cit, result.independt$BA.adj_p_cit, thres)
    
    
    # result.causal <- readRDS('scripts/simulation/result_yeast.data_generate_causal3.rds') #readRDS('scripts/simulation/result_yeast.data_generate_causal1.rds')
    # str(result.causal)
    # result.causal$model <- 1
    # thres <- 0.05
    # result.causal$predict.cit <- get_cit_direction(result.causal$AB.p_cit, result.causal$BA.p_cit, thres)
    # result.causal$predict.ecit <- get_cit_direction(result.causal$AB.adj_p_cit, result.causal$BA.adj_p_cit, thres)
    # rownames(result.causal) <- seq(1:nrow(result.causal))
    # 
    tp <- result.independt %>% group_by(n, noisea, noiseb, r_ab, r_za) %>% 
      dplyr::summarise(CIT = sum(as.logical(model == predict.cit)),
                       ECIT = sum(as.logical(model == predict.ecit)))
    
    library(reshape2)
    tp.melt <- melt(data = tp, id.vars = c('n', 'noisea', 'noiseb', 'r_ab', 'r_za'),
                    variable.name = 'method',
                    value.name = 'tp')
    
    tp.melt$lab <- 'Power'
    tp.melt$proportion.predicted <- tp.melt$tp / 100
    tp.melt$proportion.predicted <- round(tp.melt$proportion.predicted, digits = 2)
    
    tp.melt$n <- as.factor(tp.melt$n)
    levels(tp.melt$n)
    levels(tp.melt$n) <- paste("n=", levels(tp.melt$n), sep = '')
    levels(tp.melt$n)
    
    
    result.causal <- readRDS('scripts/simulation/result_yeast.data_generate_causal3.rds')
    result.causal <- result.causal %>% filter(!abs(r_ab - 0.0) < 0.001)
    thres <- 0.05
    result.causal$predict.cit <- get_cit_direction(result.causal$AB.p_cit, result.causal$BA.p_cit, thres)
    result.causal$predict.ecit <- get_cit_direction(result.causal$AB.adj_p_cit, result.causal$BA.adj_p_cit, thres)
    table(result.causal$predict.cit)
    table(result.causal$predict.ecit)
    
    fp <- result.causal %>% group_by(n, noisea, noiseb, r_ab, r_za) %>% 
      dplyr::summarise(CIT = sum(as.logical(predict.cit == 3)),
                       ECIT = sum(as.logical(predict.ecit == 3)))
    
    fp.melt <- melt(data = fp, id.vars = c('n', 'noisea', 'noiseb', 'r_ab', 'r_za'),
                    variable.name = 'method',
                    value.name = 'fp')
    
    fp.melt$lab <- 'Type I error'
    fp.melt$proportion.predicted <- fp.melt$fp / 100
    fp.melt$proportion.predicted <- round(fp.melt$proportion.predicted, digits = 2)
    
    fp.melt$n <- as.factor(fp.melt$n)
    levels(fp.melt$n)
    levels(fp.melt$n) <- paste("n=", levels(fp.melt$n), sep = '')
    levels(fp.melt$n)
    
    discovery <- rbind(tp.melt %>% select(n, noisea, noiseb, r_ab, r_za, method, lab, proportion.predicted),
                       fp.melt %>% select(n, noisea, noiseb, r_ab, r_za, method, lab, proportion.predicted))
    
    discovery$n <- as.factor(discovery$n)
    
    
    unique(result.independt$r_ab)^2
    plot.r_ab <- sqrt(0.4)
    
    unique(result.independt$noiseb)^2
    plot.noiseb <- sqrt(0.4)
    
    p.title <- paste0('Independent Model \n ', expression(sigma[et]^2), '=', round(plot.noiseb^2, 2))
    
    p.tpr.fpr <- discovery %>% filter(abs(r_ab - plot.r_ab) < 0.001, abs(noiseb - plot.noiseb) < 0.001) %>% 
      ggplot(aes(x = noisea^2, y = round(proportion.predicted, 2), fill = method)) + 
      geom_bar(stat = 'identity', position = 'dodge') +
      facet_grid(lab ~ n) + 
      labs(y="Discovery rate", x = expression(sigma[eg]^2), 
           title = p.title) 
    
    p.tpr.fpr <- p.tpr.fpr +
      theme_bw() + 
      theme(text = element_text(size=15), panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
    
    print(p.tpr.fpr)
    
    for (exten in c('.pdf', '.jpeg')) {
      units = "in"
      height = 4.5
      width = 5.7
      # ggsave(plot = p.tpr.fpr,
      #        units = units,
      #        height = height,
      #        width = width,
      #        filename = paste("/data/users/cs18s008/projects/ecit/manuscript/output/p.tpr.fpr.simulation.indep", exten, sep = ''))
      
      ggsave(plot = p.tpr.fpr,
             units = units,
             height = height,
             width = width,
             filename = paste("output/p.tpr.fpr.simulation.indep", exten, sep = ''))
      
    }
    
    
  }
  
  plotting.tpr.fpr.v2()
  
  delta.tpr.matrix <- function() {
    rm(list = ls())
    library(RColorBrewer)
    source("scripts/modelstage2/adj_cit1_4.yeast.R")
    result.causal <- readRDS('scripts/simulation/result_yeast.data_generate_causal3.rds') #readRDS('scripts/simulation/result_yeast.data_generate_causal1.rds')
    str(result.causal)
    result.causal$model <- 1
    thres <- 0.05
    result.causal$predict.cit <- get_cit_direction(result.causal$AB.p_cit, result.causal$BA.p_cit, thres)
    result.causal$predict.ecit <- get_cit_direction(result.causal$AB.adj_p_cit, result.causal$BA.adj_p_cit, thres)
    rownames(result.causal) <- seq(1:nrow(result.causal))
    
    tp <- result.causal %>% group_by(n, noisea, noiseb, r_ab, r_za) %>% 
      dplyr::summarise(CIT = sum(as.logical(model == predict.cit)),
                       ECIT = sum(as.logical(model == predict.ecit)))
    tpr.confusion <- tp %>% filter(n == 500 & abs(r_ab - sqrt(0.4)) < 0.001)
    
    tpr.confusion$delta <- tpr.confusion$ECIT - tpr.confusion$CIT
    
    # tpr.confusion.matrix <- matrix(nrow = length(unique(tpr.confusion$noisea)), ncol = length(unique(tpr.confusion$noiseb)))
    # 
    # rownames(tpr.confusion.matrix) <- unique(tpr.confusion$noiseb^2)
    # 
    # colnames(tpr.confusion.matrix) <- unique(tpr.confusion$noisea^2)
    # 
    # tpr.confusion.matrix
    # 
    # for (row in rownames(tpr.confusion.matrix)) {
    #   for (col in colnames(tpr.confusion.matrix)) {
    #     
    #     # cat(row,',' ,col, ' ')
    #     tpr.confusion.matrix[row, col] <- tpr.confusion$delta[abs(tpr.confusion$noiseb^2 - as.numeric(row)) < 0.001 &
    #                                                             abs(tpr.confusion$noisea^2 - as.numeric(col)) < 0.001]
    #     cat(temp, ',')
    #     
    #   }
    #   cat('\n')
    # }
    # 
    # plot(tpr.confusion.matrix)
    
    tpr.confusion$delta <- tpr.confusion$delta/100
    
    tpr.confusion.plot <- tpr.confusion %>% ggplot(aes(noisea^2, noiseb^2, fill = delta)) + 
      geom_tile(color='gray') +
      geom_text(aes(label = delta)) +
      scale_fill_gradient(low = "white", high = "red") +
      labs(y=expression(sigma[et]^2), x = expression(sigma[eg]^2), fill = expression(Delta[TPR]))
    
    tpr.confusion.plot <- tpr.confusion.plot + 
      theme_bw() + 
      theme(text = element_text(size=15), panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
    
    print(tpr.confusion.plot)
    
    filename <- 'tpr.matrix.simulation'
    
    for (path in c("output/")) {
      plot.save.manuscript(tpr.confusion.plot, path,
                           filename = filename)
      
    }
    
    # for (exten in c('.pdf', '.jpeg')) {
    #   units = "in"
    #   height = 4.5
    #   width = 5.7
    #   ggsave(plot = tpr.confusion.plot,
    #          units = units,
    #          height = height,
    #          width = width,
    #          filename = paste("/data/users/cs18s008/projects/ecit/manuscript/output/tpr.matrix.simulation", exten, sep = ''))
    #   
    #   ggsave(plot = tpr.confusion.plot,
    #          units = units,
    #          height = height,
    #          width = width,
    #          filename = paste("/data/users/cs18s008/projects/ecit/ECIT-ManuscriptFirst-Draft/figures/tpr.matrix.simulation", exten, sep = ''))
    #   
    # }
    
    temp <- 1  
  }
  delta.tpr.matrix()
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
  
  delta.type1.error.matrix <- function() {
    rm(list = ls())
    library(RColorBrewer)
    source("scripts/modelstage2/adj_cit1_4.yeast.R")
    result.causal <- readRDS('scripts/simulation/result_yeast.data_generate_causal3.rds') #readRDS('scripts/simulation/result_yeast.data_generate_causal1.rds')
    str(result.causal)
    result.causal$model <- 1
    thres <- 0.05
    result.causal$predict.cit <- get_cit_direction(result.causal$AB.p_cit, result.causal$BA.p_cit, thres)
    result.causal$predict.ecit <- get_cit_direction(result.causal$AB.adj_p_cit, result.causal$BA.adj_p_cit, thres)
    rownames(result.causal) <- seq(1:nrow(result.causal))
    
    tp <- result.causal %>% group_by(n, noisea, noiseb, r_ab, r_za) %>% 
      dplyr::summarise(CIT = sum(as.logical(model == predict.cit)),
                       ECIT = sum(as.logical(model == predict.ecit)))
    
    type1error <- tp %>% filter(abs(r_ab - 0) < 0.001 & n == 500)
    
    type1error$delta <- (type1error$ECIT - type1error$CIT)/100
    
    type1error.plot <- type1error %>% ggplot(aes(x = noisea^2, y = noiseb^2, fill = delta)) +
      geom_tile(color='gray') +
      geom_text(aes(label = delta)) +
      scale_fill_gradient(low = "white", high = "red") +
      labs(y=expression(sigma[et]^2), x = expression(sigma[eg]^2), fill = expression(Delta[Type1error]))
    
    type1error.plot <- type1error.plot + 
      theme_bw() + 
      theme(text = element_text(size=15), panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color="black"), 
            axis.text.y = element_text(colour = "black"), , legend.position = "bottom", legend.direction = "horizontal")
    
    
    print(type1error.plot)
    for (path in c("output/")) {
      plot.save.manuscript(type1error.plot, path,
                           filename = 'type1error.matrix')
      
    }
    
    
    temp <- 1
    
  }
  
  delta.type1.error.matrix()
  
 }
 independent.model.plot.analysis()
 