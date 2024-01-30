#Sankowski et al -functions and plots

#colors
colors_many <- toupper(c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe","#008080","#e6beff","#aa6e28","#fffac8","#800000","#aaffc3","#808000","#ffd8b1","#000080","#808080","#FFFFFF","#000000"))
colors <- toupper(c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999'))
colors_pat <- toupper(c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5'))
colors_fig <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3', "light grey", "grey", "dark grey", "#696969")

signature_genes <- data.frame("monocytes"=c('CCR2', 'CLEC12A', 'PLAC8', 'FCN1', 'S100A9'),
                              "macrophages"=c("MRC1", "MS4A7", "CD163", "LYVE1", "STAB1"),
                              "microglia"= c('P2RY12','TMEM119', 'SLC2A5',NA,NA),
                              "tcell"=c('TRAC', 'TRBC2', 'CD52', 'IL32', NA),
                              "nk"=c("NKG7","KLRB1","PRF1","GZMB","GZMA"),
                              "myeloid"=c('ITGAM',  'MS4A6A', 'TYROBP', 'CD14', NA),
                              "oligodendrocyte"=c('MBP',  'MOG', 'MAG', 'PLP1', NA),
                              "bcells"=c('CD79A', 'IGHG4', 'IGLL5', NA, NA),
                              "apc"=c("CD74", "CD80", "CD86", "HLA-DRA", "CD40"),
                              "pDC"=c("IL3RA","LILRA4","TCF4","SELL","LTB"),
                              "cDC1"=c("XCR1", "CLEC9A","CADM1", "IRF8","BATF3"),
                              "cDC2"=c("FCER1A", "CLEC10A", "CD1C","CST7","CCR6"),
                              "MigDCs"= c("CCR7", "LAMP3", "SAMSN1",NA,NA),
                              "Endothelial_Cells"=c("HSPG2","PLVAP","FLT1","VWF","CD34"),
                              "astrocyte"=c("GFAP", "HEPACAM","SOX9","AQP4",NA),
                              stringsAsFactors = F)

#plot expression seurat
plot_expmap_seurat <- function(features, object=all, reduction = "umap", dims=c(1,2), point_size=1, logsc=FALSE, line_width=0, .retain_cl = retain_cl) {
  
  dims <- paste0(Key(object = object[[reduction]]), dims)
  data <- FetchData(object = object, vars = c(dims, "ident", features),  slot = "data")
  
  if (ncol(data) > 4) {
    data2 <- data.frame(data[,1:3], rowSums(data[, 4:ncol(data)]))
  } else {
    data2 <- data
  }
  
  l <- data2[[4]][which(data$ident %in% .retain_cl)]
  mi <- min(l)
  ma <- max(l)
  ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
  ColorLevels <- seq(mi, ma, length = length(ColorRamp))
  v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
  
  kk <- bind_cols(data.frame('l'=l), data[, dims][which(data$ident %in% .retain_cl),]) %>% arrange(l)
  colnames(kk)[2:3] <- c("UMAP_1", "UMAP_2")
  
  if(logsc) {
    plot <- ggplot(kk, aes(UMAP_1, UMAP_2, color = log(l+0.1))) +
      geom_point(size = point_size, pch = 19) +
      scale_color_gradientn('', colors = ColorRamp) +
      theme_void() +
      labs(title = paste(features, collapse = ',')) 
    return(plot)
  }
  else {
    plot <- ggplot(kk, aes(UMAP_1, UMAP_2, color = l)) +
      geom_point(size = point_size, pch = 19) +
      scale_color_gradientn('', colors = ColorRamp) +
      theme_void() +
      labs(title = paste(features, collapse = ','))
    return(plot)
  }
  
}

#Marimekko plot without stats
mosaicGG2 <- function(data, X, FILL, colors = colors_many, rect_col = 'white', line_width = 0.25) {
  require(dplyr)
  require(reshape2)
  #require(ggthemes)
  # Proportions in raw data
  DF <- as.data.frame.matrix(table(data[[X]], data[[FILL]]))
  DF$groupSum <- rowSums(DF)
  DF$xmax <- cumsum(DF$groupSum)
  DF$xmin <- DF$xmax - DF$groupSum
  DF$X <- row.names(DF)
  DF$groupSum <- NULL
  DF_melted <- melt(DF, id = c("X", "xmin", "xmax"), variable.name = "FILL")
  DF_melted <- DF_melted %>%
    group_by(X) %>%
    mutate(ymax = cumsum(value/sum(value)),
           ymin = ymax - value/sum(value))
  
  # Chi-sq test
  results <- chisq.test(table(data[[FILL]], data[[X]])) # fill and then x
  resid <- reshape2::melt(results$residuals)
  names(resid) <- c("FILL", "X", "residual")
  
  # Merge data
  DF_all <- merge(DF_melted, resid)
  
  # Positions for labels
  DF_all$xposn <- DF_all$xmin + (DF_all$xmax - DF_all$xmin)/2
  index <- DF_all$xmax == max(DF_all$xmax)
  #DF_all$yposn <- DF_all$ymin[index] + (DF_all$ymax[index] - DF_all$ymin[index])/2
  yposn = 0
  # Plot
  g <- ggplot(DF_all, aes(ymin = ymin,  ymax = ymax, xmin = xmin,
                          xmax = xmax, fill = FILL)) +
    geom_rect(col = rect_col, lwd = line_width) +
    geom_text(aes(x = xposn, label = X),
              y = 1, size = 3, angle = 90, hjust = 1, show.legend = FALSE,check_overlap = T) +
    geom_text(aes(x = max(xmax),  y = yposn, label = FILL),
              size = 3, hjust = 1, show.legend = FALSE,check_overlap = T) +
    scale_fill_manual(FILL, values = colors) +
    scale_x_continuous(X, expand = c(0,0)) +
    scale_y_continuous("Proportion", expand = c(0,0)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(g)
}

hyper_test_n <- function(data = df, var1 = "Cluster", var2 = "Region") {
  require(tidyverse) 
  require(broom)
  
  .df <- data.frame()
  for (i in unique(data[[var2]])) {
    data2 <- data
    data2[[var2]] <- factor(ifelse(data2[[var2]] == i, i, paste0("non_",i)), levels = c(i, paste0("non_",i)))
    clusters <- as_tibble(table(data2[[var1]]), .name_repair = 'unique')
    colnames(clusters) <- c(var1, 'cluster_size')
    vars <- as_tibble(table(data2[[var1]], data2[[var2]]), .name_repair = 'unique')
    colnames(vars) <- c(var1, var2, "freq_var2")
    vars_wide <- spread(vars, var2, freq_var2)
    
    vars_df <- vars_wide %>%
      left_join(clusters)
    
    
    #hypergeometric test
    #option a
    test_df<- data.frame(q=vars_df[,i], 
                         m=sum(vars_df[,i]), 
                         n=sum(vars_df[,paste0("non_",i)]),
                         k=vars_df[,4])
    
    colnames(test_df)[1] <- "q"
    
    p_hyper <- apply(test_df, MARGIN = 1, function(x) 1-phyper(max(0,x[[1]]-1), x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
    
    test_df$p_hyper <- p_hyper
    test_df$Cluster <- vars_df[[`var1`]]
    test_df$enrichment_var <- i
    .df <- .df %>%
      bind_rows(test_df[,c("q","m","n","cluster_size","p_hyper","Cluster","enrichment_var")])
  }
  
  
  .df$padj <- p.adjust(.df$p_hyper, method="BH")
  .df$Significance <- ifelse(.df$padj<0.05 & .df$padj>0.01, '*',
                             ifelse(.df$padj<0.01 & .df$padj>0.001, '**',
                                    ifelse(.df$padj<0.001, '***','n.s.')))
  
  return(.df)
}

go_term_analysis_seurat <- function(.df = df, ontogeny = "BP", .sc = all, organism = 'org.Hs.eg.db') {
  require(clusterProfiler)
  require(organism,character.only = TRUE)
  require(tidyverse)
  require(viridis)
  require(pheatmap)
  
  back_genes <- rownames(.sc@assays$RNA@counts)[which(rowMeans(as.matrix(.sc@assays$RNA@counts)) > 0)]
  back_genes <- gsub('_.*', '', back_genes)
  
  background <- bitr(back_genes, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                     OrgDb = organism)
  background <- background[!duplicated(background$ENTREZID),]
  
  #define empty data frame to collect data
  enrich_up <- data.frame(matrix(ncol = 10))
  colnames(enrich_up) <- c('ID','Description', 'GeneRatio', 'BgRatio' ,'pvalue', 'p.adjust', 'qvalue', 'geneID','Count' , 'Cluster')
  
  for (i in unique(.df$cluster))  {
    
    tryCatch({
      gene <- .df$gene[.df$cluster == i]
      gene.df <- bitr(gene, fromType = "SYMBOL",
                      toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                      OrgDb = organism)
      
      
      ggo <- groupGO(gene     = gene.df[,3],
                     OrgDb    = organism,
                     ont      = ontogeny,
                     level    = 3,
                     readable = TRUE)
      
      
      ego <- enrichGO(gene          = gene.df[,3],
                      universe      = background[,3],
                      OrgDb         = organism,
                      minGSSize     = 1,
                      ont           = ontogeny,
                      pool          = TRUE,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
      
      
      ego_simpl <- my_simplify(res = ego,
                               semData = godata(ont = ego@ontology))
      ego_simpl <- ego_simpl[!duplicated(ego_simpl$geneID),]
      
      ego_simpl$Cluster <- rep(i, nrow(ego_simpl))
      
      enrich_up <- rbind(enrich_up, ego_simpl)
      
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  return(na.omit(enrich_up))
  
}

my_simplify <- function(res=ego,
                        measure = 'Wang',
                        semData = godata(ont = ego@ontology),
                        by="pvalue",
                        cutoff=0.7,
                        select_fun=min) {
  
  require(GOSemSim)
  
  sim <- mgoSim(res$ID, res$ID,
                semData = semData,
                measure=measure,
                combine=NULL)
  
  ## to satisfy codetools for calling gather
  go1 <- go2 <- similarity <- NULL
  
  sim.df <- as.data.frame(sim)
  sim.df$go1 <- row.names(sim.df)
  sim.df <- gather(sim.df, go2, similarity, -go1)
  
  sim.df <- sim.df[!is.na(sim.df$similarity),]
  
  ## feature 'by' is attached to 'go1'
  sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")
  sim.df$go2 <- as.character(sim.df$go2)
  
  ID <- res$ID
  
  GO_to_remove <- character()
  for (i in seq_along(ID)) {
    ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
    ## if length(ii) == 1, then go1 == go2
    if (length(ii) < 2)
      next
    
    sim_subset <- sim.df[ii,]
    
    jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))
    
    ## sim.df <- sim.df[-ii[-jj]]
    GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
  }
  
  enrich_go <- res[!res$ID %in% GO_to_remove, ]
}

plot_index <- function(gene, .index=index_all, point_size=5, log=T) {
  l <- .index[[gene]] + 0.1
  mi <- min(l, na.rm = T)
  ma <- max(l, na.rm = T)
  ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
  ColorLevels <- seq(mi, ma, length = length(ColorRamp))
  v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
  
  kk <- bind_cols(data.frame('l'=l), .index[,c('UMAP_1', 'UMAP_2')]) %>% arrange(l)
  
  if (log){ 
    plot <- ggplot(na.omit(kk), aes(UMAP_1, UMAP_2, color = log(l))) +
      geom_point(size = point_size, pch = 19, stroke=0.25) +
      scale_color_gradientn('', colors = ColorRamp) +
      theme_void() +
      labs(title = paste(gene, collapse = ','))}
  
  if (!log){ 
    plot <- ggplot(na.omit(kk), aes(UMAP_1, UMAP_2, color = l)) +
      geom_point(size = point_size, pch = 19, stroke=0.25) +
      scale_color_gradientn('', colors = ColorRamp) +
      theme_void() +
      labs(title = paste(gene, collapse = ','))}
  return(plot)
}

plot_continuous <- function(param, .index=.ent, point_size=4) {
  l <- .index[[param]] + 0.1
  .l <- (l-min(l))/(max(l)-min(l))
  #mi <- min(.l, na.rm = T)
  #ma <- max(.l, na.rm = T)
  ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
  ColorLevels <- seq(0, 1, length = length(ColorRamp))
  
  kk <- bind_cols(data.frame('l'=l), .index[,c('UMAP_1', 'UMAP_2')]) %>% arrange(l)
  
  plot <- ggplot(na.omit(kk), aes(UMAP_1, UMAP_2, color = l)) +
    geom_point(size = point_size, pch = 19, stroke=0.25) +
    scale_color_gradientn('', colors = ColorRamp) +
    theme_void() +
    labs(title = paste(param, collapse = ','))
  return(plot)
}

plotheatmap2 <- function (x, xpart = NULL, xcol = NULL, xlab = TRUE, xgrid = FALSE, 
                          ypart = NULL, ycol = NULL, ylab = TRUE, ygrid = FALSE, cex = 1) 
{
  mi <- min(x, na.rm = TRUE)
  ma <- max(x, na.rm = TRUE)
  pardefault <- par()
  layout(matrix(data = c(1, 2), nrow = 1, ncol = 2), widths = c(5, 
                                                                1), heights = c(5, 1))
  ColorRamp <- rev(colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100))
  ColorLevels <- seq(mi, ma, length = length(ColorRamp))
  if (mi == ma) {
    ColorLevels <- seq(0.99 * mi, 1.01 * ma, length = length(ColorRamp))
  }
  par(mar = c(3, 5, 2.5, 2))
  image(t(as.matrix(1)), col = ColorRamp, axes = FALSE, ylim = c(-0.02, 1))
  box()
  set.seed(20)
  if (!is.null(xpart)) {
    tmp <- c()
    width <- (1/length(xpart))/2
    k <- (0:(length(xpart) - 1)/(length(xpart) - 1))
    rect(k - width, rep(-0.02, length(xpart)), k + width, 
         rep(-0.005, length(xpart)), col = xcol[xpart], border = NA)
    for (u in unique(xpart)) {
      ol <- (0:(length(xpart) - 1)/(length(xpart) - 1))[xpart == 
                                                          u]
      tmp <- append(tmp, mean(ol))
      delta <- 0.5/(length(xpart) - 1)
      if (xgrid & max(ol) < 1) 
        abline(v = max(ol) + delta, col = "grey", lty = 2)
    }
    if (xlab) 
      axis(1, at = tmp, labels = unique(xpart), cex.axis = cex)
  }
  set.seed(20)
  if (!is.null(ypart)) {
    tmp <- c()
    for (u in unique(ypart)) {
      ol <- (0:(length(ypart) - 1)/(length(ypart) - 1))[ypart == 
                                                          u]
      if (!is.null(ycol)) 
        points(rep(0, length(ol)), ol, col = ycol[u + 
                                                    1], pch = 15, cex = 0.75)
      tmp <- append(tmp, mean(ol))
      delta <- 0.5/(length(ypart) - 1)
      if (ygrid & max(ol) < 1) 
        abline(a = max(ol) + delta, b = 0, col = "grey", 
               lty = 2)
    }
    if (ylab) 
      axis(2, at = tmp, labels = unique(ypart), cex.axis = cex, 
           las = 1)
  }
  par(mar = c(10, 2, 2.5, 2))
  image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
                               nrow = 1), col = ColorRamp, xlab = "", ylab = "", xaxt = "n")
  layout(1)
  par(mar = pardefault$mar)
}

my_plot_cpdb2 <- function (cell_type1, cell_type2, scdata, idents, means, pvals, 
                           deconvoluted, p.adjust.method = NULL, keep_significant_only = TRUE, 
                           split.by = NULL, standard_scale = TRUE, separator = NULL, 
                           gene_symbol_mapping = NULL, frac = 0.2, remove_self = TRUE, 
                           desiredInteractions = NULL, interaction_grouping = NULL, 
                           edge_group_colors = NULL, node_group_colors = NULL, ...) 
{
  if (length(separator) > 0) {
    sep = separator
  }
  else {
    sep = ">@<"
  }
  cpdb_int = plot_cpdb(cell_type1 = cell_type1, cell_type2 = cell_type2, 
                       scdata = scdata, idents = idents, gene.family = "costimulatory",#split.by = split.by, 
                       means = means, pvals = pvalues, keep_significant_only = keep_significant_only, 
                       standard_scale = standard_scale, ...)
  lr_interactions <- cpdb_int$data
  subset_clusters <- unique(unlist(lapply(as.list(lr_interactions$group), 
                                          strsplit, sep)))
  sce_subset <- scdata[, colData(scdata)[, idents] %in% subset_clusters]
  interactions <- means[, c("interacting_pair", "gene_a", "gene_b", 
                            "partner_a", "partner_b")]
  interactions$converted <- gsub("-", " ", interactions$interacting_pair)
  interactions$converted <- gsub("_", "-", interactions$converted)
  interactions_subset <- interactions[interactions$converted %in% 
                                        lr_interactions$Var1, ]
  if (any(grepl("ENSG", interactions_subset$gene_a))) {
    tm0 <- do.call(c, lapply(as.list(interactions_subset$interacting_pair), 
                             strsplit, "_"))
    tm0 = data.frame(t(matrix(unlist(tm0), 2, length(unlist(tm0))/2)))
    colnames(tm0) <- c("id_a", "id_b")
    interactions_subset <- cbind(interactions_subset, tm0)
    dictionary_a <- interactions_subset$id_a
    names(dictionary_a) <- interactions_subset$gene_a
    dictionary_b <- interactions_subset$id_b
    names(dictionary_b) <- interactions_subset$gene_b
    dictionary <- c(as.list(dictionary_a), as.list(dictionary_b))
    if (length(which(names(dictionary) == "")) > 0) {
      dictionary <- dictionary[-which(names(dictionary) == 
                                        "")]
    }
  }
  if (!is.null(interaction_grouping)) {
    if ((class(interaction_grouping) == "data.frame")) {
      interactions_subset$group = interaction_grouping[, 
                                                       2][match(interactions_subset$interacting_pair, 
                                                                interaction_grouping[, 1])]
    }
  }
  geneid = unique(c(interactions_subset$gene_a, interactions_subset$gene_b))
  rmg = which(geneid == "")
  if (length(rmg) > 0) {
    geneid = geneid[-which(geneid == "")]
  }
  sce_subset_tmp <- sce_subset[geneid %in% rownames(scdata), ]
  if (class(scdata) %in% c("SingleCellExperiment", "SummarizedExperiment")) {
    requireNamespace("SummarizedExperiment")
    requireNamespace("SingleCellExperiment")
    meta <- as.data.frame(colData(sce_subset_tmp))
  }
  else if (class(scdata) == "Seurat") {
    stop("sorry not supported yet")
  }
  sce_list <- list()
  sce_list_alt <- list()
  for (x in unique(meta[, split.by])) {
    sce_list[[x]] <- list()
    sce_list_alt[[x]] <- list()
  }
  for (n in names(sce_list)) {
    for (x in unique(meta[, idents])) {
      sce_list[[n]][[x]] <- sce_subset_tmp[, meta[, idents] == 
                                             x & meta[, split.by] == n]
      sce_list_alt[[n]][[x]] <- sce_subset[, meta[, idents] == 
                                             x & meta[, split.by] == n]
    }
  }
  cellTypeMeans <- function(x) {
    cm <- Matrix::rowMeans(counts(x))
    return(cm)
  }
  cellTypeFraction <- function(x) {
    cm <- Matrix::rowMeans(counts(x) > 0)
    return(cm)
  }
  sce_list2 <- lapply(sce_list, function(y) {
    z <- lapply(y, cellTypeMeans)
    return(z)
  })
  sce_list3 <- lapply(sce_list, function(y) {
    z <- lapply(y, cellTypeFraction)
    return(z)
  })
  sce_list2 <- lapply(sce_list2, function(x) do.call(cbind, 
                                                     x))
  sce_list3 <- lapply(sce_list3, function(x) do.call(cbind, 
                                                     x))
  for (n in names(sce_list2)) {
    colnames(sce_list2[[n]]) <- paste0(paste0(n, "_"), colnames(sce_list2[[n]]))
    colnames(sce_list3[[n]]) <- paste0(paste0(n, "_"), colnames(sce_list3[[n]]))
  }
  sce_list2 <- do.call(cbind, sce_list2)
  sce_list3 <- do.call(cbind, sce_list3)
  humanreadablename = c()
  for (i in row.names(sce_list2)) {
    humanreadablename = c(humanreadablename, as.character(unlist(dictionary[i])))
  }
  rownames(sce_list2) <- humanreadablename
  rownames(sce_list3) <- humanreadablename
  findComplex <- function(interaction) {
    idxa <- which(interaction$gene_a == "")
    idxb <- which(interaction$gene_b == "")
    complexa <- gsub("complex:", "", interaction$partner_a[idxa])
    complexb <- gsub("complex:", "", interaction$partner_b[idxb])
    if (length(complexa) > 0) {
      if (length(complexb) > 0) {
        res <- c(complexa, complex_b)
      }
      else {
        res <- complexa
      }
    }
    else if (length(complexb) > 0) {
      res <- complexb
    }
    else {
      res <- NULL
    }
    return(res)
  }
  decon_subset <- deconvoluted[deconvoluted$complex_name %in% 
                                 findComplex(interactions_subset), ]
  if (nrow(decon_subset) > 0) {
    decon_subset <- split(decon_subset, decon_subset$complex_name)
    decon_subset_expr <- lapply(decon_subset, function(x) {
      x <- x[, colnames(sce_list2)]
      x <- colMeans(x)
      return(x)
    })
    cellTypeFraction_complex <- function(sce_, genes, gene_symbol_mapping = NULL) {
      scex <- tryCatch(sce_[genes, ], error = function(e) {
        if (!is.null(gene_symbol_mapping)) {
          sce_[which(rowData(sce_)[, gene_symbol_mapping] %in% 
                       genes), ]
        }
        else {
          sce_[which(rowData(sce_)[, "index"] %in% genes), 
          ]
        }
      })
      cm <- mean(Matrix::rowMeans(counts(scex) > 0))
      return(cm)
    }
    decon_subset_fraction <- lapply(decon_subset, function(x) {
      x <- unique(x$gene_name)
      test <- lapply(sce_list_alt, function(y) {
        return(lapply(y, cellTypeFraction_complex, x, 
                      gene_symbol_mapping))
      })
      return(test)
    })
    decon_subset_fraction <- lapply(decon_subset_fraction, 
                                    function(x) {
                                      y <- lapply(x, function(z) do.call(cbind, z))
                                      for (i in 1:length(y)) {
                                        colnames(y[[i]]) <- paste0(names(y[i]), "_", 
                                                                   colnames(y[[i]]))
                                      }
                                      y <- do.call(cbind, y)
                                      return(y)
                                    })
    decon_subset_expr <- do.call(rbind, decon_subset_expr)
    decon_subset_fraction <- do.call(rbind, decon_subset_fraction)
    row.names(decon_subset_fraction) <- row.names(decon_subset_expr)
    expr_df <- rbind(sce_list2, decon_subset_expr)
    fraction_df <- rbind(sce_list3, decon_subset_fraction)
  }
  else {
    expr_df <- sce_list2
    fraction_df <- sce_list3
  }
  if (!is.null(desiredInteractions)) {
    if (class(desiredInteractions) == "list") {
      desiredInteractions_ <- c(desiredInteractions, lapply(desiredInteractions, 
                                                            rev))
      cell_type_grid <- as.data.frame(do.call(rbind, desiredInteractions_))
    }
    else if ((class(desiredInteractions) == "data.frame")) {
      cell_type_grid <- desiredInteractions
    }
    cells_test = unique(unlist(desiredInteractions))
  }
  else {
    cells_test <- unique(droplevels(meta[, idents]))
    cell_type_grid <- expand.grid(cells_test, cells_test)
  }
  if (remove_self) {
    rm_idx <- which(cell_type_grid[, 1] == cell_type_grid[, 
                                                          2])
    if (length(rm_idx) > 0) {
      cell_type_grid <- cell_type_grid[-rm_idx, ]
    }
  }
  ligand <- interactions_subset$id_a
  receptor <- interactions_subset$id_b
  pair <- interactions_subset$interacting_pair
  converted_pair <- interactions_subset$converted
  producers <- as.character(cell_type_grid[, 1])
  receivers <- as.character(cell_type_grid[, 2])
  generateDf <- function(ligand, receptor, pair, converted_pair, 
                         producers, receivers, cell_type_means, cell_type_fractions, 
                         splitted = NULL) {
    if (!is.null(splitted)) {
      pp <- paste0(splitted, "_", producers)
      rc <- paste0(splitted, "_", receivers)
    }
    else {
      pp <- producers
      rc <- receivers
    }
    df_ <- data.frame(ligand = ligand, receptor = receptor, 
                      pair = pair, producer = pp, receiver = rc, producer_expression = as.numeric(cell_type_means[ligand, 
                                                                                                                  pp]), producer_fraction = as.numeric(cell_type_fractions[ligand, 
                                                                                                                                                                           pp]), receiver_expression = as.numeric(cell_type_means[receptor, 
                                                                                                                                                                                                                                  rc]), receiver_fraction = as.numeric(cell_type_fractions[receptor, 
                                                                                                                                                                                                                                                                                           rc]))
    df_$from = paste0(df_$producer, "_", df_$ligand)
    df_$to = paste0(df_$receiver, "_", df_$receptor)
    if (!is.null(splitted)) {
      df_$producer_ = df_$producer
      df_$receiver_ = df_$receiver
      df_$from = gsub(paste0(splitted, "_"), "", df_$from)
      df_$to = gsub(paste0(splitted, "_"), "", df_$to)
      df_$producer = gsub(paste0(splitted, "_"), "", df_$producer)
      df_$receiver = gsub(paste0(splitted, "_"), "", df_$receiver)
      df_$barcode = paste0(df_$producer_, "-", df_$receiver_, 
                           sep, converted_pair)
    }
    else {
      df_$barcode = paste0(df_$producer, "-", df_$receiver, 
                           sep, converted_pair)
    }
    return(df_)
  }
  dfx <- list()
  if (!is.null(split.by)) {
    for (i in unique(meta[, split.by])) {
      dfx[[i]] = generateDf(ligand, receptor, pair, converted_pair, 
                            producers, receivers, expr_df, fraction_df, i)
    }
  }
  else {
    dfx[[1]] = generateDf(ligand, receptor, pair, converted_pair, 
                          producers, receivers, expr_df, fraction_df)
  }
  df0 <- lapply(dfx, function(x) x[x$producer_fraction > frac & 
                                     x$receiver_fraction > frac, ])
  constructGraph <- function(el, el0, unique_id, interactions_df, 
                             plot_cpdb_out, edge_group = FALSE, edge_group_colors = NULL, 
                             node_group_colors = NULL) {
    require(igraph)
    celltypes <- unique(c(as.character(el$producer), as.character(el$receiver)))
    el1 <- data.frame(from = "root", to = celltypes, barcode_1 = NA, 
                      barcode_2 = NA, barcode_3 = NA)
    el2 <- data.frame(from = celltypes, to = paste0(celltypes, 
                                                    "_", "ligand"), barcode_1 = NA, barcode_2 = NA, barcode_3 = NA)
    el3 <- data.frame(from = celltypes, to = paste0(celltypes, 
                                                    "_", "receptor"), barcode_1 = NA, barcode_2 = NA, 
                      barcode_3 = NA)
    el4 <- do.call(rbind, lapply(celltypes, function(x) {
      cell_ligands <- grep(x, el$from, value = TRUE)
      cell_ligands_idx <- grep(x, el$from)
      if (length(cell_ligands) > 0) {
        df <- data.frame(from = paste0(x, "_", "ligand"), 
                         to = cell_ligands, barcode_1 = el$barcode[cell_ligands_idx], 
                         barcode_2 = el$pair[cell_ligands_idx], barcode_3 = paste0(el$from[cell_ligands_idx], 
                                                                                   sep, el$to[cell_ligands_idx]))
      }
      else {
        df = NULL
      }
    }))
    el5 <- do.call(rbind, lapply(celltypes, function(x) {
      cell_ligands <- grep(x, el$to, value = TRUE)
      cell_ligands_idx <- grep(x, el$to)
      if (length(cell_ligands) > 0) {
        df <- data.frame(from = paste0(x, "_", "receptor"), 
                         to = cell_ligands, barcode_1 = el$barcode[cell_ligands_idx], 
                         barcode_2 = el$pair[cell_ligands_idx], barcode_3 = paste0(el$from[cell_ligands_idx], 
                                                                                   sep, el$to[cell_ligands_idx]))
      }
      else {
        df = NULL
      }
    }))
    gr_el <- do.call(rbind, list(el1, el2, el3, el4, el5))
    plot_cpdb_out$barcode <- paste0(plot_cpdb_out$Var2, sep, 
                                    plot_cpdb_out$Var1)
    mean_col <- grep("means$", colnames(plot_cpdb_out), value = TRUE)
    means <- plot_cpdb_out[match(gr_el$barcode_1, plot_cpdb_out$barcode), 
                           mean_col]
    pval_col <- grep("pvals", colnames(plot_cpdb_out), value = TRUE)
    pvals <- plot_cpdb_out[match(gr_el$barcode_1, plot_cpdb_out$barcode), 
                           pval_col]
    gr_el <- cbind(gr_el, means, pvals)
    if (edge_group) {
      groups <- interactions_df$group[match(gr_el$barcode_2, 
                                            interactions_df$interacting_pair)]
    }
    gr <- graph_from_edgelist(as.matrix(gr_el[, 1:2]))
    E(gr)$interaction_score <- means
    E(gr)$pvals <- pvals
    if (edge_group) {
      E(gr)$group <- groups
    }
    E(gr)$name <- gr_el$barcode_3
    V(gr)$type <- NA
    V(gr)$type[V(gr)$name %in% el4$to] <- "ligand"
    V(gr)$type[V(gr)$name %in% el5$to] <- "receptor"
    from = match(el0$from, V(gr)$name)
    to = match(el0$to, V(gr)$name)
    dat = data.frame(from = el0$from, to = el0$to)
    dat$barcode = paste0(dat$from, sep, dat$to)
    interaction_score = E(gr)$interaction_score[match(dat$barcode, 
                                                      gr_el$barcode_3)]
    pval = E(gr)$pvals[match(dat$barcode, gr_el$barcode_3)]
    pval[is.na(pval)] <- 1
    pval <- range01(-log10(pval))
    if (edge_group) {
      group = E(gr)$group[match(dat$barcode, gr_el$barcode_3)]
    }
    ligand_expr <- data.frame(cell_mol = el$from, expression = el$producer_expression, 
                              fraction = el$producer_fraction)
    recep_expr <- data.frame(cell_mol = el$to, expression = el$receiver_expression, 
                             fraction = el$receiver_fraction)
    expression <- rbind(ligand_expr, recep_expr)
    V(gr)$expression = NA
    V(gr)$expression[match(expression$cell_mol, V(gr)$name)] <- expression$expression
    V(gr)$fraction = 0
    V(gr)$fraction[match(expression$cell_mol, V(gr)$name)] <- expression$fraction
    V(gr)$celltype <- NA
    for (x in cells_test) {
      idx <- grepl(x, V(gr)$name)
      V(gr)$celltype[idx] <- x
    }
    V(gr)$name[!V(gr)$name %in% c(el0$from, el0$to)] <- NA
    for (x in unique_id) {
      V(gr)$name <- gsub(paste0(x, "_"), "", V(gr)$name)
    }
    library(ggraph)
    library(ggrepel)
    if (!is.null(edge_group_colors)) {
      edge_group_colors = edge_group_colors
    }
    else {
      nn = length(unique(E(gr)$group))
      edge_group_colors = kelvinny::gg_color_hue(nn)
    }
    if (!is.null(node_group_colors)) {
      node_group_colors = node_group_colors
    }
    else {
      nn = length(unique(meta[, idents]))
      node_group_colors = kelvinny::gg_color_hue(nn)
    }
    pl <- ggraph(gr, layout = "dendrogram", circular = TRUE) + 
      geom_conn_bundle(data = get_con(from = from, to = to, 
                                      group = group, `-log10(sig)` = pval, interaction_score = interaction_score), 
                       aes(colour = group, alpha = interaction_score, 
                           width = `-log10(sig)`), tension = 0.5) + scale_edge_width(range = c(1, 
                                                                                               3)) + scale_edge_alpha(range = c(0, 3), limits = c(0, 
                                                                                                                                                  1)) + scale_edge_color_manual(values = edge_group_colors) + 
      geom_node_point(pch = 19, aes(size = fraction, filter = leaf, 
                                    color = celltype, alpha = type)) + theme_void() + 
      coord_fixed() + scale_size_continuous(limits = c(0, 
                                                       1)) + scale_shape_manual(values = c(ligand = 19, 
                                                                                           receptor = 15)) + scale_color_manual(values = node_group_colors) + 
      geom_text_repel(aes(x = x, y = y, label = name), 
                      segment.square = TRUE, segment.inflect = TRUE, 
                      segment.size = 0.2, force = 0.5, size = 2, force_pull = 0) + 
      scale_alpha_manual(values = c(ligand = 0.5, receptor = 1)) + 
      small_legend(keysize = 0.5)
    return(pl)
  }
  gl <- list()
  if (!is.null(split.by)) {
    edge_group = TRUE
  }
  else {
    edge_group = FALSE
  }
  for (i in 1:length(dfx)) {
    gl[[i]] <- constructGraph(dfx[[i]], df0[[i]], cells_test, 
                              interactions_subset, lr_interactions, edge_group, 
                              edge_group_colors, node_group_colors)
  }
  if (length(gl) > 1) {
    return(gl)
  }
  else {
    return(gl[[1]])
  }
}

computeNetSimilarity <- function(object, slot.name = "netP", type = c("functional","structural"), k = NULL, thresh = NULL) {
  type <- match.arg(type)
  prob = methods::slot(object, slot.name)$prob
  if (is.null(k)) {
    if (dim(prob)[3] <= 25) {
      k <- ceiling(sqrt(dim(prob)[3]))
    } else {
      k <- ceiling(sqrt(dim(prob)[3])) + 1
    }
  }
  if (!is.null(thresh)) {
    prob[prob < quantile(c(prob[prob != 0]), thresh)] <- 0
  }
  if (type == "functional") {
    # compute the functional similarity
    D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
    S2 <- D_signalings; S3 <- D_signalings;
    for (i in 1:(dim(prob)[3]-1)) {
      for (j in (i+1):dim(prob)[3]) {
        Gi <- (prob[ , ,i] > 0)*1
        Gj <- (prob[ , ,j] > 0)*1
        S3[i,j] <- sum(Gi * Gj)/sum(Gi+Gj-Gi*Gj,na.rm=TRUE)
      }
    }
    # define the similarity matrix
    S3[is.na(S3)] <- 0; S3 <- S3 + t(S3); diag(S3) <- 1
    # S_signalings <- S1 *S2
    S_signalings <- S3
  } else if (type == "structural") {
    # compute the structure distance
    D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
    for (i in 1:(dim(prob)[3]-1)) {
      for (j in (i+1):dim(prob)[3]) {
        Gi <- (prob[ , ,i] > 0)*1
        Gj <- (prob[ , ,j] > 0)*1
        D_signalings[i,j] <- computeNetD_structure(Gi,Gj)
      }
    }
    # define the structure similarity matrix
    D_signalings[is.infinite(D_signalings)] <- 0
    D_signalings[is.na(D_signalings)] <- 0
    D_signalings <- D_signalings + t(D_signalings)
    S_signalings <- 1-D_signalings
  }
  # smooth the similarity matrix using SNN
  SNN <- buildSNN(S_signalings, k = k, prune.SNN = 1/15)
  Similarity <- as.matrix(S_signalings*SNN)
  rownames(Similarity) <- dimnames(prob)[[3]]
  colnames(Similarity) <- dimnames(prob)[[3]]
  comparison <- "single"
  comparison.name <- paste(comparison, collapse = "-")
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$matrix)) {
    methods::slot(object, slot.name)$similarity[[type]]$matrix <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]] <- Similarity
  return(object)
}
#' Compute signaling network similarity for any pair of datasets
#'
#' @param object A merged CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param comparison a numerical vector giving the datasets for comparison
#' @param k the number of nearest neighbors
#' @param thresh the fraction (0 to 0.25) of interactions to be trimmed before computing network similarity
#' @importFrom methods slot
#'
#' @return
#' @export
#'
computeNetSimilarityPairwise <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, k = NULL, thresh = NULL) {
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("Compute signaling network similarity for datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")
  net <- list()
  signalingAll <- c()
  object.net.nameAll <- c()
  # 1:length(setdiff(names(methods::slot(object, slot.name)), "similarity"))
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, slot.name)[[comparison[i]]]
    object.net.name <- names(methods::slot(object, slot.name))[comparison[i]]
    object.net.nameAll <- c(object.net.nameAll, object.net.name)
    net[[i]] = object.net$prob
    signalingAll <- c(signalingAll, paste0(dimnames(net[[i]])[[3]], "--", object.net.name))
    # signalingAll <- c(signalingAll, dimnames(net[[i]])[[3]])
  }
  names(net) <- object.net.nameAll
  net.dim <- sapply(net, dim)[3,]
  nnet <- sum(net.dim)
  position <- cumsum(net.dim); position <- c(0,position)
  if (is.null(k)) {
    if (nnet <= 25) {
      k <- ceiling(sqrt(nnet))
    } else {
      k <- ceiling(sqrt(nnet)) + 1
    }
  }
  if (!is.null(thresh)) {
    for (i in 1:length(net)) {
      neti <- net[[i]]
      neti[neti < quantile(c(neti[neti != 0]), thresh)] <- 0
      net[[i]] <- neti
    }
  }
  if (type == "functional") {
    # compute the functional similarity
    S3 <- matrix(0, nrow = nnet, ncol = nnet)
    for (i in 1:nnet) {
      for (j in 1:nnet) {
        idx.i <- which(position - i >= 0)[1]
        idx.j <- which(position - j >= 0)[1]
        net.i <- net[[idx.i-1]]
        net.j <- net[[idx.j-1]]
        Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
        Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
        S3[i,j] <- sum(Gi * Gj)/sum(Gi+Gj-Gi*Gj,na.rm=TRUE)
      }
    }
    # define the similarity matrix
    S3[is.na(S3)] <- 0;  diag(S3) <- 1
    S_signalings <- S3
  } else if (type == "structural") {
    # compute the structure distance
    D_signalings <- matrix(0, nrow = nnet, ncol = nnet)
    for (i in 1:nnet) {
      for (j in 1:nnet) {
        idx.i <- which(position - i >= 0)[1]
        idx.j <- which(position - j >= 0)[1]
        net.i <- net[[idx.i-1]]
        net.j <- net[[idx.j-1]]
        Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
        Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
        D_signalings[i,j] <- computeNetD_structure(Gi,Gj)
      }
    }
    # define the structure similarity matrix
    D_signalings[is.infinite(D_signalings)] <- 0
    D_signalings[is.na(D_signalings)] <- 0
    S_signalings <- 1-D_signalings
  }
  # smooth the similarity matrix using SNN
  SNN <- buildSNN(S_signalings, k = k, prune.SNN = 1/15)
  Similarity <- as.matrix(S_signalings*SNN)
  rownames(Similarity) <- signalingAll
  colnames(Similarity) <- rownames(Similarity)
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$matrix)) {
    methods::slot(object, slot.name)$similarity[[type]]$matrix <- NULL
  }
  # methods::slot(object, slot.name)$similarity[[type]]$matrix <- Similarity
  methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]] <- Similarity
  return(object)
}
#' Manifold learning of the signaling networks based on their similarity
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param comparison a numerical vector giving the datasets for comparison. No need to define for a single dataset. Default are all datasets when object is a merged object
#' @param k the number of nearest neighbors in running umap
#' @param pathway.remove a range of the number of patterns
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' 
#' 
#' 
#' 
#' 
netEmbedding <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, pathway.remove = NULL, k = NULL) {
  if (object@options$mode == "single") {
    comparison <- "single"
    cat("Manifold learning of the signaling networks for a single dataset", '\n')
  } else if (object@options$mode == "merged") {
    if (is.null(comparison)) {
      comparison <- 1:length(unique(object@meta$datasets))
    }
    cat("Manifold learning of the signaling networks for datasets", as.character(comparison), '\n')
  }
  comparison.name <- paste(comparison, collapse = "-")
  Similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
  if (is.null(pathway.remove)) {
    pathway.remove <- rownames(Similarity)[which(colSums(Similarity) == 1)]
  }
  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(rownames(Similarity) %in% pathway.remove)
    Similarity <- Similarity[-pathway.remove.idx, -pathway.remove.idx]
  }
  if (is.null(k)) {
    k <- ceiling(sqrt(dim(Similarity)[1])) + 1
  }
  options(warn = -1)
  # dimension reduction
  Y <- runUMAP(Similarity, min.dist = 0.3, n.neighbors = k)
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$dr)) {
    methods::slot(object, slot.name)$similarity[[type]]$dr <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]] <- Y
  return(object)
}
#' Classification learning of the signaling networks
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param comparison a numerical vector giving the datasets for comparison. No need to define for a single dataset. Default are all datasets when object is a merged object
#' @param k the number of signaling groups when running kmeans
#' @param methods the methods for clustering: "kmeans" or "spectral"
#' @param do.plot whether showing the eigenspectrum for inferring number of clusters; Default will save the plot
#' @param fig.id add a unique figure id when saving the plot
#' @param do.parallel whether doing parallel when inferring the number of signaling groups when running kmeans
#' @param nCores number of workers when doing parallel
#' @param k.eigen the number of eigenvalues used when doing spectral clustering
#' @importFrom methods slot
#' @importFrom future nbrOfWorkers plan
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @return
#' @export
#'
#' @examples

netClustering <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, k = NULL, methods = "kmeans", do.plot = TRUE, fig.id = NULL, do.parallel = TRUE, nCores = 4, k.eigen = NULL) {
  type <- match.arg(type)
  if (object@options$mode == "single") {
    comparison <- "single"
    cat("Classification learning of the signaling networks for a single dataset", '\n')
  } else if (object@options$mode == "merged") {
    if (is.null(comparison)) {
      comparison <- 1:length(unique(object@meta$datasets))
    }
    cat("Classification learning of the signaling networks for datasets", as.character(comparison), '\n')
  }
  comparison.name <- paste(comparison, collapse = "-")
  
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  Y[is.na(Y)] <- 0
  data.use <- Y
  if (methods == "kmeans") {
    if (!is.null(k)) {
      clusters = kmeans(data.use,k,nstart=10)$cluster
    } else {
      N <- nrow(data.use)
      kRange <- seq(2,min(N-1, 10),by = 1)
      if (do.parallel) {
        future::plan("multiprocess", workers = nCores)
        options(future.globals.maxSize = 1000 * 1024^2)
      }
      my.sapply <- ifelse(
        test = future::nbrOfWorkers() == 1,
        yes = pbapply::pbsapply,
        no = future.apply::future_sapply
      )
      results = my.sapply(
        X = 1:length(kRange),
        FUN = function(x) {
          idents <- kmeans(data.use,kRange[x],nstart=10)$cluster
          clusIndex <- idents
          #adjMat0 <- as.numeric(outer(clusIndex, clusIndex, FUN = "==")) - outer(1:N, 1:N, "==")
          adjMat0 <- Matrix::Matrix(as.numeric(outer(clusIndex, clusIndex, FUN = "==")), nrow = N, ncol = N)
          return(list(adjMat = adjMat0, ncluster = length(unique(idents))))
        },
        simplify = FALSE
      )
      adjMat <- lapply(results, "[[", 1)
      CM <- Reduce('+', adjMat)/length(kRange)
      res <- computeEigengap(as.matrix(CM))
      numCluster <- res$upper_bound
      clusters = kmeans(data.use,numCluster,nstart=10)$cluster
      if (do.plot) {
        gg <- res$gg.obj
        ggsave(filename= paste0("estimationNumCluster_",fig.id,"_",type,"_dataset_",comparison.name,".pdf"), plot=gg, width = 3.5, height = 3, units = 'in', dpi = 300)
      }
    }
  } else if (methods == "spectral") {
    A <- as.matrix(data.use)
    D <- apply(A, 1, sum)
    L <- diag(D)-A                       # unnormalized version
    L <- diag(D^-0.5)%*%L%*% diag(D^-0.5) # normalized version
    evL <- eigen(L,symmetric=TRUE)  # evL$values is decreasing sorted when symmetric=TRUE
    # pick the first k first k eigenvectors (corresponding k smallest) as data points in spectral space
    plot(rev(evL$values)[1:30])
    Z <- evL$vectors[,(ncol(evL$vectors)-k.eigen+1):ncol(evL$vectors)]
    clusters = kmeans(Z,k,nstart=20)$cluster
  }
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$group)) {
    methods::slot(object, slot.name)$similarity[[type]]$group <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]] <- clusters
  return(object)
}
