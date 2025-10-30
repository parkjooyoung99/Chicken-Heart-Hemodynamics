#!/usr/bin/env Rscript
# /programs/R-4.2.3/bin/R

setwd('/workdir/jp2626/chickenheart/jy/compareRL/Ventricle_sub/napari/noery')
set.seed(1234)

library(dplyr)
library(tidyr)
library(fgsea)
library(msigdbr)
library(simplifyEnrichment)
library(rrvgo)
library(circlize)
library(GO.db)
library("org.Gg.eg.db")
library(ggplot2)
library(stringr)

files = list.files(path = '/workdir/jp2626/chickenheart/jy/compareRL/Ventricle_sub/napari/noery/',pattern = '\\LV_new.csv$')


if (!dir.exists("GOBP_ventricle")) {
  dir.create("GOBP_ventricle")
}

msigdbr_df <- msigdbr(species = "Gallus gallus", collection = "C5", subcollection = 'GO:BP')
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
pathwaytoterm = msigdbr_df %>% dplyr::select('gs_name','gs_exact_source') %>% distinct()
enrichmate <- function(goid, plot.title, font.size = 15){
  notavail = c()
  ## 1. simplifyenrichmate
  #### 1.1. Calculate GO similarity
  sim_mat <- GO_similarity(goid)
  
  #### 1.2. Cluster GO terms base on the GO similarity matrix
  df <- simplifyGO(sim_mat, plot = F)
  df <- df %>% arrange(cluster)
  
  
  ## 2. Revigo - Define representatives
  representatives_per_cluster <- data.frame(t(rep(NA, 2)))
  representatives_per_cluster <- representatives_per_cluster[-1, ]
  for (c in unique(df$cluster)){
    tryCatch({
      if (length(df[df$cluster==c, 'id'])==1){
        tmp_id <- df[df$cluster==c, 'id']
        tmp_res <- data.frame(cluster = c,
                              parentTerm = na.omit(go_id[go_id$GOID==tmp_id, 'TERM']))
        representatives_per_cluster <- as.data.frame(rbind(representatives_per_cluster, tmp_res))
      }
      else{
        simMatrix_tmp <- GO_similarity(df[df$cluster==c, 'id'])
        
        reducedTerms_tmp <- rrvgo::reduceSimMatrix(simMatrix_tmp,
                                                   scores = NULL,
                                                   threshold=0.7,
                                                   orgdb="org.Gg.eg.db")
        
        reducedTerms_tmp[, c('cluster', 'term', 'parentTerm')] %>% arrange(cluster)
        tmp_res <- data.frame(cluster = rep(c, length(unique(reducedTerms_tmp$parentTerm))),
                              parentTerm = unique(reducedTerms_tmp$parentTerm))
        representatives_per_cluster <- as.data.frame(rbind(representatives_per_cluster, tmp_res))
      }
    }, error=function(e){print(c) })
  }
  representatives_per_cluster <- unique(representatives_per_cluster)
  representatives_per_cluster <- na.omit(representatives_per_cluster)
  
  ## 3. Heatmap
  #### 3.1. Decorate heatmap
  require(ComplexHeatmap)
  require(colorRamp2)
  require(colorspace)
  
  #### 3.1.1. Construct the text box
  representatives_per_cluster_list <- list()
  for (c in unique(representatives_per_cluster$cluster)){
    print(c)
    representatives_per_cluster_list[[c]] <- representatives_per_cluster[representatives_per_cluster$cluster==c, 'parentTerm']
  }
  
  
  representatives_per_cluster_list <- lapply(unique(representatives_per_cluster$cluster), function(x){
    df = data.frame(text = representatives_per_cluster_list[[x]])
    df$col <- rep(darken(rand_color(1)), nrow(df))
    df$fontsize <- font.size
    df
  })
  
  df = df[df$cluster %in% unique(representatives_per_cluster$cluster),]
  split <- df$cluster
  names(representatives_per_cluster_list) = unique(split)
  
  #### 3.1.2. Define split on a heatmap
  simMatrix_ord <- sim_mat[df$id, df$id] #### Sim. Matirx ordered
  split <- df$cluster
  names(representatives_per_cluster_list) = unique(split)
  
  #### 3.1.3. Color of heatmap
  col = c("white", 'red')
  col_fun = circlize::colorRamp2(seq(0, quantile(sim_mat, 0.95), length = length(col)), col, space = 'RGB')
  f2 = colorRamp2(seq(min(sim_mat), max(sim_mat), length = 3), c("white", "#EEEEEE", "red"), space = "RGB")
  
  if (min(sim_mat)==max(sim_mat)){
    ht <- Heatmap(simMatrix_ord, name = "semantic similarity", cluster_columns = F,cluster_rows = FALSE, 
                  column_title = plot.title,
                  show_row_names = T, row_split = split, column_split = split,
                  show_column_names = F, row_title = NULL,
                  row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"), border = TRUE,
                  right_annotation = rowAnnotation(textbox = anno_textbox(split, representatives_per_cluster_list,
                                                                          by='anno_block', add_new_line = T, word_wrap=F,
                                                                          background_gp = gpar(fill = 'white', col = 'gray'),
                                                                          round_corners = TRUE,
                                                                          line_space = unit(5, "mm"),
                                                                          max_width = unit(15, "cm")
                  )
                  )
    )  
  } else{
    ht <- Heatmap(simMatrix_ord, name = "semantic similarity", cluster_columns = F,cluster_rows = FALSE, col = col,
                  column_title = plot.title,
                  show_row_names = F, row_split = split, column_split = split,
                  show_column_names = F, row_title = NULL,
                  row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"), border = TRUE,
                  right_annotation = rowAnnotation(textbox = anno_textbox(split, representatives_per_cluster_list,
                                                                          by='anno_block', add_new_line = T, word_wrap=F,
                                                                          background_gp = gpar(fill = 'white', col = 'gray'),
                                                                          round_corners = TRUE,
                                                                          line_space = unit(5, "mm"),
                                                                          max_width = unit(15, "cm")
                  )
                  )
    )  
  }
  final_result <- list(cluster = df, representative = representatives_per_cluster, plot = ht)
  return(final_result)
}

for (f in files){
  print(f)
  
  prefix = f %>% str_remove_all('.csv')
  print(prefix)
  
  # Get DEG info
  res = read.csv(f, row.names=1)
  res_L = res %>% filter(group %>% str_detect('_L'))
  res_L_sig = res_L %>% filter(pvals_adj < 0.05)
  
  statdf <- res_L_sig %>%
    as.data.frame() %>%
    dplyr::select(names, scores, logfoldchanges, pvals_adj) %>%
    arrange(desc( scores ), desc(logfoldchanges), pvals_adj) %>%
    dplyr::select(names, scores)
  # Break the tie
  epsilon <- 1e-8
  statdf <- statdf %>%
    group_by(scores) %>%
    mutate(
      # only add epsilon if more than 1 row in this group
      stat_tiebreak = if(n() > 1) {
        scores + (row_number() - 1) * epsilon
      } else {
        scores
      }
    ) %>%
    ungroup()
  ranks <- tibble::deframe(statdf %>% dplyr::select(names, stat_tiebreak))
  
  # Map the ENGS to gene symbol
  for (n in 1:length(names(ranks))){
    i = names(ranks)[n]
    if (! i %in% msigdbr_df$gene_symbol){
      if (i %in% msigdbr_df$ensembl_gene){
        new_i = msigdbr_df[msigdbr_df$ensembl_gene == i,] %>% drop_na %>% .[,'gene_symbol'] %>% unique
        
        if(length(new_i) != 0){
          names(ranks)[n] = new_i[1]
        }
      }
    }
  }
  
  # Run fgsea
  fgseaRes <- fgsea(pathways=pathwaysH, stats=ranks, nproc = 1,nPermSimple = 100000)
  fgseaResTidy_p <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))

  fgseaResTidy_p_filt = fgseaResTidy_p[fgseaResTidy_p$padj < 0.05,]
  fgseaResTidy_p_filt_clean <- fgseaResTidy_p_filt %>% mutate(across(where(is.list), ~sapply(., function(x) paste(x, collapse = ";"))))
  write.csv(fgseaResTidy_p_filt_clean,paste0('./GOBP_ventricle/GOBP_',prefix,'_filt_LV_new_0.05.csv'))
  
  # Enrichmate 
  gobp_tmp_pos = fgseaResTidy_p[fgseaResTidy_p$NES > 0 & fgseaResTidy_p$padj < 0.05 ,]
  if ((gobp_tmp_pos %>% dim %>% .[1]) > 5 ){
      pathwaytoterm_use = pathwaytoterm[pathwaytoterm$gs_name %in% gobp_tmp_pos$pathway,]
      tb = pathwaytoterm_use; rm(pathwaytoterm_use) 
      colnames(tb)[2] = 'GOID'
      go_id = AnnotationDbi::select(GO.db, keys=tb$GOID, columns = c('GOID', 'TERM', 'ONTOLOGY'), keytype = 'GOID')
      id = go_id$GOID
      test <- enrichmate(id, plot.title = 'LAL enriched',font.size = 14)
      p = test$plot
      saveRDS(test,paste0('./GOBP_ventricle_LR/GOBP_',prefix,'_enrichmate_L_LV_new_0.05.rds'))
  
      png(paste0('./GOBP_ventricle/GOBP_',prefix,'_enrichmate_L_LV_new_0.05.png'), width  = 20*100,height = 20*100)
      print(p)
      dev.off()
  }
  gobp_tmp_neg = fgseaResTidy_p[fgseaResTidy_p$NES < 0 & fgseaResTidy_p$padj < 0.05 ,]
  if ((gobp_tmp_neg %>% dim %>% .[1]) > 5 ){
      pathwaytoterm_use = pathwaytoterm[pathwaytoterm$gs_name %in% gobp_tmp_neg$pathway,]
      tb = pathwaytoterm_use; rm(pathwaytoterm_use) 
      colnames(tb)[2] = 'GOID'
      go_id = AnnotationDbi::select(GO.db, keys=tb$GOID, columns = c('GOID', 'TERM', 'ONTOLOGY'), keytype = 'GOID')
      id = go_id$GOID
      test <- enrichmate(id, plot.title = 'Sham enriched',font.size = 14)
      p = test$plot

      saveRDS(test,paste0('./GOBP_ventricle/GOBP_',prefix,'_enrichmate_S_LV_new_0.05.rds'))
      png(paste0('./GOBP_ventricle/GOBP_',prefix,'_enrichmate_S_LV_new_0.05.png'), width  = 20*100,height = 20*100)
      print(p)
      dev.off()
  }

                                                                                             
}