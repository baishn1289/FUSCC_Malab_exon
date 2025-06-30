
custom_round <- function(x, digits) {
  posneg <- sign(x)   
  z <- abs(x)*10^digits  
  z <- z + 0.5   
  z <- trunc(z)   
  z <- z/10^digits 
  z * posneg  
}


format_frequency <- function(freq) {
  rounded_value <- custom_round(freq * 100, 1)
  paste0(format(rounded_value, nsmall = 1), "%")
  
}


format_Pval <- function(p) {
  ifelse(p < 0.0001,
         "<0.0001",
         case_when(
           p<0.001 ~ sprintf("%.5f", round(p*100000)/100000),
           p<0.01 ~sprintf("%.4f", round(p*10000)/10000),
           p<0.1 ~ sprintf("%.3f", round(p*1000)/1000),
           TRUE ~ sprintf("%.2f", round(p*100)/100)
         )
  )
}


uniquegenes <- function(m1, m2, m1Name = NULL, m2Name = NULL,
                        minMut = 5){
  m1.gs <- getGeneSummary(x = m1)
  m2.gs <- getGeneSummary(x = m2)
  
  
  m1.genes = as.character(m1.gs[AlteredSamples >= minMut, 
                                Hugo_Symbol])
  m2.genes = as.character(m2.gs[AlteredSamples >= minMut, 
                                Hugo_Symbol])
  uniquegenes =  intersect(m1.genes, m2.genes)
  
  return(uniquegenes)
}


mutation_data <- rt@data[, 2:30]
dim(mutation_data)

gene_mutation_counts <- mutation_data %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  reframe(Mutation_Count = n()) %>%
  pivot_wider(names_from = Hugo_Symbol, 
              values_from = Mutation_Count,
              values_fill = list(Mutation_Count = 0)) %>% 
  left_join(rbind(TMB_LOCRC_mutload,TMB_EOCRC_mutload) %>% dplyr::select(Tumor_Sample_Barcode,TMB), by = 'Tumor_Sample_Barcode') %>%
  left_join(rt@clinical.data[,.(Tumor_Sample_Barcode,Age_class)], by = 'Tumor_Sample_Barcode')

# Generate gene lists for each sequencing assay/panel based on the panel-specific gene information 
# provided by AACR and the hg19 reference genome, and load them into the environment
assays_genes <- read_xlsx("panel_genelist_modified.xlsx")%>% 
  na.omit() %>% 
  filter(Panel %in% rt@clinical.data$Panel) %>% 
  group_by(Panel) %>% 
  filter(!duplicated(Hugo_Symbol)) %>% 
  ungroup() %>% 
  as.data.frame()

gene_panel_map <- assays_genes %>% 
  group_by(Hugo_Symbol) %>% 
  reframe(Genes = list(Panel)) %>%
  deframe()  

logstic_mafcompare <- function(m1, m2, m1Name = NULL, m2Name = NULL, Clinical.data,Gene_logistic) {
  
  Alter_m1_Gene <- as.data.table(m1@data)[Hugo_Symbol == Gene_logistic]
  Alter_m1_Gene <- unique(Alter_m1_Gene, by = "Tumor_Sample_Barcode")
  
  Alter_m2_Gene <- as.data.table(m2@data)[Hugo_Symbol == Gene_logistic]
  Alter_m2_Gene <- unique(Alter_m2_Gene, by = "Tumor_Sample_Barcode")
  
  Gene_panels <- Clinical.data %>% 
    filter(Panel %in% gene_panel_map[[Gene_logistic]]) %>% 
    pull(Panel)
  
  mut_samples_m1 <- length(unique(Alter_m1_Gene$Tumor_Sample_Barcode))
  panels_samples_m1 <- length(unique(m1@clinical.data[m1@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  m1_freq <- mut_samples_m1 / panels_samples_m1
  
  mut_samples_m2 <- length(unique(Alter_m2_Gene$Tumor_Sample_Barcode))
  panels_samples_m2 <- length(unique(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  m2_freq <- mut_samples_m2 / panels_samples_m2
  
  if ( m1_freq < 0.01 | m2_freq<0.01 ) {
    print(paste(Gene_logistic,
                "Low mutation frequencies, returning from function.",
                sep = '____'))
    return()
  }
  
  imp_list_Gene_clinical_info <-lapply(imp_list, function(df) {
    dt <- as.data.table(df)
    dt <- dt[Tumor_Sample_Barcode %in% Clinical.data$Tumor_Sample_Barcode]
    dt[, Gene_status := as.integer(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode, Alter_m2_Gene$Tumor_Sample_Barcode))]
    dt[, Gene_status := factor(Gene_status, levels = c(0, 1))]
    dt <- dt[Panel %in% Gene_panels]
    return(dt)
  })
  
  valid_vars <- "Age_class"
  covariates <- c( "Panel", "Sex", "Tumor_site", "Sample_type", "Histology_type")
  for (cov in covariates) {
    cov_exists <- all(sapply(imp_list_Gene_clinical_info, function(df) cov %in% colnames(df)))
    if (!cov_exists) next
    sufficient_levels <- all(sapply(imp_list_Gene_clinical_info, function(df) {
      length(unique(df[[cov]])) >= 2
    }))
    if (sufficient_levels) {
      valid_vars <- c(valid_vars, cov)
    }
  }
  
  if (length(valid_vars) == 0) {
    formula_base <- "Gene_status ~ 1" 
  } else {
    formula_base <- paste("Gene_status ~", paste(valid_vars, collapse = " + "))
  }
  
  logistic_models <- lapply(imp_list_Gene_clinical_info, function(df) {
    tryCatch({
      glm(as.formula(formula_base), family = binomial(link = "logit"), data = df)
    }, error = function(e) {
      return(NULL)
    })
  })
  
  logistic_models <- Filter(Negate(is.null), logistic_models)
  
  if (length(logistic_models) == 0) {
    return()
  }
  
  
  summary_pooled <- summary(pool(logistic_models))
  
  locrc_eocrc_results <- summary_pooled %>% 
    filter(term == "Age_classEOCRC") %>% 
    mutate(conf.low = estimate - 1.959964 * std.error,
           conf.high = estimate + 1.959964 * std.error)%>% 
    mutate(OR = exp(estimate),
           `95% lower limit` = exp(conf.low),
           `95% upper limit` = exp(conf.high))
  
  each_gene_mutation_counts = gene_mutation_counts %>% 
    dplyr::select(all_of(Gene_logistic), Tumor_Sample_Barcode,Age_class, TMB) %>%
    filter(Tumor_Sample_Barcode %in% imp_list_Gene_clinical_info[[1]]$Tumor_Sample_Barcode)  %>% 
    mutate(mutation_rate = .data[[Gene_logistic]] / TMB)
  
  each_gene_mutation_EO = each_gene_mutation_counts %>% 
    filter(Age_class == 'EOCRC') %>% pull(mutation_rate)
  each_gene_mutation_LO = each_gene_mutation_counts%>% 
    filter(Age_class == 'LOCRC')%>% pull(mutation_rate)
  
  EOCRC_rate = mean(each_gene_mutation_EO)
  LOCRC_rate = mean(each_gene_mutation_LO)
  
  
  mutation_rate_pvalue <- wilcox.test(each_gene_mutation_EO, each_gene_mutation_LO)$p.value
  
  df_Gene_test <- data.table::data.table(
    Hugo_Symbol = Gene_logistic,
    EOCRC_total_case= panels_samples_m1,
    LOCRC_total_case= panels_samples_m2,
    EOCRC = mut_samples_m1, 
    LOCRC = mut_samples_m2, 
    EOCRC_freq = m1_freq,
    LOCRC_freq = m2_freq,
    pval = locrc_eocrc_results$p.value, 
    OR = locrc_eocrc_results$OR,
    `95% upper limit` = locrc_eocrc_results$`95% upper limit`,
    `95% lower limit` = locrc_eocrc_results$`95% lower limit`,
    EOCRC_mutation_rate = EOCRC_rate,
    LOCRC_mutation_rate = LOCRC_rate,
    mutation_rate_pval = mutation_rate_pvalue,
    Formula = formula_base
  )
  return(list(df_Gene_test))
  
}


overall_genes <- uniquegenes(clin.EOCRC,clin.LOCRC,minMut = 5)

logistic_overall <- list()
cl <- makeCluster(12)
registerDoParallel(cl)

logistic_overall <-  foreach(gene = overall_genes, .combine =  c, 
                             .packages = c("data.table", "dplyr","mice","stats")) %dopar% {
                               logstic_mafcompare(m1 = clin.EOCRC, m2 = clin.LOCRC, 
                                                  Clinical.data = rt@clinical.data, Gene_logistic = gene)
                             }

logistic_overall <- rbindlist(logistic_overall)
stopCluster(cl)
registerDoSEQ()
gc()

result_logistic <- logistic_overall  %>% 
  mutate(mean_freq = (EOCRC_freq+LOCRC_freq)/2) %>% 
  filter(mean_freq > 0.01,
         Formula == 'Gene_status ~ Age_class + Panel + Sex + Tumor_site + Sample_type + Histology_type' ) %>% 
  mutate(adjPval = p.adjust(p = pval, method = "BH"),
         mutation_rate_adjPval =  p.adjust(p = mutation_rate_pval, method = "BH"))

result_logistic_publish <- result_logistic %>% 
  mutate(across(c(OR, `95% upper limit`, `95% lower limit`), ~ round(., digits = 2)),
         across(c(pval, adjPval, mutation_rate_pval,mutation_rate_adjPval), ~ format_Pval(.))
  ) 

nrow(result_logistic[adjPval<0.05])

rownames(result_logistic) <- seq_len(nrow(result_logistic))

save(logistic_overall,
     file=paste0(dir_data,'3.0_overall_logistic_result.Rdata'))

write.table(result_logistic, file=paste0(dir,'table/3.0_overall_EOCRC_vs_LOCRC_freq.csv'), 
            quote=FALSE, row.names=FALSE, sep="\t")
write.csv(result_logistic_publish, 
          file=paste0(dir,'table/3.0_publish_overall_EOCRC_vs_LOCRC_freq.csv'))


#overall_forest
{
  filtered_df_overall <- fread(paste0(dir,'table/3.0_overall_EOCRC_vs_LOCRC_freq.csv')) %>% 
    filter(!is.infinite(OR) & !is.infinite(`95% lower limit`) & !is.infinite(`95% upper limit`)) %>% 
    filter(adjPval <= 0.05) %>%
    mutate( Freq_average = (LOCRC_freq+EOCRC_freq )/2  ) %>%
    arrange(desc(Freq_average)) 
  
  write.csv(filtered_df_overall,
            file = paste0(dir,'table/3.0_overall_differential_gene.csv'),
            row.names = F)
  
  filtered_df_overall <-filtered_df_overall %>% head(20)
  
  # add new labels
  plot <- filtered_df_overall %>% 
    mutate(EOCRC_freq = format_frequency(EOCRC_freq),
           LOCRC_freq = format_frequency(LOCRC_freq),
           sample_size = paste0(LOCRC_freq, " vs. ", EOCRC_freq),
           OR_CI = paste0(format(round(OR, 2), nsmall = 2), " (", format(round(`95% lower limit`, 2), nsmall = 2), " to ", format(round(`95% upper limit`, 2), nsmall = 2), ")"),
           adjPval = format_Pval(adjPval),
           EOCRC_mutation_rate = format_frequency(EOCRC_mutation_rate),
           LOCRC_mutation_rate = format_frequency(LOCRC_mutation_rate),
           mutation_rate = paste0(LOCRC_mutation_rate, " vs. ", EOCRC_mutation_rate),
           mutation_rate_adjPval = format_Pval(mutation_rate_adjPval)
    ) %>% 
    mutate(across(everything(), as.character)) %>% 
    bind_rows(data.frame(Hugo_Symbol = "Gene", 
                         sample_size = "Frequency\n(LOCRC vs. EOCRC)", 
                         adjPval = "adjPval", 
                         OR_CI = "OR (95% CI)",
                         mutation_rate = 'mutation_rate\n(LOCRC vs. EOCRC)',
                         mutation_rate_adjPval = 'mutation_rate_adjPval')) %>% 
    mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
  
  # middle part
  p_mid <- filtered_df_overall %>% 
    ggplot(aes(y = fct_rev(Hugo_Symbol))) +
    theme_classic() +
    geom_point(aes(x = OR, size = OR), shape = 23, fill = "black", show.legend = F) +
    geom_errorbarh(aes(xmin = `95% lower limit`, xmax = `95% upper limit`), height = 0.3) +
    labs(x = "Odds Ratio") +
    scale_x_continuous( trans = 'log10',breaks = c(0.25,0.5,1,2,4),
                        labels = scales::number_format(accuracy = 0.05)) +
    coord_cartesian(ylim = c(1, 21), xlim = c(0.25,4)) +
    geom_vline(xintercept = c(1#,
                              # mean(filtered_df_overall %>% filter(OR<1) %>% pull(OR)),
                              # mean(filtered_df_overall %>% filter(OR>1) %>% pull(OR))
    ), linetype = 1) +  
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(color = "black"))
  
  p_left <- plot %>%
    ggplot(aes(y = Hugo_Symbol)) +
    geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
    geom_text(aes(x = 1, label = sample_size), hjust = 0.5, size = 4,vjust=0.5,
              fontface = ifelse(plot$sample_size == "Frequency\n(LOCRC vs. EOCRC)", "bold", "plain")) +
    theme_void() + coord_cartesian(xlim = c(0, 1.5))
  
  # right_part
  p_right <- plot %>% ggplot() +
    geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
              fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
    geom_text(aes(x = 1.2, y = Hugo_Symbol, label = adjPval), hjust = 0,
              fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
    geom_text(aes(x = 2.3, y = Hugo_Symbol, label = mutation_rate), hjust = 0.5,
              fontface = ifelse(plot$mutation_rate == 'mutation_rate\n(LOCRC vs. EOCRC)', "bold", "plain"), size = 4) +
    geom_text(aes(x = 3, y = Hugo_Symbol, label = mutation_rate_adjPval), hjust = 0,
              fontface = ifelse(plot$mutation_rate_adjPval == "mutation_rate_adjPval", "bold", "plain"), size = 4) +
    theme_void() + coord_cartesian(xlim = c(0, 4))
  
  
  layout <- c(
    patchwork::area(t = 0, l = 0, b = 30, r = 6), 
    patchwork::area(t = 0, l = 7, b = 30, r = 13), 
    patchwork::area(t = 0, l = 14.5, b = 30, r = 26.5))
  
}


pdf(paste0(dir,'picture/3_overall_forestPlot.pdf'), 
    width = 14, height = 8)
print(p_left + p_mid + p_right + plot_layout(design = layout))

dev.off()
