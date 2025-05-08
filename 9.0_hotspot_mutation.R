
logstic_mafcompare_hotspot <- function(m1, m2, Clinical.data, Gene_logistic, aachange, TMB_Status) {
  
  m1_AA_changed_samples <- as.data.table(m1@data)[Hugo_Symbol == Gene_logistic & HGVSp_Short == aachange,]
  m1_AA_changed_samples <- unique(m1_AA_changed_samples, by = "Tumor_Sample_Barcode")
  
  m2_AA_changed_samples <- as.data.table(m2@data)[Hugo_Symbol == Gene_logistic & HGVSp_Short == aachange,]
  m2_AA_changed_samples <- unique(m2_AA_changed_samples, by = "Tumor_Sample_Barcode")
  
  Gene_panels <- Clinical.data  %>% 
    filter(Panel %in% gene_panel_map[[Gene_logistic]]) %>% 
    pull(Panel)
  
  panels_samples_m1 <- length(unique(m1@clinical.data[m1@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode)) 
  m1_freq <- nrow(m1_AA_changed_samples) / panels_samples_m1
  
  
  panels_samples_m2 <- length(unique(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
  m2_freq <- nrow(m2_AA_changed_samples) / panels_samples_m2
  
  
  imp_list_Gene_clinical_info <-lapply(imp_list, function(df) {
    dt <- as.data.table(df)
    dt <- dt[Tumor_Sample_Barcode %in% Clinical.data$Tumor_Sample_Barcode]
    dt[, Gene_status := as.integer(Tumor_Sample_Barcode %in% c(m1_AA_changed_samples$Tumor_Sample_Barcode, 
                                                               m2_AA_changed_samples$Tumor_Sample_Barcode))]
    dt[, Gene_status := factor(Gene_status, levels = c(0, 1))]
    dt <- dt[Panel %in% Gene_panels]
    return(dt)
  })
  
  
  each_gene_mutation_counts = gene_mutation_counts %>% 
    dplyr::select(Tumor_Sample_Barcode,Age_class, TMB) %>% 
    filter(Tumor_Sample_Barcode %in% imp_list_Gene_clinical_info[[1]]$Tumor_Sample_Barcode) %>%
    mutate(hotspot = ifelse(Tumor_Sample_Barcode %in% c(m1_AA_changed_samples$Tumor_Sample_Barcode,
                                                        m2_AA_changed_samples$Tumor_Sample_Barcode), 1, 0)) %>% 
    mutate(mutation_rate = hotspot / TMB)
  
  each_gene_mutation_EO = each_gene_mutation_counts %>% 
    filter(Age_class == 'EOCRC') %>% pull(mutation_rate)
  each_gene_mutation_LO = each_gene_mutation_counts%>% 
    filter(Age_class == 'LOCRC')%>% pull(mutation_rate)
  
  EOCRC_rate = mean(each_gene_mutation_EO)
  LOCRC_rate = mean(each_gene_mutation_LO)
  
  mutation_rate_pvalue <- wilcox.test(each_gene_mutation_EO, each_gene_mutation_LO)$p.value
  
  logistic_models <-  lapply(imp_list_Gene_clinical_info, function(df) {
    glm(as.formula("Gene_status ~ Age_class + Sex + Tumor_site + Panel+ Histology_type+Sample_type"), family = binomial(link = "logit"),
        data = df)
    
  })
  
  summary_pooled <- summary(pool(logistic_models))
  
  
  locrc_eocrc_results <- summary_pooled %>% 
    filter(term == "Age_classEOCRC") %>% 
    mutate(conf.low = estimate - 1.959964 * std.error,
           conf.high = estimate + 1.959964 * std.error)%>% 
    mutate(OR = exp(estimate),
           `95% lower limit` = exp(conf.low),
           `95% upper limit` = exp(conf.high))
  
  
  df_Gene_test <- data.table::data.table(
    Hugo_Symbol = Gene_logistic,
    Hotspot = aachange,
    EOCRC_total_case= panels_samples_m1,
    LOCRC_total_case= panels_samples_m2,
    EOCRC = nrow(m1_AA_changed_samples), 
    LOCRC = nrow(m2_AA_changed_samples), 
    EOCRC_freq = m1_freq,
    LOCRC_freq = m2_freq,
    pval = locrc_eocrc_results$p.value, 
    OR = locrc_eocrc_results$OR,
    `95% upper limit` = locrc_eocrc_results$`95% upper limit`,
    `95% lower limit` = locrc_eocrc_results$`95% lower limit`,
    EOCRC_mutation_rate = EOCRC_rate,
    LOCRC_mutation_rate = LOCRC_rate,
    mutation_rate_pval = mutation_rate_pvalue,
    TMB_status = TMB_Status
  )
  
  return(df_Gene_test)
  
}


mut_freq <- table(rt@data[Hugo_Symbol == 'TP53']$HGVSp_Short, useNA = 'always')
sorted_mut_freq <- sort(mut_freq, decreasing = TRUE)
head(sorted_mut_freq,10)

TMB_statuses <- c("_overall", "_hyper", "_nonhyper")

clin.EOCRC_overall <- clin.EOCRC
clin.LOCRC_overall <- clin.LOCRC

BRAF_lists <- lapply(TMB_statuses, function(status) {
  clin_data <- copy(switch(status,
                           "_overall" = rt@clinical.data,
                           "_hyper" = rbind(clin.EOCRC_hyper@clinical.data, clin.LOCRC_hyper@clinical.data),
                           "_nonhyper" = rbind(clin.EOCRC_nonhyper@clinical.data, clin.LOCRC_nonhyper@clinical.data)))
  
  
  logstic_mafcompare_hotspot(
    m1 = get(paste0("clin.EOCRC", tolower(status))), 
    m2 = get(paste0("clin.LOCRC", tolower(status))),
    Clinical.data = clin_data,
    Gene_logistic = 'BRAF',
    aachange = 'p.V600E',
    TMB_Status = status
  )
})


KRAS_12D_lists <- lapply(TMB_statuses, function(status) {
  clin_data <- copy(switch(status,
                           "_overall" = rt@clinical.data,
                           "_hyper" = rbind(clin.EOCRC_hyper@clinical.data, clin.LOCRC_hyper@clinical.data),
                           "_nonhyper" = rbind(clin.EOCRC_nonhyper@clinical.data, clin.LOCRC_nonhyper@clinical.data)))
  
  
  logstic_mafcompare_hotspot(
    m1 = get(paste0("clin.EOCRC", tolower(status))), 
    m2 = get(paste0("clin.LOCRC", tolower(status))),
    Clinical.data = clin_data,
    Gene_logistic = 'KRAS',
    aachange = 'p.G12D',
    TMB_Status = status
  )
})


TP53_R175H_lists <- lapply(TMB_statuses, function(status) {
  clin_data <- copy(switch(status,
                           "_overall" = rt@clinical.data,
                           "_hyper" = rbind(clin.EOCRC_hyper@clinical.data, clin.LOCRC_hyper@clinical.data),
                           "_nonhyper" = rbind(clin.EOCRC_nonhyper@clinical.data, clin.LOCRC_nonhyper@clinical.data)))
  
  
  logstic_mafcompare_hotspot(
    m1 = get(paste0("clin.EOCRC", tolower(status))), 
    m2 = get(paste0("clin.LOCRC", tolower(status))),
    Clinical.data = clin_data,
    Gene_logistic = 'TP53',
    aachange = 'p.R175H',
    TMB_Status = status
  )
})


TP53_R248Q_lists <- lapply(TMB_statuses, function(status) {
  clin_data <- copy(switch(status,
                           "_overall" = rt@clinical.data,
                           "_hyper" = rbind(clin.EOCRC_hyper@clinical.data, clin.LOCRC_hyper@clinical.data),
                           "_nonhyper" = rbind(clin.EOCRC_nonhyper@clinical.data, clin.LOCRC_nonhyper@clinical.data)))
  
  
  logstic_mafcompare_hotspot(
    m1 = get(paste0("clin.EOCRC", tolower(status))), 
    m2 = get(paste0("clin.LOCRC", tolower(status))),
    Clinical.data = clin_data,
    Gene_logistic = 'TP53',
    aachange = 'p.R248Q',
    TMB_Status = status
  )
})

KRAS_G12C_lists <- lapply(TMB_statuses, function(status) {
  clin_data <- copy(switch(status,
                           "_overall" = rt@clinical.data,
                           "_hyper" = rbind(clin.EOCRC_hyper@clinical.data, clin.LOCRC_hyper@clinical.data),
                           "_nonhyper" = rbind(clin.EOCRC_nonhyper@clinical.data, clin.LOCRC_nonhyper@clinical.data)))
  
  
  logstic_mafcompare_hotspot(
    m1 = get(paste0("clin.EOCRC", tolower(status))), 
    m2 = get(paste0("clin.LOCRC", tolower(status))),
    Clinical.data = clin_data,
    Gene_logistic = 'KRAS',
    aachange = 'p.G12C',
    TMB_Status = status
  )
})


BRAF_lists <- rbindlist(BRAF_lists)
KRAS_12D_lists <- rbindlist(KRAS_12D_lists)
TP53_R175H_lists <- rbindlist(TP53_R175H_lists)
TP53_R248Q_lists <- rbindlist(TP53_R248Q_lists)
KRAS_G12C_lists <- rbindlist(KRAS_G12C_lists)

result_df <- rbind(BRAF_lists,KRAS_12D_lists,KRAS_G12C_lists,TP53_R175H_lists)
result_hotspot =  result_df %>% 
  mutate(TMB_status = case_when(
    TMB_status =="_overall"~ "Overall",
    TMB_status =="_hyper" ~ 'Hypermutated',
    TMB_status == "_nonhyper" ~ 'Non-hypermutated'
  )) %>% 
  group_by(TMB_status) %>% 
  filter(EOCRC >0 & LOCRC >0) %>% 
  mutate(adjPval = p.adjust(p = pval, method = "BH"),
         mutation_rate_adjPval =  p.adjust(p = mutation_rate_pval, method = "BH")) %>% 
  ungroup()  %>% 
  mutate(Hugo_Symbol = case_when(
    TMB_status == 'Non-hypermutated' ~ paste(Hugo_Symbol, ' '),  
    TMB_status == 'Hypermutated' ~ paste(Hugo_Symbol, '  '),
    TRUE ~ Hugo_Symbol 
  )) %>% 
  mutate( Hugo_Symbol = case_when(
    Hotspot == 'p.G12C' ~ paste(Hugo_Symbol, '    '),  
    # Hotspot == 'p.G12C' ~ paste(Hugo_Symbol, '  '),
    TRUE ~ Hugo_Symbol  
  )
  )

result_hotspot_publish <- result_hotspot %>% 
  mutate(across(c(OR, `95% upper limit`, `95% lower limit`), ~ round(., digits = 2)),
         across(c(pval, adjPval, mutation_rate_pval,mutation_rate_adjPval), ~ format_Pval(.))
  ) 

write.csv(result_df,paste0(dir,'table/10.hotspot_original_data.csv'))
write.csv(result_hotspot,paste0(dir,'table/10.hotspot.csv'))
write.csv(result_hotspot_publish,paste0(dir,'table/10_publish_hotspot.csv'))

plot <- result_hotspot %>% 
  mutate(EOCRC_freq = format_frequency(EOCRC_freq),
         LOCRC_freq = format_frequency(LOCRC_freq),
         sample_size = paste0(LOCRC_freq, " vs. ", EOCRC_freq),
         OR_CI = paste0(format(round(OR, 2), nsmall = 2), " (", 
                        format(round(`95% lower limit`, 2), nsmall = 2), " to ", 
                        format(round(`95% upper limit`, 2), nsmall = 2), ")"),
         adjPval = format_Pval(adjPval),
         EOCRC_mutation_rate = format_frequency(EOCRC_mutation_rate),
         LOCRC_mutation_rate = format_frequency(LOCRC_mutation_rate),
         mutation_rate = paste0(LOCRC_mutation_rate, " vs. ", EOCRC_mutation_rate),
         mutation_rate_adjPval = format_Pval(mutation_rate_adjPval)
  ) %>% 
  mutate(across(everything(), as.character)) %>% 
  bind_rows(data.frame(Hugo_Symbol = "Gene", 
                       Hotspot = "Hotspot",
                       TMB_status = "TMB_status",
                       sample_size = "Frequency\n(LOCRC vs. EOCRC)", 
                       adjPval = "adjPval", 
                       OR_CI = "OR (95% CI)",
                       mutation_rate = 'mutation_rate\n(LOCRC vs. EOCRC)',
                       mutation_rate_adjPval = 'mutation_rate_adjPval')) %>% 
  mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))

# middle part
p_mid <- result_hotspot %>% 
  ggplot(aes(y = fct_rev(Hugo_Symbol))) +
  
  theme_classic() +
  geom_point(aes(x = OR, size = OR), shape = 23, fill = "black", show.legend = F) +
  geom_errorbarh(aes(xmin = `95% lower limit`, xmax = `95% upper limit`), height = 0.3) +
  labs(x = "Odds Ratio") +
  scale_x_continuous(trans = 'log10',breaks = c(0.01,0.1,1,2,4),
                     labels = scales::number_format(accuracy = 0.01)) +
  coord_cartesian(ylim = c(1,13), xlim = c(0.01,4)) +
  geom_vline(xintercept = 1, linetype = 2) +  
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black"))

# left part
p_left <- plot %>%
  ggplot(aes(y = Hugo_Symbol)) +
  
  geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
  geom_text(aes(x = 0.8, label = Hotspot), hjust = 0, fontface = "bold", size = 4) +
  geom_text(aes(x = 1.6, label = TMB_status), hjust = 0, fontface = "bold", size = 4) +
  geom_text(aes(x = 3, label = sample_size), hjust = 0.5, size = 4,vjust=0.5,
            fontface = ifelse(plot$sample_size == "Frequency\n(LOCRC vs. EOCRC)", "bold", "plain")) +
  theme_void() + coord_cartesian(xlim = c(0, 3.5))

# right_part
p_right <- plot %>% ggplot() +
  
  geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
            fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
  geom_text(aes(x = 1.1, y = Hugo_Symbol, label = adjPval), hjust = 0,
            fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
  geom_text(aes(x = 2.3, y = Hugo_Symbol, label = mutation_rate), hjust = 0.5,
            fontface = ifelse(plot$mutation_rate == 'mutation_rate\n(LOCRC vs. EOCRC)', "bold", "plain"), size = 4) +
  geom_text(aes(x = 3, y = Hugo_Symbol, label = mutation_rate_adjPval), hjust = 0,
            fontface = ifelse(plot$mutation_rate_adjPval == "mutation_rate_adjPval", "bold", "plain"), size = 4) +
  theme_void() + coord_cartesian(xlim = c(0, 4))

layout <- c(
  patchwork::area(t = 0, l = 0, b = 30, r = 10), 
  patchwork::area(t = 0, l = 11, b = 30, r = 16), 
  patchwork::area(t = 0, l = 17, b = 30, r = 28))

p_hotspot =p_left + p_mid + p_right + plot_layout(design = layout)

pdf(paste(dir,'picture/10.hotspot_forestplot.pdf',sep = ""), 
    width = 16, height = 8)

p_hotspot
dev.off()

######BRAF V600E proportion calculation
proportion_calculation_df <- data.frame(
  TMB_Status = character(),
  BRAF_EO = numeric(),
  BRAF_LO =  numeric(),
  BRAF_V600E_EO =  numeric(),
  BRAF_V600E_LO =  numeric(),
  BRAF_V600E_EO_rate =  numeric(),
  BRAF_V600E_LO_rate =  numeric(),
  Pval =  numeric()
)


hotspot_proportion_calculation <- function(clin_data_1, clin_data_2, 
                                           gene = 'BRAF', mutation = 'p.V600E',
                                           status) {
  
  count_1_overall <- nrow(clin_data_1@data[Hugo_Symbol == gene]%>% 
                            filter(!duplicated(Tumor_Sample_Barcode)))
  count_1_mutation <- nrow(clin_data_1@data[Hugo_Symbol == gene & HGVSp_Short == mutation]%>% 
                             filter(!duplicated(Tumor_Sample_Barcode)))
  
  
  count_2_overall <- nrow(clin_data_2@data[Hugo_Symbol == gene]%>% 
                            filter(!duplicated(Tumor_Sample_Barcode)))
  count_2_mutation <- nrow(clin_data_2@data[Hugo_Symbol == gene & HGVSp_Short == mutation]%>% 
                             filter(!duplicated(Tumor_Sample_Barcode)))
  
  
  rate_1 <- count_1_mutation / count_1_overall
  rate_2 <- count_2_mutation / count_2_overall
  cat("Mutation rate for group 1 (", mutation, "): ", rate_1, "\n")
  cat("Mutation rate for group 2 (", mutation, "): ", rate_2, "\n")
  
  
  table <- matrix(c(
    count_1_mutation, count_1_overall - count_1_mutation,
    count_2_mutation, count_2_overall - count_2_mutation
  ), nrow = 2, byrow = TRUE)
  
  if(count_1_mutation <=5  | count_2_mutation <=5){
    chisq_result <- fisher.test(table)
    print('fisher')
  }else{
    chisq_result <- chisq.test(table)
    print('chisq')
  }
  
  
  print(chisq_result)
  
  df <- data.frame(
    TMB_Status = status,
    BRAF_EO = count_1_overall,
    BRAF_LO = count_2_overall,
    BRAF_V600E_EO = count_1_mutation,
    BRAF_V600E_LO = count_2_mutation,
    BRAF_V600E_EO_rate = rate_1,
    BRAF_V600E_LO_rate = rate_2,
    
    pval = chisq_result$p.value
  )
  
  proportion_calculation_df  <<- rbind(proportion_calculation_df,df)
}


hotspot_proportion_calculation(clin_data_1 = clin.EOCRC, 
                               clin_data_2 = clin.LOCRC,
                               status = 'Overall')
hotspot_proportion_calculation(clin_data_1 = clin.EOCRC_hyper, 
                               clin_data_2 = clin.LOCRC_hyper,
                               status = 'Hyper')
hotspot_proportion_calculation(clin_data_1 = clin.EOCRC_nonhyper, 
                               clin_data_2 = clin.LOCRC_nonhyper,
                               status = 'Non-hyper')

proportion_calculation_df <- proportion_calculation_df %>% 
  mutate(adjPval = p.adjust(p = pval, method = "BH")) %>% 
  mutate(across(c(pval, adjPval), ~ format_Pval(.))
  ) 

write.csv(proportion_calculation_df,
          file= paste0(dir,'table/10.V600E_proportion_calculation_df.csv'))
