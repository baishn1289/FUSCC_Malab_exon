
clin.EOCRC_nonhyper <- subsetMaf(maf=rt, tsb=TMB_EOCRC_nonhyper_sample, isTCGA=FALSE)
clin.LOCRC_nonhyper <- subsetMaf(maf=rt, tsb=TMB_LOCRC_nonhyper_sample, isTCGA=FALSE)

Clinical_data_nonhyper <- rt@clinical.data[Tumor_Sample_Barcode %in% c(TMB_EOCRC_nonhyper_sample,
                                                                       TMB_LOCRC_nonhyper_sample)]


save(clin.EOCRC_nonhyper,clin.LOCRC_nonhyper,Clinical_data_nonhyper,
     file = paste0(dir_data,'nonhyper_data.Rdata'))

load(paste0(dir_data,'nonhyper_data.Rdata'))

nonhyper_genes <- uniquegenes(clin.EOCRC_nonhyper,clin.LOCRC_nonhyper,minMut = 5)


logistic_nonhyper <- list()
cl <- makeCluster(12)
registerDoParallel(cl)

logistic_nonhyper <-  foreach(gene = nonhyper_genes, .combine =  c, 
                              .packages = c("data.table", "dplyr","mice","stats")) %dopar% {
                                logstic_mafcompare(m1 = clin.EOCRC_nonhyper, m2 = clin.LOCRC_nonhyper, 
                                                   Clinical.data = Clinical_data_nonhyper, Gene_logistic = gene)
                              }

logistic_nonhyper <- rbindlist(logistic_nonhyper)
stopCluster(cl) 
gc()


result_logistic_nonhyper <- logistic_nonhyper  %>% 
  mutate(mean_freq = (EOCRC_freq+LOCRC_freq)/2) %>% 
  filter(mean_freq > 0.01 ,
         Formula == 'Gene_status ~ Age_class + Panel + Sex + Tumor_site + Sample_type + Histology_type') %>% 
  mutate(adjPval = p.adjust(p = pval, method = "BH"),
         mutation_rate_adjPval =  p.adjust(p = mutation_rate_pval, method = "BH"))


result_logistic_publish_nonhyper <- result_logistic_nonhyper %>% 
  mutate(across(c(OR, `95% upper limit`, `95% lower limit`), ~ round(., digits = 2)),
         across(c(pval, adjPval, mutation_rate_pval,mutation_rate_adjPval), ~ format_Pval(.))
  ) 

nrow(result_logistic_nonhyper[adjPval<0.05])

rownames(result_logistic_nonhyper) <- seq_len(nrow(result_logistic_nonhyper))

save(logistic_nonhyper,
     file=paste0(dir_data,'7_nonhyper_logistic_result.Rdata'))

write.csv(result_logistic_nonhyper,
          file=paste0(dir,'table/7_nonhyper_EOCRC_vs_LOCRC_freq.csv'),
          row.names=FALSE)

write.csv(result_logistic_publish_nonhyper, 
          file=paste0(dir,'table/7_publish_nonhyper_EOCRC_vs_LOCRC_freq.csv'))

filtered_df_nonhyper <- fread(paste0(dir,'table/7_nonhyper_EOCRC_vs_LOCRC_freq.csv')) %>% 
  filter(!is.infinite(OR) & !is.infinite(`95% lower limit`) & !is.infinite(`95% upper limit`)) %>% 
  filter(adjPval <= 0.05) %>%
  mutate( Freq_average = (LOCRC_freq+EOCRC_freq )/2  ) %>%
  arrange(desc(Freq_average)) 

write.csv(filtered_df_nonhyper,
          file = paste0(dir,'table/7.nonhyper_differential_gene.csv'),
          row.names = F)

filtered_df_nonhyper <-filtered_df_nonhyper %>% head(20)


# add new labels
plot <- filtered_df_nonhyper %>% 
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


p_mid <- filtered_df_nonhyper %>% 
  ggplot(aes(y = fct_rev(Hugo_Symbol))) +
  theme_classic() +
  geom_point(aes(x = OR, size = OR), shape = 23, fill = "black", show.legend = F) +
  geom_errorbarh(aes(xmin = `95% lower limit`, xmax = `95% upper limit`), height = 0.3) +
  labs(x = "Odds Ratio") +
  scale_x_continuous(trans = 'log10',breaks = c(0.25,0.5,2.0),
                     labels = scales::number_format(accuracy = 0.05)) +
  coord_cartesian(ylim = c(1, 11),xlim = c(0.25,2)) +
  # scale_x_continuous(trans = 'log10', breaks = c(1.0, 1.5, 2.0),
  #                    labels = scales::number_format(accuracy = 0.01)) +
  # coord_cartesian(ylim = c(1, 3),xlim = c(1,2)) +
  
  geom_vline(xintercept = c(1,mean(filtered_df_nonhyper %>% filter(OR<1) %>% pull(OR))), linetype = c(1,2)) +  
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
  patchwork::area(t = 0, l = 14, b = 30, r = 26))

p_nonhyper =p_left + p_mid + p_right + plot_layout(design = layout)

pdf(file=paste0(dir,'picture/7_nonhyper_forestplot.pdf'), 
    width = 14, height =5)

p_nonhyper
dev.off()
