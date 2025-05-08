
Gene_list_country_LO <- merge(
  setDT(clin.LOCRC@clinical.data[, .(Country, Panel)]), 
  setDT(assays_genes),   by = "Panel",   allow.cartesian = TRUE) %>% 
  distinct(Country, Hugo_Symbol) %>%  
  group_by(Country) %>% 
  summarize(Gene_list = list(unique(Hugo_Symbol)), .groups = "drop")

gene_lists <- lapply(Gene_list_country_LO$Gene_list, unlist)
common_genes <- Reduce(intersect, gene_lists)

mutation_table <- data.table(
  gene = common_genes,
  Canada = NA_real_,
  China = NA_real_,
  US = NA_real_,
  Netherlands = NA_real_,
  Nigeria = NA_real_,
  France = NA_real_,
  Spain = NA_real_
)

countries_LO <- unique(rt@clinical.data$Country)[1:7]

for (i in common_genes) {
  generate_country_frequency(countries = countries_LO, 
                             maf = clin.LOCRC, Gene_common = i)
}

mutation_table_LO <- mutation_table %>% 
  filter(rowSums(across(c("Canada", "China", "US", 
                          "Netherlands", "Nigeria", 
                          "France", "Spain"), ~ .x > 0)) >= 7) %>% 
  mutate(Korea = 0) %>% 
  mutate(Age_class = 'LOCRC')

write.csv(mutation_table_LO,
          file = paste0(dir,'table/5.0_LOCRC_Country_commongenes_mutation_table_LO.csv'))


common_EO_LO_gene <- intersect(mutation_table_EO$gene,
                               mutation_table_LO$gene)

common_EO_LO_df <- rbind(mutation_table_EO[mutation_table_EO$gene %in%common_EO_LO_gene,],
                         mutation_table_LO[mutation_table_LO$gene %in%common_EO_LO_gene,]) %>% 
  pivot_longer(cols = Canada:Korea, 
               names_to = "Country", 
               values_to = "Frequency")


palette <- c('#1f78b4','#fc8d62')

plot_common_EO_LO <- common_EO_LO_df %>% 
  ggplot(aes(Frequency, Country, fill = Age_class)) +
  geom_col(position = position_dodge(1), width = 0.6) + 
  labs(y = NULL) +
  facet_wrap('gene', scales = "free", nrow  = 3) + 
  theme_few() +
  scale_fill_manual(values = palette) +
  theme(
    axis.text.y = element_text(color = "black", size = 8), # y轴文本显示
    # axis.text.x = element_text(color = "black", size = 8),
    axis.text.x   = element_text(color = "black", size = 8, hjust=1,angle = 45),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    axis.ticks.y = element_blank(), 
    legend.position = "right", 
    plot.background = element_blank(), 
    panel.background = element_blank(), 
    strip.background = element_blank(), 
    strip.placement = "outside", 
    strip.text.y.left = element_text(angle = 0) 
  )+
  coord_flip()

pdf(paste0(dir,"picture/4_EOCRC_LOCRC_common_genes_mainplot.pdf"), 
    width = 10, height = 5.5)
plot_common_EO_LO
dev.off()


palette <- c('Canada' = '#b2182b',
             'China' = '#d6604d',
             'France' = '#f4a582',
             'Korea' = '#fddbc7',
             'Netherlands' = '#d1e5f0',
             'Nigeria' = '#92c5de',
             'Spain' = '#4393c3',
             'US'='#2166ac') 

mutation_table_LO_sup <- mutation_table_LO %>% 
  filter(!(gene %in% common_EO_LO_df$gene)) %>% 
  dplyr::select(-Age_class)


plotdata_LO <- data.frame(t(mutation_table_LO_sup))
colnames(plotdata_LO) <- plotdata_LO[1,]
plotdata_LO = plotdata_LO[-1,]
plotdata_LO$Country = rownames(plotdata_LO)

mergedata <- reshape2::melt(plotdata_LO[,], id = 'Country') 
mergedata$value = round(as.numeric(mergedata$value), 3)

Common_country_LO = ggplot(mergedata, aes(x = Country, y = value, fill = Country)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ variable, scales = "free", nrow =3) +
  theme_few() +
  scale_fill_manual(values = palette) +
  labs(subtitle = "LOCRC")+
  theme(
    plot.title    = element_text(color = "black", size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = "black", size   = 16,hjust = 0.5),
    plot.caption  = element_text(color = "black", size   = 16,face = "italic", hjust = 1),
    axis.text.x   = element_text(color = "black", size = 16, hjust=1,angle = 45),
    axis.text.y   = element_text(color = "black", size = 16, angle = 0),
    axis.title.x  = element_text(color = "black", size = 20, angle = 0),
    axis.title.y  = element_text(color = "black", size = 20, angle = 90),
    legend.title  = element_text(color = "black", size  = 16),
    legend.text   = element_text(color = "black", size   = 16),
    strip.text = element_text(size = 14),
    axis.line.y = element_line(color = "black", linetype = "solid"),
    axis.line.x = element_line (color = "black",linetype = "solid"),
    panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA) )

Common_country_LO
ggsave(paste0(dir,"picture/4_LOCRC_common_genes_supplementary.pdf"), 
       width = 24, height = 10, dpi = 300)


mutation_table_EO_sup <- mutation_table_EO %>% 
  filter(!(gene %in% common_EO_LO_df$gene)) %>% 
  dplyr::select(-Age_class)


plotdata_EO <- data.frame(t(mutation_table_EO_sup))

colnames(plotdata_EO) <- plotdata_EO[1,]
plotdata_EO = plotdata_EO[-1,]
plotdata_EO$Country = rownames(plotdata_EO)

mergedata <- reshape2::melt(plotdata_EO[,], id = 'Country') 
mergedata$value = round(as.numeric(mergedata$value), 3)


Common_country_EO = ggplot(mergedata, aes(x = Country, y = value, fill = Country)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ variable, scales = "free", nrow =1) +
  theme_few() +
  scale_fill_manual(values = palette) +
  labs(subtitle = "EOCRC")+
  theme(
    plot.title    = element_text(color = "black", size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = "black", size   = 16,hjust = 0.5),
    plot.caption  = element_text(color = "black", size   = 16,face = "italic", hjust = 1),
    axis.text.x   = element_text(color = "black", size = 16, hjust=1,angle = 45),
    axis.text.y   = element_text(color = "black", size = 16, angle = 0),
    axis.title.x  = element_text(color = "black", size = 20, angle = 0),
    axis.title.y  = element_text(color = "black", size = 20, angle = 90),
    legend.title  = element_text(color = "black", size  = 16),
    legend.text   = element_text(color = "black", size   = 16),
    strip.text = element_text(size = 14),
    axis.line.y = element_line(color = "black", linetype = "solid"),
    axis.line.x = element_line (color = "black",linetype = "solid"),
    panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA) )

Common_country_EO
ggsave(paste0(dir,"picture/4_EOCRC_common_genes_supplementary.pdf"), 
       width = 11, height = 4, dpi = 300)


  