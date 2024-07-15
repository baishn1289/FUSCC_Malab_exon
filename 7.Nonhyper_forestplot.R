
{
  # TMB nonhypermutated
  clin.EOCRC_nonhyper <- subsetMaf(maf=rt, tsb=TMB_EOCRC_nonhyper_sample, isTCGA=FALSE)
  clin.LOCRC_nonhyper <- subsetMaf(maf=rt, tsb=TMB_LOCRC_nonhyper_sample, isTCGA=FALSE)
  
  fvsm <- mafCompare(m1=clin.EOCRC_nonhyper, m2=clin.LOCRC_nonhyper, m1Name="EOCRC", m2Name="LOCRC", minMut=5)
  write.table(fvsm$results, file=paste(dir,'/nonhyperEOCRC_vs_LOCRC_counts.tsv',sep = ""), quote=FALSE, row.names=FALSE, sep="\t")
  
  fvsm$SampleSummary$SampleSize[1]
  op = fvsm$results
  op$EOCRC_freq = round(op$EOCRC/fvsm$SampleSummary$SampleSize[1], 4)
  op$LOCRC_freq = round(op$LOCRC/fvsm$SampleSummary$SampleSize[2], 4)
  genes = op$Hugo_Symbol[op$pval < 0.05] # adj
  genes
  
  write.table(op, file='./nonhyper_EOCRC_vs_LOCRC_freq.tsv', 
              quote=FALSE, row.names=FALSE, sep="\t")
  
  
  {
    
    
    df_nonhyper <- read_tsv('nonhyper_EOCRC_vs_LOCRC_freq.tsv')
    df_nonhyper <- df_nonhyper %>% 
      filter(!is.infinite(or) & !is.infinite(ci.low) & !is.infinite(ci.up))
    
    # set P_value and adjustPval if needed
    #p_value_threshold <- 0.05
    adjustPval_threshold <- 0.05
    
    filtered_df_nonhyper <- df_nonhyper %>% 
      filter(#pval <= p_value_threshold,
             adjPval <= adjustPval_threshold) %>%
      arrange(desc(EOCRC_freq)) 
    
    write.csv(filtered_df_nonhyper,
              file = paste(dir,'/7.nonhyper_differential_gene.csv',sep=''),
              row.names = F)
    
    filtered_df_nonhyper <-filtered_df_nonhyper %>% head(20)
    
    # add new labels
    plot <- filtered_df_nonhyper %>% 
      mutate(EOCRC_freq = format_frequency(EOCRC_freq),
             LOCRC_freq = format_frequency(LOCRC_freq),
             sample_size = paste0(EOCRC_freq, " vs. ", LOCRC_freq),
             OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
             adjPval = format_adjPval(adjPval)) %>% 
      mutate(across(everything(), as.character)) %>% 
      bind_rows(data.frame(Hugo_Symbol = "Gene", sample_size = "Frequency (EOCRC vs LOCRC)", adjPval = "adjPval", OR_CI = "OR (95% CI)")) %>% 
      mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
    
    
    p_mid <- filtered_df_nonhyper %>% 
      ggplot(aes(y = fct_rev(Hugo_Symbol))) +
      theme_classic() +
      geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
      geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
      labs(x = "Odds Ratio") +
      scale_x_continuous(trans = 'log10', breaks = c(0.2, 0.3, 1,2), 
                         labels = scales::number_format(accuracy = 0.01)) +
      coord_cartesian(ylim = c(1, 21),xlim = c(0.2, 2)) + 
      
      geom_vline(xintercept = 1, linetype = "dashed") +  
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(color = "black"))
    
    
    p_left <- plot %>% 
      ggplot(aes(y = Hugo_Symbol)) + 
      geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
      geom_text(aes(x = 1, label = sample_size), hjust = 0, size = 4,
                fontface = ifelse(plot$sample_size == "Frequency (EOCRC vs LOCRC)", "bold", "plain")) +
      theme_void() + coord_cartesian(xlim = c(0, 3))
    
   
    p_right <- plot %>% ggplot() +
      geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
                fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
      geom_text(aes(x = 1, y = Hugo_Symbol, label = adjPval), hjust = 0,
                fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
      theme_void() + coord_cartesian(xlim = c(0, 2))
    
    
    layout <- c(
      area(t = 0, l = 0, b = 30, r = 10), 
      area(t = 0, l = 11, b = 30, r = 17), 
      area(t = 0, l = 18.5, b = 30, r = 25))
    
  }
  
  
}



pdf(paste(dir,'/4_nonhyper_forestPlot.pdf',sep = ""), 
    width = 11, height = 6)

print(p_left + p_mid + p_right + plot_layout(design = layout))

dev.off()




filtered_df_overall = read.csv(paste(dir,'/3.overall_differential_gene.csv',sep=''))
filtered_df_nonhyper = read.csv(paste(dir,'/7.nonhyper_differential_gene.csv',sep=''))
filtered_df_hyper = read.csv(paste(dir,'/6.hyper_differential_gene.csv',sep=''))
insect_gene = intersect(filtered_df_hyper$Hugo_Symbol,
                        filtered_df_nonhyper$Hugo_Symbol)
insect_gene = intersect(insect_gene,filtered_df_overall$Hugo_Symbol)


bubble_table = data.frame(Hugo_Symbol=character(), or=numeric(), 
                          Fold=numeric(), TMB=character(), Country=character())

total_bubble_table <- data.frame(Hugo_Symbol=character(), or=numeric(), 
                                 Fold=numeric(), TMB=character(), Country=character())

generate_bubble_table <- function(kokka){
  
  country_sample <- rt@clinical.data$Tumor_Sample_Barcode[rt@clinical.data$Country == kokka]
 
   if (length(intersect(clin.EOCRC@clinical.data$Tumor_Sample_Barcode, country_sample)) > 0 &
      length(intersect(clin.LOCRC@clinical.data$Tumor_Sample_Barcode, country_sample)) > 0){
  rt_country_EO <- subsetMaf(clin.EOCRC,tsb=country_sample)
  rt_country_LO <- subsetMaf(clin.LOCRC,tsb=country_sample)
  country_comp <- mafCompare(rt_country_EO,rt_country_LO) $results
  country_comp$Fold = country_comp$M1 / country_comp$M2
  country_comp$TMB = 'Overall'
  country_comp$Country = kokka
  country_comp = country_comp[country_comp$Hugo_Symbol %in% insect_gene]
      
 # for (gene in insect_gene){
  #  if (!(gene %in% country_comp$Hugo_Symbol)){
   #   test.data =   data.frame(Hugo_Symbol=gene, or=1,Fold=NA,
    #                           TMB='Overall', 
     #                          Country=kokka)
      
     # country_comp = bind_rows(test.data,country_comp)
    #}
    
    
 # }
}else{
  
  country_comp <- data.frame(Hugo_Symbol=character(), or=numeric(), 
                                   Fold=numeric(), TMB=character(), Country=character())
}


  if (length(intersect(clin.EOCRC_hyper@clinical.data$Tumor_Sample_Barcode, country_sample)) > 0 &
      length(intersect(clin.LOCRC_hyper@clinical.data$Tumor_Sample_Barcode, country_sample)) > 0) {
  rt_country_hyper_EO <- subsetMaf(clin.EOCRC_hyper,tsb=country_sample)
  rt_country_hyper_LO <- subsetMaf(clin.LOCRC_hyper,tsb=country_sample)
  country_hyper_comp <- mafCompare(rt_country_hyper_EO,rt_country_hyper_LO) $results
  country_hyper_comp$Fold = country_hyper_comp$M1 / country_hyper_comp$M2
  country_hyper_comp$TMB = 'Hypermutated'
  country_hyper_comp$Country = kokka
  country_hyper_comp = country_hyper_comp[country_hyper_comp$Hugo_Symbol %in% insect_gene]
 
  #for (gene in insect_gene){
  #  if (!(gene %in% country_hyper_comp$Hugo_Symbol)){
 #     test.data =   data.frame(Hugo_Symbol=gene, or=1,Fold=NA,
   #                            TMB='Hypermutated', 
    #                           Country=kokka)
      
     # country_hyper_comp = bind_rows(test.data,country_hyper_comp)
    #}
#  }
  
   }else{
    
    country_hyper_comp <- data.frame(Hugo_Symbol=character(), or=numeric(), 
                                     Fold=numeric(), TMB=character(), Country=character())
    
   } 
  
  if (length(intersect(clin.EOCRC_nonhyper@clinical.data$Tumor_Sample_Barcode, country_sample)) > 0 &
      length(intersect(clin.LOCRC_nonhyper@clinical.data$Tumor_Sample_Barcode, country_sample)) > 0) {
  rt_country_nonhyper_EO <- subsetMaf(clin.EOCRC_nonhyper,tsb=country_sample)
  rt_country_nonhyper_LO <- subsetMaf(clin.LOCRC_nonhyper,tsb=country_sample)
  country_nonhyper_comp <- mafCompare(rt_country_nonhyper_EO,rt_country_nonhyper_LO) $results
  country_nonhyper_comp$Fold = country_nonhyper_comp$M1 / country_nonhyper_comp$M2
  country_nonhyper_comp$TMB = 'Nonhypermutated'
  country_nonhyper_comp$Country = kokka
  country_nonhyper_comp = country_nonhyper_comp[country_nonhyper_comp$Hugo_Symbol %in% insect_gene]
  
 # for (gene in insect_gene){
  #  if (!(gene %in% country_nonhyper_comp$Hugo_Symbol)){
   #   test.data =   data.frame(Hugo_Symbol=gene, or=1,Fold=NA,
    #                           TMB='Nonhypermutated', 
     #                          Country=kokka)
      
    #  country_nonhyper_comp = bind_rows(test.data,country_nonhyper_comp)
 #   }
#  }

  
  
  } else {
    country_nonhyper_comp <- data.frame(Hugo_Symbol=character(), or=numeric(), 
                                        Fold=numeric(), TMB=character(), Country=character())
  }
  
  
  bubble_table = bind_rows(bubble_table,country_nonhyper_comp,
                       country_hyper_comp,country_comp)
  
  
  return(bubble_table)
  
}

                         
                         
for (kokka in unique(rt@clinical.data$Country) ){
  country_bubble_table <- generate_bubble_table(kokka)
  total_bubble_table <- bind_rows(total_bubble_table, country_bubble_table)
    
}

names(total_bubble_table)
total_bubble_table = total_bubble_table[,c(1:5)]

#total_bubble_table$Fold = 2^(total_bubble_table$Fold)

names(total_bubble_table)[names(total_bubble_table) =='Hugo_Symbol'] = 'Gene'
names(total_bubble_table)[names(total_bubble_table) =='TMB'] = 'TMB_Category'
names(total_bubble_table)[names(total_bubble_table) =='or'] = 'OR'

#total_bubble_table = total_bubble_table %>%
 # filter(!(Country == "Nigeria" & TMB_Category == "Nonhypermutated"))



p1 <- ggplot(data = total_bubble_table, 
            aes(x = Gene, y = `TMB_Category`)) +
  geom_point(aes(size = OR, color = Fold), alpha = 0.7) +
  facet_grid(Country ~ ., scales = "free_y", space = "free") +  
  scale_y_discrete(position = "left") +  
  scale_size_continuous(limits = c(0,3),breaks = c(0,1,2,3)) +  
  scale_color_gradient2(midpoint = 1,high='red',mid = 'white',low  = 'blue' 
                          ) +  
  theme_few() +
  theme(
    strip.background = element_blank(), 
    strip.text.y.right = element_text(size = 10, angle = 0),
    panel.spacing = unit(0, "lines"),   
    panel.border = element_rect(fill = NA, color = "black", size = 0.5) ,
    #axis.text.x   = element_blank(),
    axis.text.y   = element_text(color = "black", size = 10, angle = 0)
  )


p1
ggsave(paste(dir,'/7.interser_gene_among_countries.pdf',sep = ''),
       width = 6,height = 6,dpi = 300)

{
  
  #p1 = forestPlot(mafCompareRes=fvsm, pVal=0.05, color=c("maroon", "royalblue"), geneFontSize=0.8)
  pdf(paste(dir,'/4_nonhyper_forestPlot_mutation_p0.001_fdr0.001.pdf',sep = ""), width = 12, height = 10)
  forestPlot(mafCompareRes = fvsm, pVal=0.001,  fdr = 0.001,
             color=c("maroon", "royalblue"), 
             titleSize = 2,
             geneFontSize = 1,
             lineWidth = 1.5
  )
  dev.off()
  
  se_gene = getGeneSummary(clin.EOCRC)[1:30]$Hugo_Symbol
  png(paste(dir,'/5_nonhyper_coOncoplot_top30.png',sep = ""),
      width = 1500, height =1000)
  #pdf(paste(dir,'/5_nonhyper_coOncoplot_top30.pdf',sep = ""),width = 15, height = 10)
  coOncoplot(m1=clin.EOCRC, m2=clin.LOCRC, 
             m1Name="EOCRC", m2Name="LOCRC", 
             genes = se_gene)
  dev.off()
  
  png(paste(dir,'/6_nonhyper_coOncoplot_diff_top30.png',sep = ""),
      width = 1500, height =1000)
  #pdf(paste(dir,'/6_nonhyper_coOncoplot_diff_top30.pdf',sep = ""), width = 15, height = 10)
  coOncoplot(m1=clin.EOCRC, m2=clin.LOCRC, 
             m1Name="EOCRC", m2Name="LOCRC", 
             genes = genes) # [1:30]
  dev.off()
}
