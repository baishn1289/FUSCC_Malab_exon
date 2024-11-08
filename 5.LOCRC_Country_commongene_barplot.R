library(ggsci)
library(ggplot2)
library(tidyverse)
library(treeio)
library(ape)
library(ggnewscale)
library(tidyr)
library(ggthemes)



  
  
  countries_LO <- unique(rt@clinical.data$Country)[1:7]
  
  Gene_list_country_LO <- left_join(clin.LOCRC@data[,.(Hugo_Symbol,Tumor_Sample_Barcode)],
                                    clin.LOCRC@clinical.data[,.(Tumor_Sample_Barcode,Country)],
                                    by = 'Tumor_Sample_Barcode') %>% 
    group_by(Country) %>% 
    reframe(
      Gene_list = list(Hugo_Symbol)
    )
  
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
  
  
  for (i in common_genes) {
    generate_country_frequency(countries = countries_LO, 
                               maf = clin.LOCRC, jiyin = i)
  }
  
  
  
  # gene_symbols <- c()
  # generate_country_sample(countries = countries , cl = cl, rt = rt, XOCRC = 'LOCRC')
  
  
  # for (country in countries) {
  #   gene_symbols <- c(gene_symbols, get(paste(country, "LOCRC", sep = "."))$Hugo_Symbol)
  # }
  # uniqueGenes <- genes[genes %in% unique(gene_symbols)]
  # 
  # generate_country_frequency(countries,'LOCRC',uniqueGenes = uniqueGenes)
  # 
  # mutation_table <- data.frame(gene = uniqueGenes,
  #                              Canada = Canada.LOCRC_commongene$mutfreq[match(uniqueGenes, Canada.LOCRC_commongene$Hugo_Symbol)],
  #                              Netherlands = Netherlands.LOCRC_commongene$mutfreq[match(uniqueGenes, Netherlands.LOCRC_commongene$Hugo_Symbol)],
  #                              Nigeria = Nigeria.LOCRC_commongene$mutfreq[match(uniqueGenes, Nigeria.LOCRC_commongene$Hugo_Symbol)],
  #                              Spain = Spain.LOCRC_commongene$mutfreq[match(uniqueGenes, Spain.LOCRC_commongene$Hugo_Symbol)],
  #                              US = US.LOCRC_commongene$mutfreq[match(uniqueGenes, US.LOCRC_commongene$Hugo_Symbol)],
  #                              China = China.LOCRC_commongene$mutfreq[match(uniqueGenes, China.LOCRC_commongene$Hugo_Symbol)],
  #                              France = France.LOCRC_commongene$mutfreq[match(uniqueGenes, France.LOCRC_commongene$Hugo_Symbol)])
  # 
  # mutation_table[is.na(mutation_table)] = 0
  # mutation_table <- na.omit(mutation_table)

  mutation_table_LO <- mutation_table %>% 
    mutate(Korea = 0) %>% 
    mutate(Age_class = 'LOCRC')
  
  
  
  #EO和LO相同基因单独做条形图，后续再做条形图进行筛选
  names(mutation_table_LO)
  names(mutation_table_EO)
  common_EO_LO_gene <- intersect(mutation_table_EO$gene,
                                 mutation_table_LO$gene)
  
  #不同的基因作分面，EO和LO做变量，Country做分类横坐标
 
  common_EO_LO_df <- rbind(mutation_table_EO[gene%in%common_EO_LO_gene,],
                           mutation_table_LO[gene%in%common_EO_LO_gene,]) %>% 
    pivot_longer(cols = Canada:Korea,  # 将国家列名替换为实际列名范围
                 names_to = "Country", # 将这些列名变成 Country 列
                 values_to = "Frequency")
 
  
  palette <- c('#1f78b4','#fc8d62')
  
  
  plot_common_EO_LO <- common_EO_LO_df %>% 
    ggplot(aes(Frequency, Country, fill = Age_class)) +
    geom_col(position = position_dodge(1), width = 0.6) + # 增大间距并调整条形宽度
    labs(y = NULL) +
    facet_wrap('gene', scales = "free", nrow  = 3) + # 将分面标题移到左侧
    theme_few() +
    # scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_fill_manual(values = palette) +
    theme(
      axis.text.y = element_text(color = "black", size = 8), # y轴文本显示
      # axis.text.x = element_text(color = "black", size = 8),
      axis.text.x   = element_text(color = "black", size = 8, hjust=1,angle = 45),
      axis.title = element_text(color = "black", size = 14, face = "bold"),
      axis.ticks.y = element_blank(), 
      legend.position = "right", 
      plot.background = element_blank(), # 去掉背景
      panel.background = element_blank(), # 去掉面板背景
      strip.background = element_blank(), # 去掉分面标题背景
      strip.placement = "outside", # 将分面标题移到左侧外部
      strip.text.y.left = element_text(angle = 0) # 调整分面标题的方向为水平
    )+
    coord_flip()
  
  pdf(paste(dir,"/4_EOCRC_LOCRC_common_genes_mainplot.pdf", sep = ""), width = 10, height = 6)
  plot_common_EO_LO
  dev.off()
  
  
  
  
  ######################
  
  
  
  
  # common_EO_LO_df$Max_Frequency <- apply(common_EO_LO_df[, 2:9], 1, max)
  # 
  # common_EO_LO_df_higherfreq <- common_EO_LO_df %>% 
  #   group_by(gene) %>% 
  #   filter(all(Max_Frequency > 0.14)) %>% 
  #   ungroup() %>% 
  #   select(-Max_Frequency) %>% 
  #   pivot_longer(cols = Canada:Korea,  # 将国家列名替换为实际列名范围
  #                          names_to = "Country", # 将这些列名变成 Country 列
  #                          values_to = "Frequency")
  # 

  
  # palette <- c('#deebf7','#3182bd')
  
  # plot_common_EO_LO_df_higherfreq <- common_EO_LO_df_higherfreq %>% 
  #   ggplot(aes(Frequency, Country, fill = Age_class)) +
  #   geom_col(position = position_dodge(1), width = 0.6) + # 增大间距并调整条形宽度
  #   labs(y = NULL) +
  #   facet_wrap('gene', scales = "free", nrow  = 1) + # 将分面标题移到左侧
  #   theme_few() +
  #   # scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  #   scale_fill_manual(values = palette) +
  #   theme(
  #     axis.text.y = element_text(color = "black", size = 8), # y轴文本显示
  #     # axis.text.x = element_text(color = "black", size = 8),
  #     axis.text.x   = element_text(color = "black", size = 8, hjust=1,angle = 45),
  #     axis.title = element_text(color = "black", size = 14, face = "bold"),
  #     axis.ticks.y = element_blank(), 
  #     legend.position = "right", 
  #     plot.background = element_blank(), # 去掉背景
  #     panel.background = element_blank(), # 去掉面板背景
  #     strip.background = element_blank(), # 去掉分面标题背景
  #     strip.placement = "outside", # 将分面标题移到左侧外部
  #     strip.text.y.left = element_text(angle = 0) # 调整分面标题的方向为水平
  #   )+
  #   coord_flip() 
  # 
  # plot_common_EO_LO_df_higherfreq
  # 
  # common_EO_LO_df_lowerfreq <- common_EO_LO_df %>% 
  #   filter(!(gene%in% common_EO_LO_df_higherfreq$gene)) %>% 
  #   select(-Max_Frequency) %>% 
  #   pivot_longer(cols = Canada:Korea,  # 将国家列名替换为实际列名范围
  #                names_to = "Country", # 将这些列名变成 Country 列
  #                values_to = "Frequency")
  # 
  # plot_common_EO_LO_df_lowerfreq <- common_EO_LO_df_lowerfreq %>% 
  #   ggplot(aes(Frequency, Country, fill = Age_class)) +
  #   geom_col(position = position_dodge(1), width = 0.6) + # 增大间距并调整条形宽度
  #   labs(y = "Country") +
  #   facet_wrap('gene', scales = "free", nrow  = 1) + # 将分面标题移到左侧
  #   theme_few() + 
  #   # scale_x_continuous(expand = c(0, 0), limits = c(0, 0.35)) +
  #   scale_fill_manual(values = palette) +
  #   theme(
  #     axis.text.y = element_text(color = "black", size = 8), # y轴文本显示
  #     # axis.text.x = element_text(color = "black", size = 8),
  #     axis.text.x   = element_text(color = "black", size = 8, hjust=1,angle = 45),
  #     axis.title = element_text(color = "black", size = 14, face = "bold"),
  #     axis.ticks.y = element_blank(),
  #     legend.position = "right",
  #     plot.background = element_blank(),
  #     panel.background = element_blank(), 
  #     strip.background = element_blank(), 
  #     strip.placement = "outside", # 将分面标题移到左侧外部
  #     strip.text.y.left = element_text(angle = 0) # 调整分面标题的方向为水平
  #   )+
  #   coord_flip() 
  # 
  # 
  # plot_common_EO_LO_df_lowerfreq
  # 
  # 
  # library(cowplot)
  # pdf(paste(dir,"/4_EOCRC_LOCRC_common_genes_mainplot.pdf", sep = ""), width = 18, height = 5)
  # plot_grid( plot_common_EO_LO_df_higherfreq,
  #            plot_common_EO_LO_df_lowerfreq, 
  #              
  #              nrow = 2,rel_heights = c(2, 2.25),
  #                       rel_widths = c(1,1))
  #   dev.off()
  
    {

    #supplementary plot
    
    # palette <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')
    
    palette <- c('Canada' = '#b2182b',
                 'China' = '#d6604d',
                 'France' = '#f4a582',
                 
                 'Korea' = '#fddbc7',
                 'Netherlands' = '#d1e5f0',
                 'Nigeria' = '#92c5de',
                 'Spain' = '#4393c3',
                 'US'='#2166ac') 
    
    
   
    #supplementary plot_LO
    mutation_table_LO_sup <- mutation_table_LO %>% 
      filter(!(gene %in% common_EO_LO_df$gene)) %>% 
      select(-Age_class)
   
  
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
  ggsave(paste(dir,"/4_LOCRC_common_genes_supplementary.pdf",sep = ""), width = 24, height = 10, dpi = 300)

  
  #supplementary plot_EO
  
  mutation_table_EO_sup <- mutation_table_EO %>% 
    filter(!(gene %in% common_EO_LO_df$gene)) %>% 
    select(-Age_class)
  
  
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
ggsave(paste(dir,"/4_EOCRC_common_genes_supplementary.pdf",sep = ""), width = 11, height = 4, dpi = 300)
  
  
# 
#   library(gridExtra)
#   pdf(paste(dir,"/4_EOCRC_LOCRC_common_genes_supplementary.pdf", sep = ""), width = 35, height = 22)
#    plot_grid(Common_country_EO, Common_country_LO, 
#                  ncol = 1, rel_heights = c(1, 3),rel_widths = c(1,8))
#   
#   dev.off()

    }
  