library(ggsci)
library(ggplot2)
library(tidyverse)
library(treeio)
library(ape)
library(ggnewscale)
library(MetBrewer)
library(tidyr)
library(ggthemes)


{#LOCRC_Country
  
  # Extract all unique Hugo_Symbol
  countries = countries[1:7]
  gene_symbols <- c()
  generate_country_sample(countries = countries , cl = cl, rt = rt, XOCRC = 'LOCRC')
  
  
  for (country in countries) {
    gene_symbols <- c(gene_symbols, get(paste(country, "LOCRC", sep = "."))$Hugo_Symbol)
  }
  uniqueGenes <- genes[genes %in% unique(gene_symbols)]
  
  generate_country_frequency(countries,'LOCRC',uniqueGenes = uniqueGenes)
  
  mutation_table <- data.frame(gene = uniqueGenes,
                               Canada = Canada.LOCRC_commongene$mutfreq[match(uniqueGenes, Canada.LOCRC_commongene$Hugo_Symbol)],
                               Netherlands = Netherlands.LOCRC_commongene$mutfreq[match(uniqueGenes, Netherlands.LOCRC_commongene$Hugo_Symbol)],
                               Nigeria = Nigeria.LOCRC_commongene$mutfreq[match(uniqueGenes, Nigeria.LOCRC_commongene$Hugo_Symbol)],
                               Spain = Spain.LOCRC_commongene$mutfreq[match(uniqueGenes, Spain.LOCRC_commongene$Hugo_Symbol)],
                               US = US.LOCRC_commongene$mutfreq[match(uniqueGenes, US.LOCRC_commongene$Hugo_Symbol)],
                               China = China.LOCRC_commongene$mutfreq[match(uniqueGenes, China.LOCRC_commongene$Hugo_Symbol)],
                               France = France.LOCRC_commongene$mutfreq[match(uniqueGenes, France.LOCRC_commongene$Hugo_Symbol)])
  
  # mutation_table[is.na(mutation_table)] = 0
  mutation_table <- na.omit(mutation_table)

  
  plotdata <- data.frame(t(mutation_table))
  colnames(plotdata) <- plotdata[1,]
  plotdata = plotdata[-1,]
  plotdata$country = rownames(plotdata)
  
  palette <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')

  palette <- c('Canada' = '#b2182b',
               'China' = '#d6604d',
               'France' = '#f4a582',
               
               #'Korea' = '#fddbc7',
               'Netherlands' = '#d1e5f0',
               'Nigeria' = '#92c5de',
               'Spain' = '#4393c3',
               'US'='#2166ac') 
  
  
  mergedata <- reshape2::melt(plotdata[,], id = 'country') 
  mergedata$value = round(as.numeric(mergedata$value), 3)
} 
  

    ggplot(mergedata, aes(x = country, y = value, fill = country)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(~ variable, scales = "free", ncol = 3) +
    theme_few() +
    scale_fill_manual(values = palette) + 
   theme(
      plot.title    = element_text(color = "black", size   = 16, hjust = 0.5),
      plot.subtitle = element_text(color = "black", size   = 16,hjust = 0.5),
      plot.caption  = element_text(color = "black", size   = 16,face = "italic", hjust = 1),
      axis.text.x   = element_text(color = "black", size = 16, hjust=1,angle = 45),
      axis.text.y   = element_text(color = "black", size = 16, angle = 0),
      axis.title.x  = element_text(color = "black", size = 16, angle = 0),
      axis.title.y  = element_text(color = "black", size = 16, angle = 90),
      legend.title  = element_text(color = "black", size  = 16),
      legend.text   = element_text(color = "black", size   = 16),
      strip.text = element_text(size = 14),
      axis.line.y = element_line(color = "black", linetype = "solid"), 
      axis.line.x = element_line (color = "black",linetype = "solid"), 
      panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA) )
 

  ggsave(paste(dir,"/4_LOCRC_common_genes.pdf",sep = ""), width = 10, height = 8, dpi = 300)
