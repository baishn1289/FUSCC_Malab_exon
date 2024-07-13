library(ggsci)
library(ggplot2)
library(tidyverse)
library(treeio)
library(ape)
library(ggnewscale)
library(MetBrewer)
library(tidyr)
library(ggthemes)

# Load unique countries
countries <- unique(rt@clinical.data$Country)

# Function to generate Country sample and corresponding gene summary
generate_country_sample <- function(countries, cl, rt, XOCRC) {
  for (kokka in countries) {
    if (XOCRC == 'LOCRC') {
      country_sample <- subset(cl, Country == kokka & Age > 50)$Tumor_Sample_Barcode
      country_LOCRC <- getGeneSummary(subsetMaf(maf = rt, tsb = country_sample, isTCGA = FALSE))
      
      assign(paste(kokka, "sample", sep = "."), country_sample, envir = .GlobalEnv)
      assign(paste(kokka, "LOCRC", sep = "."), country_LOCRC, envir = .GlobalEnv)
    } else {
      country_sample <- subset(cl, Country == kokka & Age <= 50)$Tumor_Sample_Barcode
      country_EOCRC <- getGeneSummary(subsetMaf(maf = rt, tsb = country_sample, isTCGA = FALSE))
      
      assign(paste(kokka, "sample", sep = "."), country_sample, envir = .GlobalEnv)
      assign(paste(kokka, "EOCRC", sep = "."), country_EOCRC, envir = .GlobalEnv)
    }
  }
}

# Function to generate Country frequency for common genes
generate_country_frequency <- function(countries, XOCRC, uniqueGenes) {
  for (kokka in countries) {
    if (XOCRC == 'LOCRC') {
      country_sample <- get(paste(kokka, "sample", sep = "."))
      country_LOCRC <- get(paste(kokka, "LOCRC", sep = "."))
      
      country_LOCRC_commongene <- subset(country_LOCRC, Hugo_Symbol %in% uniqueGenes)
      country_LOCRC_commongene$mutfreq <- country_LOCRC_commongene$MutatedSamples / length(country_sample)
      
      assign(paste(kokka, "LOCRC_commongene", sep = "."), country_LOCRC_commongene, envir = .GlobalEnv)
      
    } else {
      country_sample <- get(paste(kokka, "sample", sep = "."))
      country_EOCRC <- get(paste(kokka, "EOCRC", sep = "."))
      
      country_EOCRC_commongene <- subset(country_EOCRC, Hugo_Symbol %in% uniqueGenes)
      country_EOCRC_commongene$mutfreq <- country_EOCRC_commongene$MutatedSamples / length(country_sample)
      
      assign(paste(kokka, "EOCRC_commongene", sep = "."), country_EOCRC_commongene, envir = .GlobalEnv)
    }
  }
}

# Extract all unique Hugo_Symbol
gene_symbols <- c()
generate_country_sample(countries = countries , cl = cl, rt = rt, XOCRC = 'EOCRC')

  for (Country in countries) {
    gene_symbols <- c(gene_symbols, get(paste(Country, "EOCRC", sep = "."))$Hugo_Symbol)
  }
  uniqueGenes <- genes[genes %in% unique(gene_symbols)]
  
  generate_country_frequency(countries,'EOCRC',uniqueGenes = uniqueGenes)

  
  mutation_table <- data.frame(gene = uniqueGenes,
                               Canada = Canada.EOCRC_commongene$mutfreq[match(uniqueGenes, Canada.EOCRC_commongene$Hugo_Symbol)],
                               Netherlands = Netherlands.EOCRC_commongene$mutfreq[match(uniqueGenes, Netherlands.EOCRC_commongene$Hugo_Symbol)],
                               Nigeria = Nigeria.EOCRC_commongene$mutfreq[match(uniqueGenes, Nigeria.EOCRC_commongene$Hugo_Symbol)],
                               Spain = Spain.EOCRC_commongene$mutfreq[match(uniqueGenes, Spain.EOCRC_commongene$Hugo_Symbol)],
                               US = US.EOCRC_commongene$mutfreq[match(uniqueGenes, US.EOCRC_commongene$Hugo_Symbol)],
                               China = China.EOCRC_commongene$mutfreq[match(uniqueGenes, China.EOCRC_commongene$Hugo_Symbol)],
                               France = France.EOCRC_commongene$mutfreq[match(uniqueGenes, France.EOCRC_commongene$Hugo_Symbol)],
                               Korea = Korea.EOCRC_commongene$mutfreq[match(uniqueGenes, Korea.EOCRC_commongene$Hugo_Symbol)])
  
 
  # mutation_table[is.na(mutation_table)] = 0
  mutation_table <- na.omit(mutation_table)

  
  
  plotdata <- data.frame(t(mutation_table))
  colnames(plotdata) <- plotdata[1,]
  plotdata = plotdata[-1,]
  plotdata$Country = rownames(plotdata)
  
  palette <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')
  
  palette <- c('Canada' = '#b2182b',
               'China' = '#d6604d',
               'France' = '#f4a582',
               
               'Korea' = '#fddbc7',
               'Netherlands' = '#d1e5f0',
               'Nigeria' = '#92c5de',
               'Spain' = '#4393c3',
               'US'='#2166ac') 
  
  mergedata <- reshape2::melt(plotdata[,], id = 'Country') #c(1:10, ncol(plotdata))
  mergedata$value = round(as.numeric(mergedata$value), 3)

  ggplot(mergedata, aes(x = Country, y = value, fill = Country)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(~ variable, scales = "free", ncol = 2) + 
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
    
  ggsave(paste(dir,"/4_EOCRC_common_genes.pdf",sep = ""), width = 10, height = 8, dpi = 300) 
  