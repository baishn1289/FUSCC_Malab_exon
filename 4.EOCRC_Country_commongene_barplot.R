library(ggsci)
library(ggplot2)
library(tidyverse)
library(treeio)
library(ape)
library(ggnewscale)
library(tidyr)
library(ggthemes)


#国家间共同基因代码逻辑梳理:
#先提取出各个国家共有基因
#然后计算这些共有基因在各个国家中的频率

# Load unique countries
countries_EO <- unique(rt@clinical.data$Country)

# kokka <- 'US'
# country_sample <- subset(cl, Country == kokka & Age > 50)$Tumor_Sample_Barcode
# country_LOCRC <- getGeneSummary(subsetMaf(maf = rt, tsb = country_sample, isTCGA = FALSE))
# 
# country_LOCRC_commongene <- subset(country_LOCRC, Hugo_Symbol %in% uniqueGenes)
# country_LOCRC_commongene$mutfreq <- country_LOCRC_commongene$MutatedSamples / length(country_sample)
# 
# 
# gene_symbols <- country_LOCRC$Hugo_Symbol
# 
# generate_country_frequency(countries,'EOCRC',uniqueGenes = uniqueGenes)

#提取出青年队列中各个国家所拥有的基因
Gene_list_country_EO <- left_join(clin.EOCRC@data[,.(Hugo_Symbol,Tumor_Sample_Barcode)],
                     clin.EOCRC@clinical.data[,.(Tumor_Sample_Barcode,Country)],
          by = 'Tumor_Sample_Barcode') %>% 
  group_by(Country) %>% 
  reframe(
    Gene_list = list(Hugo_Symbol)
  )

#提取出各个国家的共有基因
gene_lists <- lapply(Gene_list_country_EO$Gene_list, unlist)
common_genes <- Reduce(intersect, gene_lists)

mutation_table <- data.table(
  gene = common_genes,
  Canada = NA_real_,
  China = NA_real_,
  US = NA_real_,
  Netherlands = NA_real_,
  Nigeria = NA_real_,
  France = NA_real_,
  Spain = NA_real_,
  Korea = NA_real_
)


# 
# maf = clin.EOCRC
# jiyin = 'TP53'
# guojia = 'France'



generate_country_frequency <- function(countries, maf,jiyin){
  
  #这里设计函数同时适用于青年和老年国家队列，用clinical.EOCRC和clinical.LOCRC代替maf即可完成队列的选择
  
  #提取出有该基因变异的所有国家的样本
  Alter_Gene = maf@data[Hugo_Symbol == jiyin,] %>% 
    filter(!duplicated(Tumor_Sample_Barcode))
  
  #给突变样本附上国家信息
  mut_samples <- merge(Alter_Gene[, c('Tumor_Sample_Barcode', 'Hugo_Symbol')],
                                           maf@clinical.data[, c('Tumor_Sample_Barcode', "Country")],
                                           by = 'Tumor_Sample_Barcode')
  
  #提取出所有包含该基因的panels
  Gene_panels <- unique(panel_hugo_symbol[Hugo_Symbol == jiyin]$Panel)
  
  
  for (guojia in countries) {
     
      
    #某个基因在某个国家的突变频率计算：该基因在该国家的突变样本数/该国家中测量该基因的panel的样本数量
    
    #根据国家提取出该国家的所有该基因突变样本
    mut_country_samples <- mut_samples %>% 
      filter(Country == guojia)  %>% 
      summarise(unique_samples_count = n_distinct(Tumor_Sample_Barcode)) %>% 
      pull(unique_samples_count)
    
    #提取出该国家中包含该突变基因的panels，并计算出panels中的样本数量
    panels_country_samples <- maf@clinical.data %>% 
      filter(Panel %in% Gene_panels,
             Country == guojia) %>% 
      summarise(unique_samples_count = n_distinct(Tumor_Sample_Barcode)) %>% 
      pull(unique_samples_count)
    
    #计算该国该基因突变频率，并计入表格中，便于后续绘图
    mutation_table[gene == jiyin, (guojia) := mut_country_samples / panels_country_samples]
    
  }
  mutation_table  <<- mutation_table
}


for (i in common_genes) {
  generate_country_frequency(countries = countries_EO, 
                             maf = clin.EOCRC, jiyin = i)
}



# # Function to generate Country sample and corresponding gene summary
# generate_country_sample <- function(countries, cl, rt, XOCRC) {
#   for (kokka in countries) {
#     if (XOCRC == 'LOCRC') {
#       country_sample <- subset(cl, Country == kokka & Age > 50)$Tumor_Sample_Barcode
#       country_LOCRC <- getGeneSummary(subsetMaf(maf = rt, tsb = country_sample, isTCGA = FALSE))
#
#       assign(paste(kokka, "sample", sep = "."), country_sample, envir = .GlobalEnv)
#       assign(paste(kokka, "LOCRC", sep = "."), country_LOCRC, envir = .GlobalEnv)
#     } else {
#       country_sample <- subset(cl, Country == kokka & Age <= 50)$Tumor_Sample_Barcode
#       country_EOCRC <- getGeneSummary(subsetMaf(maf = rt, tsb = country_sample, isTCGA = FALSE))
#
#       assign(paste(kokka, "sample", sep = "."), country_sample, envir = .GlobalEnv)
#       assign(paste(kokka, "EOCRC", sep = "."), country_EOCRC, envir = .GlobalEnv)
#     }
#   }
# }
#
# # Function to generate Country frequency for common genes
# generate_country_frequency <- function(countries, XOCRC, uniqueGenes) {
#   for (kokka in countries) {
#     if (XOCRC == 'LOCRC') {
#       country_sample <- get(paste(kokka, "sample", sep = "."))
#       country_LOCRC <- get(paste(kokka, "LOCRC", sep = "."))
#
#       country_LOCRC_commongene <- subset(country_LOCRC, Hugo_Symbol %in% uniqueGenes)
#       #某个基因在某个国家的突变频率计算：该基因在该国家的突变患者数/该国家中测量该基因的panel的患者数量
#       country_LOCRC_commongene$mutfreq <- country_LOCRC_commongene$MutatedSamples / length(country_sample)
#
#       assign(paste(kokka, "LOCRC_commongene", sep = "."), country_LOCRC_commongene, envir = .GlobalEnv)
#
#     } else {
#       country_sample <- get(paste(kokka, "sample", sep = "."))
#       country_EOCRC <- get(paste(kokka, "EOCRC", sep = "."))
#
#       country_EOCRC_commongene <- subset(country_EOCRC, Hugo_Symbol %in% uniqueGenes)
#       country_EOCRC_commongene$mutfreq <- country_EOCRC_commongene$MutatedSamples / length(country_sample)
#
#       assign(paste(kokka, "EOCRC_commongene", sep = "."), country_EOCRC_commongene, envir = .GlobalEnv)
#     }
#   }
# }

#
# # Extract all unique Hugo_Symbol
# gene_symbols <- c()
# generate_country_sample(countries = countries , cl = cl, rt = rt, XOCRC = 'EOCRC')
#
# #在 generate_country_sample 函数执行后，遍历每个国家的 EOCRC 数据，提取所有Hugo_Symbol并合并到 gene_symbols 列表中
#   for (Country in countries) {
#     gene_symbols <- c(gene_symbols, get(paste(Country, "EOCRC", sep = "."))$Hugo_Symbol)
#   }
#
#
#   uniqueGenes <- genes[genes %in% unique(gene_symbols)]
#
#   generate_country_frequency(countries,'EOCRC',uniqueGenes = uniqueGenes)
#
#
#   mutation_table <- data.frame(gene = uniqueGenes,
#                                Canada = Canada.EOCRC_commongene$mutfreq[match(uniqueGenes, Canada.EOCRC_commongene$Hugo_Symbol)],
#                                Netherlands = Netherlands.EOCRC_commongene$mutfreq[match(uniqueGenes, Netherlands.EOCRC_commongene$Hugo_Symbol)],
#                                Nigeria = Nigeria.EOCRC_commongene$mutfreq[match(uniqueGenes, Nigeria.EOCRC_commongene$Hugo_Symbol)],
#                                Spain = Spain.EOCRC_commongene$mutfreq[match(uniqueGenes, Spain.EOCRC_commongene$Hugo_Symbol)],
#                                US = US.EOCRC_commongene$mutfreq[match(uniqueGenes, US.EOCRC_commongene$Hugo_Symbol)],
#                                China = China.EOCRC_commongene$mutfreq[match(uniqueGenes, China.EOCRC_commongene$Hugo_Symbol)],
#                                France = France.EOCRC_commongene$mutfreq[match(uniqueGenes, France.EOCRC_commongene$Hugo_Symbol)],
#                                Korea = Korea.EOCRC_commongene$mutfreq[match(uniqueGenes, Korea.EOCRC_commongene$Hugo_Symbol)])


  # # mutation_table[is.na(mutation_table)] = 0
  # mutation_table <- na.omit(mutation_table)


  mutation_table_EO <- mutation_table %>%
    mutate(Age_class = 'EOCRC')

  # plotdata_EO <- data.frame(t(mutation_table))
  # 
  # 
  # 
  # colnames(plotdata_EO) <- plotdata_EO[1,]
  # plotdata_EO = plotdata_EO[-1,]
  # plotdata_EO$Country = rownames(plotdata_EO)
# 
# 
#   palette <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')
# 
#   palette <- c('Canada' = '#b2182b',
#                'China' = '#d6604d',
#                'France' = '#f4a582',
# 
#                'Korea' = '#fddbc7',
#                'Netherlands' = '#d1e5f0',
#                'Nigeria' = '#92c5de',
#                'Spain' = '#4393c3',
#                'US'='#2166ac')
# 
#   mergedata <- reshape2::melt(plotdata_EO[,], id = 'Country') #c(1:10, ncol(plotdata_EO))
#   mergedata$value = round(as.numeric(mergedata$value), 3)
# 
# common_EO = ggplot(mergedata, aes(x = Country, y = value, fill = Country)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     facet_wrap(~ variable, scales = "free", nrow = 2) +
#     theme_few() +
#     scale_fill_manual(values = palette) +
#     labs(title = "Common genes among countries",
#          subtitle = "EOCRC" )+
#   theme(
#     plot.title    = element_text(color = "black", size   = 20, hjust = 0.5),
#     plot.subtitle = element_text(color = "black", size   = 16,hjust = 0.5),
#     plot.caption  = element_text(color = "black", size   = 16,face = "italic", hjust = 1),
#     axis.text.x   = element_text(color = "black", size = 16, hjust=1,angle = 45),
#     axis.text.y   = element_text(color = "black", size = 16, angle = 0),
#     axis.title.x  = element_blank(),#element_text(color = "black", size = 16, angle = 0),
#     axis.title.y  =  element_text(color = "black", size = 20, angle = 90),
#     legend.title  = element_text(color = "black", size  = 16),
#     legend.text   = element_text(color = "black", size   = 16),
#     strip.text = element_text(size = 14),
#     axis.line.y = element_line(color = "black", linetype = "solid"),
#     axis.line.x = element_line (color = "black",linetype = "solid"),
#     panel.border = element_rect(linetype = "solid", linewidth = 1.2,fill = NA) )
# 
# common_EO
#   ggsave(paste(dir,"/4_EOCRC_common_genes.pdf",sep = ""), width = 30, height = 7, dpi = 300)


