library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
library(maftools)

maf.mutload <- read.csv('TMB_mutload.csv', row.names = 1)

# 仅保留需要的列
mutation_data <- as.data.table(rt@data[, 2:30]) 
dim(mutation_data)

# 计算每个样本的基因突变计数
gene_mutation_counts <- mutation_data %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  summarise(Mutation_Count = n(), .groups = "drop") %>%
  spread(key = Hugo_Symbol, value = Mutation_Count, fill = 0) %>% as.data.frame

rownames(gene_mutation_counts) <- gene_mutation_counts$Tumor_Sample_Barcode
gene_mutation_counts <- gene_mutation_counts[,-1]

gene_mutation_counts <- gene_mutation_counts[, colnames(gene_mutation_counts) %in% genelist]

panel_gene_map <- panel %>% 
  group_by(SEQ_ASSAY_ID) %>% 
  summarize(Genes = list(Hugo_Symbol), .groups = "drop") %>%
  deframe()

results <- data.frame(Gene = character(), 
                      Overall_pval = numeric(), EOCRC_rate = numeric(), LOCRC_rate = numeric(),
                      Hypermutated_pval = numeric(), Hypermutated_EOCRC_rate = numeric(), Hypermutated_LOCRC_rate = numeric(),
                      Nonhypermutated_pval = numeric(), Nonhypermutated_EOCRC_rate = numeric(), Nonhypermutated_LOCRC_rate = numeric(),
                      stringsAsFactors = FALSE)

dim(gene_mutation_counts)

for(i in colnames(gene_mutation_counts)) { #  gene
  print(i)
  for(j in rownames(gene_mutation_counts)) { # patient id
    refpanel <- maf.mutload$Panel[match(j, maf.mutload$Tumor_Sample_Barcode)]
    if(! i %in% panel_gene_map[[refpanel]])
      gene_mutation_counts[j, i] = NA
  }
  
  final_data <- data.frame(x1 = rownames(gene_mutation_counts), 
                   x2 = gene_mutation_counts[,i])
  
  colnames(final_data) <- c('Tumor_Sample_Barcode', i)
  
  final_data <- maf.mutload %>% left_join(final_data, by = "Tumor_Sample_Barcode")
  
  x1 <- c(na.omit(final_data[final_data$Age_class == 'EOCRC', i]))
  x2 <- c(na.omit(final_data[final_data$Age_class == 'LOCRC', i]))
  
  if(length(x1) >= 10 & length(x2) >= 10) {
    EOCRC_rate <- mean(x1)
    LOCRC_rate <- mean(x2)
    
    # Overall p 值
    overall_test <- wilcox.test(x1, x2)
    
    # Hypermutated
    x1 <- c(na.omit(final_data[final_data$Age_class == 'EOCRC' & final_data$TMB > 10, i]))
    x2 <- c(na.omit(final_data[final_data$Age_class == 'LOCRC' & final_data$TMB > 10, i]))
    Hypermutated_EOCRC_rate <- mean(x1)
    Hypermutated_LOCRC_rate <- mean(x2)
    hyper_test <- wilcox.test(x1, x2)
    
    # Nonhypermutated
    x1 <- c(na.omit(final_data[final_data$Age_class == 'EOCRC' & final_data$TMB <= 10, i]))
    x2 <- c(na.omit(final_data[final_data$Age_class == 'LOCRC' & final_data$TMB <= 10, i]))
    Nonhypermutated_EOCRC_rate <- mean(x1)
    Nonhypermutated_LOCRC_rate <- mean(x2)
    nonhyper_test <- wilcox.test(x1, x2)
    
    results <- rbind(results, data.frame(
      Gene = i,
      Overall_pval = overall_test$p.value, EOCRC_rate = EOCRC_rate, LOCRC_rate = LOCRC_rate,
      Hypermutated_pval = hyper_test$p.value, Hypermutated_EOCRC_rate = Hypermutated_EOCRC_rate, Hypermutated_LOCRC_rate = Hypermutated_LOCRC_rate,
      Nonhypermutated_pval = nonhyper_test$p.value, Nonhypermutated_EOCRC_rate = Nonhypermutated_EOCRC_rate, Nonhypermutated_LOCRC_rate = Nonhypermutated_LOCRC_rate
    ))
  }
}

write.csv(results, "gene_mutation_rate_pvalue_filtered.csv", row.names = FALSE)
