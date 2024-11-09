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

# genelist <- genelist[genelist %in% colnames(gene_mutation_counts)]
# gene_mutation_counts <- gene_mutation_counts[, c('Tumor_Sample_Barcode', genelist)]

rownames(gene_mutation_counts) <- gene_mutation_counts$Tumor_Sample_Barcode
gene_mutation_counts <- gene_mutation_counts[,-1]

gene_mutation_counts <- gene_mutation_counts[, colnames(gene_mutation_counts) %in% genelist]

final <- maf.mutload %>%
  left_join(gene_mutation_counts, by = "Tumor_Sample_Barcode")

colnames(final)[1:20]
final <- final[!is.na(final$TMB),]

gene_panel_map <- panel %>% 
  group_by(Hugo_Symbol) %>% 
  summarize(Genes = list(SEQ_ASSAY_ID), .groups = "drop") %>%
  deframe()  # 将数据框转换为列表，便于查找

results <- data.frame(Gene = character(), 
                      Overall_pval = numeric(), EOCRC_rate = numeric(), LOCRC_rate = numeric(),
                      Hypermutated_pval = numeric(), Hypermutated_EOCRC_rate = numeric(), Hypermutated_LOCRC_rate = numeric(),
                      Nonhypermutated_pval = numeric(), Nonhypermutated_EOCRC_rate = numeric(), Nonhypermutated_LOCRC_rate = numeric(),
                      stringsAsFactors = FALSE)

dim(gene_mutation_counts)

for(i in colnames(gene_mutation_counts)) { #  gene
  print(i)
  final_data <- final %>% filter(Panel %in% gene_panel_map[[i]])
  
  final_data[,i] <- final_data[,i]/final_data$TMB
  
  x1 <- final_data[final_data$Age_class == 'EOCRC', i]
  x2 <- final_data[final_data$Age_class == 'LOCRC', i]
  
  if(length(x1) >= 10 & length(x2) >= 10) {
    EOCRC_rate <- mean(x1)
    LOCRC_rate <- mean(x2)
    
    # Overall
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
