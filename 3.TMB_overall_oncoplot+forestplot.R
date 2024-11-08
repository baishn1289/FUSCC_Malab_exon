library(tidyverse)
library(patchwork)
library(scales)
library(foreach)



format_pvalue <- function(p) {
  ifelse(p < 0.001, format(p, scientific = TRUE, digits = 3), 
         format(p, digits = 3, nsmall = 3))
}

custom_round <- function(x, digits) {
  posneg <- sign(x)
  z <- abs(x)*10^digits
  z <- z + 0.5
  z <- trunc(z)
  z <- z/10^digits
  z * posneg
}

# 定义百分数格式化函数
format_frequency <- function(freq) {
  # 将数值乘以100，并使用自定义四舍五入函数保留一位小数
  rounded_value <- custom_round(freq * 100, 1)
  
  paste0(format(rounded_value, nsmall = 1), "%")
  
}

format_adjPval <- function(p) {
  format(p, digits = 3, scientific = TRUE)
}



{
  
  #uniquegenes()函数主要是筛选出两个队列共有的基因（并对该突变基因的最低样本数量做要求）
  uniquegenes <- function(m1, m2, m1Name = NULL, m2Name = NULL,
                          minMut = 5){
    
    # numMutatedSamples = maf[!Variant_Type %in% 'CNV', .(MutatedSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
    # numAlteredSamples = maf[, .(AlteredSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
    # numAlteredSamples = merge(numMutatedSamples, numAlteredSamples, by = 'Hugo_Symbol', all = TRUE)
    
    m1.gs <- getGeneSummary(x = m1)
    m2.gs <- getGeneSummary(x = m2)
    
    if (is.null(m1Name)) {
      m1Name = "M1"
    }
    if (is.null(m2Name)) {
      m2Name = "M2"
    }
    
    #有可能出现虽然总样本中，老年和青年都有5个以上突变数
    #但是细分到具体的panel中，可能panel中只有1个年龄类别
    #解决办法：对于1个基因，如果某个panel只有一个age_class，弃去该基因在该panel中的样本
    #但注意，国家、种族、性别可能也会出现这样的问题
    
    
    #筛选出突变样本数大于等于minMut的基因，并生成两个数据集的基因列表
    #相当于筛掉人群中突变数太少的基因
    m1.genes = as.character(m1.gs[AlteredSamples >= minMut, 
                                  Hugo_Symbol])
    m2.genes = as.character(m2.gs[AlteredSamples >= minMut, 
                                  Hugo_Symbol])
    uniquegenes =  intersect(m1.genes, m2.genes)
    
    return(uniquegenes)
    
  }
  
  
  
  cl$Age_class <-  factor(cl$Age_class,levels = c('LOCRC','EOCRC'))
  rt@clinical.data$Age_class <-  factor(rt@clinical.data$Age_class
                                        ,levels = c('LOCRC','EOCRC'))
  

# 
#   Gene_test = 'TP53'
#   m1=clin.EOCRC
#   m2= clin.LOCRC
#   Clinical.data=rt@clinical.data
#   
  #将各个基因单独提取出来进行讨论和逻辑回归
  
  logstic_mafcompare <- function(m1, m2, m1Name = NULL, m2Name = NULL, Clinical.data,Gene_test) {
    
    
    # 提取有变异的样本
    # https://github.com/PoisonAlien/maftools/blob/master/R/summarizeMaf.R
    # if(!is.data.frame(maf)){
    #   #try to coerce into data.frame
    #   maf = data.table::as.data.table(maf)
    # }
    # numMutatedSamples = maf[!Variant_Type %in% 'CNV', .(MutatedSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
    # numAlteredSamples = maf[, .(AlteredSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
    # numAlteredSamples = merge(numMutatedSamples, numAlteredSamples, by = 'Hugo_Symbol', all = TRUE)
    
    #从maf@data中提取出来的Tumor_Sample_Barcode要去重，因为一个样本可能在一个基因上有多个突变外显子，这导致了一个TSB可能会被多次记录，不能代表样本数
    Alter_m1_Gene = m1@data[Hugo_Symbol == Gene_test,] %>% 
      filter(!duplicated(Tumor_Sample_Barcode))
    
    Alter_m2_Gene = m2@data[Hugo_Symbol == Gene_test,] %>% 
      filter(!duplicated(Tumor_Sample_Barcode))
    
    #提取包含该基因对应的panel，划定逻辑回归计算的患者范围
    #每个panel至少包含一名青年和一名老年患者，否则逻辑回归模型报错/拟合不佳
    #Age_class中可能因未匹配上出现NA值，会影响Age_class_count ==2的计算，须事先去除
    Gene_panels <- panel_hugo_symbol[Hugo_Symbol == Gene_test] %>%
      left_join(Clinical.data[, c('Tumor_Sample_Barcode', 'Age_class', 'Sex', 'Country', 'Race')], 
                by = "Tumor_Sample_Barcode") %>%
      filter(!duplicated(Tumor_Sample_Barcode),
             !is.na(Age_class)) %>%
      group_by(Panel) %>%
      reframe(
        Age_class_count = n_distinct(Age_class)
      ) %>%
      filter(Age_class_count ==2 ) %>% 
      pull(Panel)
  
  
    
    
    if (length(Gene_panels) == 0) {
      print(paste(Gene_test, '没有一个中心有两种性别', sep = ' '))
      
      singlesex_panel_record_log <<- rbind(singlesex_panel_record_log, 
                                           data.frame(Gene_test = Gene_test,
                                                      stringsAsFactors = FALSE))
      
      return()
    }
    
   
    
    #基因突变频率计算
    #突变频率按照样本突变频率计算：突变样本数量/突变基因所在panel的样本总数
   
 
    #提出该基因变异的样本数量
    mut_samples_m1 <- length(unique(Alter_m1_Gene$Tumor_Sample_Barcode))
    
    
    #提出该基因所在的panels的所有样本数量
    panels_samples_m1 <- length(unique(m1@clinical.data[m1@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
    
    
    m1_freq <- mut_samples_m1 / panels_samples_m1
    
    # ------------
    mut_samples_m2 <- length(unique(Alter_m2_Gene$Tumor_Sample_Barcode))
    
    
    panels_samples_m2 <- length(unique(m2@clinical.data[m2@clinical.data$Panel %in% Gene_panels, ]$Tumor_Sample_Barcode))
    
    
    m2_freq <- mut_samples_m2 / panels_samples_m2

    # 更新临床信息以备逻辑回归：
    #1.筛选所在panels；
    #2.以肿瘤样本数为准做逻辑回归，不筛除多样本的患者；
    #3.标记某个样本是否含有该基因变异
    
    Gene_clinical_info  <-  Clinical.data[, Gene_status := ifelse(Tumor_Sample_Barcode %in% c(Alter_m1_Gene$Tumor_Sample_Barcode,
                                                                                              Alter_m2_Gene$Tumor_Sample_Barcode), 1, 0)]%>% 
      left_join(TMB_mutload[,.(Tumor_Sample_Barcode,TMB)],by = 'Tumor_Sample_Barcode') %>% 
      mutate(Age_class = factor(Age_class,levels = c('LOCRC','EOCRC'))) %>% 
      filter(Panel %in% Gene_panels) 
    
    
    # 逻辑回归模型
    
    if (length(unique(Gene_clinical_info$Race)) > 1) {
      formula_base <- paste(formula_base, "+ Race")
    }
    if (length(unique(Gene_clinical_info$Panel)) > 1) {
      formula_base <- paste(formula_base, "+ Panel")
    }
    if (length(unique(Gene_clinical_info$Country)) > 1) {
      formula_base <- paste(formula_base, "+ Country")
    }
    if (length(unique(Gene_clinical_info$Sex)) > 1) {
      formula_base <- paste(formula_base, "+ Sex")
    }
    
    
    # 执行逻辑回归
    fit <- glm(as.formula(formula_base), family = binomial(link = "logit"), data = Gene_clinical_info)
    
    
    # summary(fit)
  
    # 提取OR值、P值和置信区间
    # summary_fit = summary(fit)
    # OR_vals = exp(coef(fit))
    # pvals = coef(summary_fit)[, "Pr(>|z|)"]
    
    # 提取 Age_classEOCRC 的 OR（Odds Ratio）
    OR_age_class <- exp(coef(fit)["Age_classEOCRC"])
    
    # 提取 Age_classEOCRC 的 p 值
    pval_age_class <- summary(fit)$coefficients["Age_classEOCRC", "Pr(>|z|)"]
    
    # confint(fit)用profile likelihood估计CI，大样本计算太费时间；改用标准误差计算CI
    # 提取 Age_classEOCRC 的置信区间
    # 如果计算时间充足，大规模运算可以使用以下代码
    # conf_interval_age_class <- exp(confint(fit, parm = "Age_classEOCRC"))
    # ci_low <- conf_interval_age_class[1]
    # ci_up  <- conf_interval_age_class[2]
    coef_age_class <- coef(fit)["Age_classEOCRC"]
    se_age_class <- summary(fit)$coefficients["Age_classEOCRC", "Std. Error"]
    # 基于标准误差的近似计算置信区间
    ci_low <- exp(coef_age_class - 1.96 * se_age_class)
    ci_up  <- exp(coef_age_class + 1.96 * se_age_class)
    
    # 构建结果数据框
    df_Gene_test <- data.table::data.table(
      Hugo_Symbol = Gene_test,
      EOCRC_total_case= panels_samples_m1,
      LOCRC_total_case= panels_samples_m2,
      EOCRC = mut_samples_m1, 
      LOCRC = mut_samples_m2, 
      EOCRC_freq = m1_freq,
      LOCRC_freq = m2_freq,
      pval = pval_age_class, 
      or = OR_age_class,
      ci.up = ci_up,
      ci.low = ci_low,
      formula = formula_base
    )
    
    
    result_df <<- rbind(result_df, df_Gene_test)  
    
    
  }
  

  #提取出所有目标基因
  #minMut暂时调整到10
  overall_genes <- uniquegenes(clin.EOCRC,clin.LOCRC,minMut = 10)
  
  #统计各个基因所拥有的panel(x)、国家(x)、种族(x)、性别(x)等临床信息
  Hugo_clinical_info <- panel_hugo_symbol[Hugo_Symbol %in% overall_genes] %>%
    merge(rt@clinical.data[, c('Tumor_Sample_Barcode', 'Age_class','Country',"Sex","Race")], 
          by = "Tumor_Sample_Barcode") %>%
    group_by(Hugo_Symbol) %>%
    summarise(
      Panel_count = n_distinct(Panel),
      Country_count = n_distinct(Country),
      Race_count = n_distinct(Race),
      Sex_count = n_distinct(Sex)
    ) %>% ungroup()

  names(Hugo_clinical_info)
  # 确保每个基因的各个临床特征因子至少有2个水平
  # 但是这边进入逻辑回归计算过程的基因，有可能只是包含了部分的panel、国家、race、sex
  Hugo_clinical_info_xxxx <- Hugo_clinical_info %>% 
    pull(Hugo_Symbol)
  
 
  
  singlesex_panel_record_log <- data.frame(Gene_test = character(),
                                           stringsAsFactors = FALSE)
  
  
  problematic_genes <- c() 
  
  
  result_df <- data.table::data.table(
    Hugo_Symbol = NA,
    EOCRC_total_case= NA,
    LOCRC_total_case= NA,
    EOCRC = NA, 
    LOCRC = NA, 
    EOCRC_freq = NA,
    LOCRC_freq = NA,
    # EOCRC_freq_normalized = NA,
    # LOCRC_freq_normalized = NA,
    # EOCRC_panel_TMB= NA,
    # LOCRC_panel_TMB= NA,
    pval = NA, 
    or = NA,
    ci.up =  NA,
    ci.low =  NA,
    formula =NA
  )
  
  formula_base <- "Gene_status ~ Age_class+ TMB"
 
  # 
  # library(future.apply)
  # plan(multisession, workers = 16)
  # options(future.globals.maxSize = 2.5 * 1024^3)
  
  
  # results_list <- future_lapply(Hugo_clinical_info_xxxx, function(i) {
  #   tryCatch({
  #     # 调用自定义函数并返回结果
  #     result <- logstic_mafcompare(m1 = clin.EOCRC, 
  #                                  m2 = clin.LOCRC, 
  #                                  Clinical.data = linchuangshuju, 
  #                                  Gene_test = i)
  #     return(result)  # 返回计算结果
  #   }, warning = function(w) {
  #     problematic_genes <<- c(problematic_genes, i)  # 记录产生警告的基因
  #     return(NULL)  # 遇到警告时返回 NULL
  #   })
  # })
  
  
  for (i in Hugo_clinical_info_xxxx) {
    result <- tryCatch({
      logstic_mafcompare(m1 = clin.EOCRC, m2 = clin.LOCRC, Clinical.data = rt@clinical.data, Gene_test = i)
      NULL  # 如果没有错误，返回NULL
    }, warning = function(w) {
      problematic_genes <<- c(problematic_genes, i)  # 记录产生警告的基因
      return(NULL)  # 返回NULL以继续循环
    })
  }
  
  
  result_df_overall <- result_df
  
  result_logistic <- result_df_overall[-1,] %>% 
    filter(!Hugo_Symbol %in% problematic_genes) %>% 
    filter(formula == 'Gene_status ~ Age_class+ TMB + Race + Panel + Country + Sex') %>% 
    mutate(adjPval = p.adjust(p = pval, method = "fdr")) 
  
  
  rownames(result_logistic) <- seq_len(nrow(result_logistic))
  
  save(result_df,problematic_genes,singlesex_panel_record_log,
       file='overall_logistic_result.Rdata')
  
  write.table(result_logistic, file='./overall_EOCRC_vs_LOCRC_freq.tsv', 
              quote=FALSE, row.names=FALSE, sep="\t")
  
  
} 


  {#overall_forest
    
    df_overall <- read_tsv('overall_EOCRC_vs_LOCRC_freq.tsv')
    
    df_overall <- df_overall %>% 
      filter(!is.infinite(or) & !is.infinite(ci.low) & !is.infinite(ci.up))
    
    #set p_value and adjustPval if needed
    #p_value_threshold <- 0.05
    adjustPval_threshold <- 0.05
    
    filtered_df_overall <- df_overall %>% 
      filter(
        adjPval <= adjustPval_threshold) %>%
      # filter(
      #   LOCRC_freq< EOCRC_freq
      # ) %>%
      mutate(
        Freq_average = (LOCRC_freq+EOCRC_freq )/2,
        Freq_diff = EOCRC_freq - LOCRC_freq  >0 ,
        Freq_normalized_diff = EOCRC_freq_normalized - LOCRC_freq_normalized >0
      ) %>%
      arrange(desc(Freq_average)) 
    
    write.csv(filtered_df_overall,
              file = paste(dir,'/3.overall_differential_gene.csv',sep=''),
              row.names = F)
    
    table(filtered_df_overall$Freq_diff != filtered_df_overall$Freq_normalized_diff)
    
    filtered_df_overall <-filtered_df_overall %>% head(20)
    
    
    # add new labels
    plot <- filtered_df_overall %>% 
      mutate(EOCRC_freq = format_frequency(EOCRC_freq),
             LOCRC_freq = format_frequency(LOCRC_freq),
             sample_size = paste0(LOCRC_freq, " v.s. ", EOCRC_freq),
             OR_CI = paste0(format(round(or, 2), nsmall = 2), " (", format(round(ci.low, 2), nsmall = 2), " to ", format(round(ci.up, 2), nsmall = 2), ")"),
             adjPval = format_adjPval(adjPval)) %>% 
      mutate(across(everything(), as.character)) %>% 
      bind_rows(data.frame(Hugo_Symbol = "Gene", sample_size = "Frequency (LOCRC v.s. EOCRC)", adjPval = "adjPval", OR_CI = "OR (95% CI)")) %>% 
      mutate(Hugo_Symbol = fct_rev(fct_relevel(Hugo_Symbol, "Gene")))
    
    # middle part
    p_mid <- filtered_df_overall %>% 
      ggplot(aes(y = fct_rev(Hugo_Symbol))) +
      theme_classic() +
      geom_point(aes(x = or, size = or), shape = 23, fill = "black", show.legend = F) +
      geom_errorbarh(aes(xmin = ci.low, xmax = ci.up), height = 0.3) +
      labs(x = "Odds Ratio") +
      # scale_x_continuous(trans = 'log10', breaks = c(1.0,1.5,2.0,2.5),
      #                    labels = scales::number_format(accuracy = 0.1)) + 
      # coord_cartesian(ylim = c(1, 21), xlim = c(1.0,2.5)) + 
      scale_x_continuous( trans = 'log10',breaks = c(0.2,0.5,1,2,4,8,16,32),
                          labels = scales::number_format(accuracy = 0.1)) +
      coord_cartesian(ylim = c(1, 21), xlim = c(0.2,32)) +
      geom_vline(xintercept = 1, linetype = "dashed") +  
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(color = "black"))
    
    p_left <- plot %>%
      ggplot(aes(y = Hugo_Symbol)) +
      geom_text(aes(x = 0, label = Hugo_Symbol), hjust = 0, fontface = "bold", size = 4) +
      geom_text(aes(x = 1, label = sample_size), hjust = 0, size = 4, vjust = 0,
                fontface = ifelse(plot$sample_size == "Frequency (LOCRC v.s. EOCRC)", "bold", "plain")) +
      theme_void() + coord_cartesian(xlim = c(0, 3))
    
    # right part
    p_right <- plot %>% ggplot() +
      geom_text(aes(x = 0, y = Hugo_Symbol, label = OR_CI), hjust = 0,
                fontface = ifelse(plot$OR_CI == "OR (95% CI)", "bold", "plain"), size = 4) +
      geom_text(aes(x = 1, y = Hugo_Symbol, label = adjPval), hjust = 0,
                fontface = ifelse(plot$adjPval == "adjPval", "bold", "plain"), size = 4) +
      theme_void() + coord_cartesian(xlim = c(0, 2))
    
    
    layout <- c(
      area(t = 0, l = 0, b = 30, r = 9), 
      area(t = 0, l = 10, b = 30, r = 17), 
      area(t = 0, l = 18, b = 30, r = 25)) 
    
    
  }


pdf(paste(dir,'/3_overall_forestPlot_freq_EOLOadded20.pdf',sep = ""), 
    width = 11, height = 6)
print(p_left + p_mid + p_right + plot_layout(design = layout))

dev.off()
