
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(ggnewscale)
library(ggtreeExtra)
library(MetBrewer)


df_data <- read_csv("df_countries_alone_tree_EO.csv")

df_data <- df_data %>%  arrange(TMB_Status)
df_data <- df_data %>%
  mutate(Gene = case_when(
    TMB_Status == "hyper" ~ paste0(Gene,'_H'),
    TMB_Status == "nonhyper" ~ paste0(Gene, "_N"),
    TRUE ~ Gene
  ))
df = df_data

genes <- df$Gene
newick_string <- paste("(", paste(genes, collapse = ","), ");", sep = "")
tree <- read.tree(text = newick_string)

expdata <- df %>% select(-TMB_Status) %>% 
  pivot_longer(-c(Gene)) #%>% select(-TMB_Status)
group <- df %>% select(Gene,TMB_Status) %>% mutate(group="group")


pdf(file=paste(dir,'/12_EO_uniquegenes_tree.pdf',sep = ''),width = 10,height = 10)
 p <- ggtree(tree, branch.length = "none", layout = "circular", linetype = 0, size = 2) +
  layout_fan(angle = 60) +  
   geom_fruit(data = expdata, geom = geom_tile,
             mapping = aes(y = Gene, x = name, fill = value),
             pwidth = 0.3, offset = 0.005,
             axis.params = list(axis = "x", text.angle = -90, text.size = 4, hjust = 0)) +
  scale_fill_gradient(low = "#ffeda0", high = "#f03b20",na.value = "white") +
  geom_tiplab(offset = 0.5, size = 2, color = "black") +
  new_scale_fill() +
  geom_fruit(data = group, geom = geom_tile,
             mapping = aes(y = Gene, x = TMB_Status, fill = TMB_Status),
             pwidth = 0.2, offset = 0.3) +
  scale_fill_manual(values = c("#85D4E3", "#7294D4", "#C18748")) + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

print(p)
dev.off()



df_data <- read_csv("df_countries_alone_tree_LO.csv")

df_data <- df_data %>%  arrange(TMB_Status)
df_data <- df_data %>%
  mutate(Gene = case_when(
    TMB_Status == "hyper" ~ paste0(Gene,'_H'),
    TMB_Status == "nonhyper" ~ paste0(Gene, "_N"),
    TRUE ~ Gene
  ))
df = df_data


genes <- df$Gene
newick_string <- paste("(", paste(genes, collapse = ","), ");", sep = "")
tree <- read.tree(text = newick_string)

expdata <- df %>% select(-TMB_Status) %>% 
  pivot_longer(-c(Gene)) #%>% select(-TMB_Status)
group <- df %>% select(Gene,TMB_Status) %>% mutate(group="group")


pdf(file=paste(dir,'/12_LO_uniquegenes_tree.pdf',sep = ''),width = 10,height = 10)
p <- ggtree(tree, branch.length = "none", layout = "circular", linetype = 0, size = 2) +
  layout_fan(angle = 60) +  
  geom_fruit(data = expdata, geom = geom_tile,
             mapping = aes(y = Gene, x = name, fill = value),
             pwidth = 0.3, offset = 0.005,
             axis.params = list(axis = "x", text.angle = -90, text.size = 4, hjust = 0)) +
  scale_fill_gradient(low = "#ffeda0", high = "#f03b20",na.value = "white") +
  geom_tiplab(offset = 0.5, size = 2, color = "black") +
  new_scale_fill() +
  geom_fruit(data = group, geom = geom_tile,
             mapping = aes(y = Gene, x = TMB_Status, fill = TMB_Status),
             pwidth = 0.2, offset = 0.3) +
  scale_fill_manual(values = c("#85D4E3", "#7294D4", "#C18748")) + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

print(p)
dev.off()

