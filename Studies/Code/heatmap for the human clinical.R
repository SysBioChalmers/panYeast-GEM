
# heastmap-superheat
library(tidyverse)
library(readxl)
library(ggplot2)
library(plyr)
require(reshape2)
require(scales)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(readr)
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## 注意缺省值为 'hjust=0' 和'rot=270'
}
#然后将default draw_colnames覆盖
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))


human_clinical <- read.table("data/human_clinical.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
#save to change the format
write.table(human_clinical, "data/human_clinical.txt", row.names = FALSE, sep = "\t")
human_clinical0 <- read_excel("data/human_clinical.xlsx")
rowname0 <- human_clinical0$X
human_clinical1 <- human_clinical0[,-1]
rownames(human_clinical1) <- rowname0

pheatmap(human_clinical1,
         scale = 'none',
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = TRUE,
         cluster_col = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize_col = 7,
         fontsize_row=5,
         border_color = "LightGrey", cellwidth = 7, cellheight = 4,
         color = colorRampPalette(c("white", "SandyBrown", "firebrick3"))(100))




# for human group
human <- read_excel("data/human_group.xlsx")
rowname0 <- human$strain
human1 <- human[,-1]
rownames(human1) <- rowname0

pheatmap(human1,
         scale = 'none',
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = TRUE,
         cluster_col = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize_col = 8,
         fontsize_row=8,
         border_color = "LightGrey", cellwidth = 16, cellheight = 7,
         color = colorRampPalette(c("white", "SandyBrown", "firebrick3"))(100))
