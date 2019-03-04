# library
library(readxl)
library(hongR)
library(tidyverse)
library(tibble)
# Load gene
pangene_and_coregene <- read_excel("data/pangene and coregene.xlsx")

number_ticks <- function(n) {
  function(limits) pretty(limits, n)
}
ggplot(pangene_and_coregene, aes(num)) + 
  geom_line(aes(y = coregene, colour = "Core gene")) + 
  geom_line(aes(y = pangene, colour = "Pan gene")) +
  geom_line(aes(y = accessory_gene, colour = "Accessory gene")) +
  
  scale_x_continuous(breaks = number_ticks(10)) +
  scale_y_continuous(limits = c(1000, 8000), breaks = number_ticks(10)) +
  labs(x="Sampled strain number",y="Number") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=15)) +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.position="none")

ggsave(out <- paste('result/','pangene and coregene','.png', sep = ""), width=8, height=6, dpi=300)







# change the format of geneMatrix
# load gene list for 1011 yeast strains
geneMatrix <- read.table('data/genesMatrix_PresenceAbsence.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
g1 <- geneMatrix[,2:7797]
geneName <- colnames(g1)
strainName <- geneMatrix[,1]
geneMatrix0 <- as.data.frame(t(g1), stringsAsFactors = FALSE)

#name
rownames(geneMatrix0) <- geneName
colnames(geneMatrix0) <- strainName
geneMatrix0 <- add_column(geneMatrix0, geneName, .before = "BFC")
write.table(geneMatrix0, "result/geneMatrix0 of 1011 yeast strains.txt", row.names = FALSE, sep = "\t")





