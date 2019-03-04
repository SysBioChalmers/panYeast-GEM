# library
library(readxl)
library(hongR)
library(tidyverse)
library(tibble)
library(readr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
#-------------------------------------------------------------------------
#reaction and gene number analysis for all the models
#-------------------------------------------------------------------------

# Load strain specific model
specificModelResultFile <- read_excel("data/specificModelResultFile.xlsx")
#reactiion number analysis
strain_specific_rxn <- select(specificModelResultFile, strain, rxns)
summary(strain_specific_rxn$rxns)
sd(strain_specific_rxn$rxns)
#strain_specific_rxn$rxns <- as.factor(strain_specific_rxn$rxns)
ggplot(strain_specific_rxn, aes(rxns)) + geom_bar(fill = "blue") +
  xlab('Reaction number') + ylab('Strain number') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size=20,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  ggsave(out <- paste('result/','rxn number for all strain','.eps', sep = ""), width=6, height=6, dpi=300)


#gene number analysis
strain_specific_gene <- select(specificModelResultFile, strain, genes)
#strain_specific_gene$genes <- as.factor(strain_specific_rxn$genes)
ggplot(strain_specific_gene, aes(genes)) + geom_bar(fill = "blue") +
  xlab('Gene number') + ylab('Strain number') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size=20,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  ggsave(out <- paste('result/','gene number for all strain','.eps', sep = ""), width=6, height=6, dpi=300)


#Metabolite number analysis
strain_specific_met <- select(specificModelResultFile, strain, mets)
#strain_specific_met$mets <- as.factor(strain_specific_met$mets)
ggplot(strain_specific_met, aes(mets)) + geom_bar(fill = "blue") +
  xlab('Metabolite number') + ylab('Strain number') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size=20,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  ggsave(out <- paste('result/','met number for all strain','.eps', sep = ""), width=6, height=6, dpi=300)




#rxn number analysis based on the strain classification
strain_classification <- read_excel("data/strain_classification.xls", sheet = "Table S1")
unique(strain_classification$Ecological_origins)
unique(strain_classification$Clades)
#reload the data
strain_specific_rxn <- read.table("data/strain specific model-rxn number", header = FALSE, sep = ";", stringsAsFactors = FALSE)
strain_specific_rxn$Ecological_origins <- getSingleReactionFormula(strain_classification$Ecological_origins,strain_classification$Standardized_name,strain_specific_rxn$V1)
strain_specific_rxn$Clades <- getSingleReactionFormula(strain_classification$Clades,strain_classification$Standardized_name,strain_specific_rxn$V1)
strain_specific_rxn$Ecological_origins <- as.factor(strain_specific_rxn$Ecological_origins)
strain_specific_rxn$V2 <- as.numeric(strain_specific_rxn$V2)
boxplot(V2~Ecological_origins,data=strain_specific_rxn, main="", col=terrain.colors(4),
        xlab="Mutation classification", ylab="Relative growth rate")


p <-  ggplot(strain_specific_rxn, aes( x = Ecological_origins, y = V2, fill = Ecological_origins))
p + geom_boxplot(show.legend = FALSE) +
  xlab('') + ylab('Reaction number') +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=18)) +
  theme(axis.text=element_text(size=15,face="bold", family="Arial"),
        axis.title=element_text(size=20,face="bold", family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.position="none")

ggsave(out <- paste('result/','Reaction number','.eps', sep = ""), width=8, height=4, dpi=300)




#-------------------------------------------------------------------------
# substrate usage analysis for all the models
#-------------------------------------------------------------------------

substrate <- read_tsv('data/Biolog_Substrate.tsv')
#substrate_use <- read.table('data/results for substrate usage.txt', stringsAsFactors = FALSE)
strain_list <- strain_classification$Standardized_name
substrate_list <- paste(substrate$Substrate, substrate$Substrate_type, sep='@')
substrate_usage <- read_csv("data/substrate_usage.csv", 
                            col_names = FALSE)
substrate_usage0 <- as.data.frame(t(substrate_usage), stringsAsFactors = FALSE)
rownames(substrate_usage0) <- strain_list
colnames(substrate_usage0) <- substrate_list
s <- substrate_usage0$`α-D-Glucose@C`
s1 <- t(substrate_usage0[1, ])
choosed_index <- which(s1 %in% "NG"== FALSE)
choosed_strain <- substrate_list[choosed_index]

substrate_usage1 <- substrate_usage0[,choosed_strain]

#change the string into number
substrate_usage1[, 1:221] <- sapply(substrate_usage1[, 1:221], as.numeric)
zero_index <- colSums(substrate_usage1)
zero_index <- zero_index[zero_index!=0]
choosed_strain1 <- names(zero_index)
substrate_usage2 <- substrate_usage0[,choosed_strain1]
#analysis
result0 <- data.frame(strain_name=strain_list, stringsAsFactors = FALSE)
used_substrate <- vector()
for (i in 1:nrow(substrate_usage2)){
  print(i)
  ss <- as.numeric(unlist(substrate_usage2[i, ]))
  if(length(which(ss != 0)) !=0){
  used_substrate[i] <- length(which(ss != 0))
  } else{
    used_substrate[i] <- 0
  }
}

result0$num_substrate_use <- used_substrate
result0$Ecological_origins <- getSingleReactionFormula(strain_classification$Ecological_origins,strain_classification$Standardized_name,result0$strain_name)
result0$Clades <- getSingleReactionFormula(strain_classification$Clades,strain_classification$Standardized_name,result0$strain_name)
result0$Ecological_origins <- as.factor(result0$Ecological_origins)
result0$Clades <- as.factor(result0$Clades)


# boxplot
p <-  ggplot(result0, aes( x = Ecological_origins, y = num_substrate_use, fill = Ecological_origins))
p + geom_boxplot(show.legend = FALSE) +
  xlab('') + ylab( 'Number of substrates used') +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  theme(axis.text=element_text(size=15,face="bold", family="Arial"),
        axis.title=element_text(size=20,face="bold", family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.position="none")

ggsave(out <- paste('result/','Number of substrates used','.eps', sep = ""), width=8, height=4, dpi=300)






write.table(result0, "result/substrate usage analysis for 1011 strains.txt", row.names = FALSE, sep = "\t")
s0 <- filter(result0, Ecological_origins=='Human')



#-------------------------------------------------------------------------
# decision tree based on the substrate usage
#-------------------------------------------------------------------------

# decision tree based on only the carbon source usage
# prepare the data
library("rpart")
library("rpart.plot")
library("ParamHelpers")
library("mlr")
substrate0 <- colnames(substrate_usage2)
index0 <- which(str_detect(substrate0, "@C"))
carbon_usage <- substrate_usage2[, substrate0[index0]]
carbon_usage[,1:59] <- sapply(carbon_usage[,1:59], as.numeric)
carbon_usage$Ecological_origins <- result0$Ecological_origins
str(carbon_usage)
indexes = sample(1011, 808)
carbon_usage_train = carbon_usage[indexes,]
carbon_usage_test = carbon_usage[-indexes,]
target = Ecological_origins ~.
tree = rpart(target, data = carbon_usage, minsplit=32, minbucket=5, cp=0.001)
rpart.plot(tree,extra= 106,legend.x=FALSE, legend.y=FALSE)

# tune the hyperparameters
colnames(carbon_usage) <- paste("s",1:60, sep = "") 
trainTask <- makeClassifTask(data = carbon_usage, target = "s60")
testTask <- makeClassifTask(data = carbon_usage, target = "s60")
#Search for hyperparameters
getParamSet("classif.rpart")

#make tree learner
makeatree <- makeLearner("classif.rpart", predict.type = "response")
#set 3 fold cross validation
set_cv <- makeResampleDesc("CV",iters = 3L)
#Search for hyperparameters
gs <- makeParamSet(
  makeIntegerParam("minsplit",lower = 10, upper = 50),
  makeIntegerParam("minbucket", lower = 5, upper = 50),
  makeNumericParam("cp", lower = 0.001, upper = 0.2)
)
#do a grid search
gscontrol <- makeTuneControlGrid()

#hypertune the parameters
stune <- tuneParams(learner = makeatree, resampling = set_cv, task = trainTask, par.set = gs, control = gscontrol, measures = acc)

target = s60 ~.
tree = rpart(target, data = carbon_usage, minsplit=32, minbucket=5, cp=0.001)
rpart.plot(tree,extra= 106,legend.x=FALSE, legend.y=FALSE)
ggsave(out <- paste('result/','Decision tree','.eps', sep = ""), width=8, height=4, dpi=300)




#-------------------------------------------------------------------------
# simple analysis based on yield of metabolite
#-------------------------------------------------------------------------
metabolite_type_for_yield <- read_excel("data/metabolite_type_for_yield.xlsx")
yield_of_metabolites <- read_csv("data/yield of metabolites.csv", col_names = FALSE)
colnames1 <- metabolite_type_for_yield$met_type
colnames1 <- append('strain_name', colnames1, after = 1)
colnames1 <- str_replace_all(colnames1,"-","_")
colnames(yield_of_metabolites) <- colnames1
yield_of_metabolites$strain_name <- str_replace_all(yield_of_metabolites$strain_name, "\\.mat","")

yield_of_metabolites$Ecological_origins <- getSingleReactionFormula(strain_classification$Ecological_origins,strain_classification$Standardized_name,yield_of_metabolites$strain_name)
yield_of_metabolites$Clades <- getSingleReactionFormula(strain_classification$Clades,strain_classification$Standardized_name,yield_of_metabolites$strain_name)
yield_of_metabolites$Ecological_origins <- as.factor(yield_of_metabolites$Ecological_origins)
yield_of_metabolites$Clades <- as.factor(yield_of_metabolites$Clades)

# vioplot
p <- ggplot(yield_of_metabolites, aes(factor(Ecological_origins), growth_minimalmedia))
p + geom_violin(aes(fill = factor(Ecological_origins)), show.legend = FALSE)+
  xlab('') + ylab( 'Yield of biomass') +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text=element_text(size=15,face="bold", family="Arial"),
        axis.title=element_text(size=20,face="bold", family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.position="none")

ggsave(out <- paste('result/','Yield of biomass','.eps', sep = ""), width=8, height=4, dpi=300)
yield_of_metabolites0 <- filter(yield_of_metabolites, Ecological_origins=='Human, clinical')

yield_analysis <- data.frame(unclass(summary(yield_of_metabolites)), stringsAsFactors = FALSE)
yield_analysis <- yield_analysis[,-c(1,28,29,30,31)]
yield_of_metabolites01 <- yield_of_metabolites[,-c(1,28,29,30,31)]
product_analysis <- data.frame(substrate=colnames(yield_of_metabolites01), stringsAsFactors = FALSE)
for(i in 1:length(colnames(yield_of_metabolites01))){
  s0 <- colnames(yield_of_metabolites01)[i] 
  print(s0)
  product_analysis$min[i] <- min(yield_of_metabolites01[,s0])
  product_analysis$mean[i] <- mean(unlist(yield_of_metabolites01[,s0]), na.rm = TRUE)
  product_analysis$max[i] <- max(yield_of_metabolites01[,s0])
}



#雷达图
devtools::install_github("ricardo-bion/ggradar", 
                         dependencies=TRUE)
library(ggradar)
library(scales)

suppressPackageStartupMessages(library(dplyr))
library(scales)
library(tibble)

#data change
product_analysis0 <- as.data.frame(t(product_analysis), stringsAsFactors = FALSE)
colnames(product_analysis0) <- product_analysis0[1,] %>% str_replace_all(.,'L_asparagine','asn') %>% str_replace_all(.,'isobutanol','butanol')%>%str_replace_all(.,'L_glutamine','gln') %>%
  str_replace_all(.,'L_','') %>%
  str_replace_all(.,'\\(S\\)_','') %>% str_extract(.,"^.{3}") %>% str_replace_all(.,'but','butanol')

product_analysis0 <- product_analysis0[-1,]
for (i in colnames(product_analysis0)){
  product_analysis0[,i] <- as.numeric(product_analysis0[,i])
}

product_analysis0 <- product_analysis0 %>%
  rownames_to_column( var = "group" )
product_analysis0 <- product_analysis0[-2,]



ggradar(product_analysis0,
        base.size=15,
        values.radar = c("0", "1.5", "3"),
        grid.min=0,  #10,
        grid.mid=1.5,  #50,
        grid.max=3, #100
        gridline.min.colour="grey",
        gridline.mid.colour="grey",
        gridline.max.colour="black",
        gridline.max.linetype="solid",
        grid.label.size=6,
        group.point.size=3,
        group.line.width=1,
        background.circle.colour='white',
        group.colours=c('OrangeRed','MediumBlue'),
        axis.line.colour="Gainsboro",
        axis.label.size=7
)







#------------------------------------------------------------------------------------------------------
# pca analysis based on yield data, the gene matrix and reactiion matrix in 1011 strains specific models
#------------------------------------------------------------------------------------------------------
library(ggfortify)
df1 <- yield_of_metabolites[,colnames1[-c(1,28,29)]]
autoplot(prcomp(df1))
autoplot(prcomp(df1), data = yield_of_metabolites, colour = 'Clades')
# pca 3D map
library(pca3d)
df2 <- df1[ , apply(df1, 2, var) != 0]
pca <- prcomp(df2, scale.=TRUE)
pca3d(pca, group=yield_of_metabolites$Clades, legend="bottomleft")
snapshotPCA3d(file="ellipses.png") # save the 3D PCA map


# pca analysis based on the reaction matrix in 1011 strains specific models
rxn_matrix <- read.delim2('data/rxn_matrix_1011_gem.txt', header = FALSE, stringsAsFactors = FALSE)
rxn_matrix0 <- as.data.frame(t(rxn_matrix), stringsAsFactors = FALSE)
rxn_matrix01 <- rxn_matrix0
rxn_matrix0$strain_name <- strain_list
rxn_matrix0$Ecological_origins <- getSingleReactionFormula(strain_classification$Ecological_origins,strain_classification$Standardized_name,rxn_matrix0$strain_name)
rxn_matrix0$Clades <- getSingleReactionFormula(strain_classification$Clades,strain_classification$Standardized_name,rxn_matrix0$strain_name)
rxn_matrix0$Ecological_origins <- as.factor(rxn_matrix0$Ecological_origins)
rxn_matrix0$Clades <- as.factor(rxn_matrix0$Clades)
autoplot(prcomp(rxn_matrix01), data = rxn_matrix0, colour = 'Ecological_origins') +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))  +
  theme(legend.text=element_text(size=15))

ggsave(out <- paste('result/','PCA-Rxn','.png', sep = ""), width=12, height=8, dpi=300)



# pca analysis based on the gene matrix  matrix in 1011 strains specific models
gene_matrix <- read.delim2('data/gene_matrix_1011_gem.txt', header = FALSE, stringsAsFactors = FALSE)
gene_matrix0 <- as.data.frame(t(gene_matrix), stringsAsFactors = FALSE)
gene_matrix01 <- gene_matrix0
gene_matrix0$strain_name <- strain_list
gene_matrix0$Ecological_origins <- getSingleReactionFormula(strain_classification$Ecological_origins,strain_classification$Standardized_name,gene_matrix0$strain_name)
gene_matrix0$Clades <- getSingleReactionFormula(strain_classification$Clades,strain_classification$Standardized_name,gene_matrix0$strain_name)
gene_matrix0$Ecological_origins <- as.factor(gene_matrix0$Ecological_origins)
gene_matrix0$Clades <- as.factor(gene_matrix0$Clades)
autoplot(prcomp(gene_matrix01), data = gene_matrix0,  colour = 'Ecological_origins') +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
      axis.title=element_text(size=24,face="bold", family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.text=element_text(size=15))
ggsave(out <- paste('result/','PCA-gene','.png', sep = ""), width=12, height=8, dpi=300)

#obtain the pca information
pca_results <- prcomp(gene_matrix01)
pca_results$rotation[,1:4]
#if we choose to remove the outliers, some data point in the right part






#-----------------------------------------------------------------------
# heatmap for the industrial strain based on the yield of 26 metabolites
#-----------------------------------------------------------------------
yield_of_metabolites0 <- filter(yield_of_metabolites, Ecological_origins=="Industrial")
yield_of_metabolites1 <- yield_of_metabolites0[,8:27]
rownames(yield_of_metabolites1) <- yield_of_metabolites0$strain_name
s1 <- c("ala","arg","asn","asp","cys","gln","glu","gly","his","ile","leu","lys","met","phe","pro","ser","thr","trp","tyr","val")
colnames(yield_of_metabolites1) <- s1
#improve the quality of graph
library(grid)
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## 注意缺省值为 'hjust=0' 和'rot=270'
}
#然后将default draw_colnames覆盖
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

pheatmap(yield_of_metabolites1,
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = FALSE,
         cluster_col = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 15,
         color = colorRampPalette(c("white", "SandyBrown", "firebrick3"))(100))

ggsave(out <- paste('result/','Yield of metabolite for industrial','.eps', sep = ""), width=8, height=6, dpi=300)












#-----------------------------------------------------------
# easy cluster
#-----------------------------------------------------------
d <- dist(gene_matrix01, method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1, cex = 0.6, hang = -1)




#-------------------------------------------------------------------------
# example1
# hclust analysis based on dendextend package
#-------------------------------------------------------------------------
iris <- datasets::iris
iris2 <- iris[,-5]
species_labels <- iris[,5]
library(colorspace) # get nice colors
species_col <- rev(rainbow_hcl(3))[as.numeric(species_labels)]

pairs(iris2, col = species_col,
      lower.panel = NULL,
      cex.labels=2, pch=19, cex = 1.2)

# Add a legend
par(xpd = TRUE)
legend(x = 0.05, y = 0.4, cex = 2,
       legend = as.character(levels(species_labels)),
       fill = unique(species_col))
par(xpd = NA)

# http://blog.safaribooksonline.com/2014/03/31/mastering-parallel-coordinate-charts-r/
par(las = 1, mar = c(4.5, 3, 3, 2) + 0.1, cex = .8)
MASS::parcoord(iris2, col = species_col, var.label = TRUE, lwd = 2)

# Add Title
title("Parallel coordinates plot of the Iris data")
# Add a legend
par(xpd = TRUE)
legend(x = 1.75, y = -.25, cex = 1,
       legend = as.character(levels(species_labels)),
       fill = unique(species_col), horiz = TRUE)
par(xpd = NA)
d_iris <- dist(iris2) # method="man" # is a bit better
hc_iris <- hclust(d_iris, method = "complete")
iris_species <- rev(levels(iris[,5]))
library(dendextend)
dend <- as.dendrogram(hc_iris)
# order it the closest we can to the order of the observations:
dend <- rotate(dend, 1:150)

#Color the branches based on the clusters:
dend <- color_branches(dend, k=3) #, groupLabels=iris_species)

# Manually match the labels, as much as possible, to the real classification of the flowers:
labels_colors(dend) <-
  rainbow_hcl(3)[sort_levels_values(
    as.numeric(iris[,5])[order.dendrogram(dend)]
  )]

# We shall add the flower type to the labels:
labels(dend) <- paste(as.character(iris[,5])[order.dendrogram(dend)],
                      "(",labels(dend),")", 
                      sep = "")
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
# dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.5)
# And plot:
par(mar = c(3,3,3,7))
plot(dend, 
     main = "Clustered Iris data set
     (the labels give the true flower species)", 
     horiz =  TRUE,  nodePar = list(cex = .007))
legend("topleft", legend = iris_species, fill = rainbow_hcl(3))
# Requires that the circlize package will be installed
par(mar = rep(0,4))
circlize_dendrogram(dend)



#-------------------------------------------------------------------------
# example2
# decision tree analysis based on dendextend package
#-------------------------------------------------------------------------
install.packages("rpart")
install.packages("rpart.plot")
library("rpart")
library("rpart.plot")
data("iris")
indexes = sample(150, 110)
iris_train = iris[indexes,]
iris_test = iris[-indexes,]
target = Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width
tree = rpart(target, data = iris_train, method = "class")
rpart.plot(tree)



#------------------------------------------------------------------
# example3
# pca base on binary method
#------------------------------------------------------------------
library(logisticPCA)
s0 <- data("house_votes84")
