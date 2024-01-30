getwd()
setwd('h:\\Users\\yrt05\\Desktop\\Zika virus drug')  # \\ abs path
# setwd("C:/Users/yrt05/Desktop/Systems Biology_Stroke-20210522T190533Z-001/Systems Biology_Stroke") # \\ abs path

#download data from GEO database, get file ".csv"
#function source: https://github.com/jmzeng1314/humanid/blob/master/R/downGSE.R
if(T){
  downGSE <- function(studyID = "GSE235210", destdir = ".") {
    library(GEOquery)
    eSet <- getGEO(studyID, destdir = destdir, getGPL = T)
    
    exprSet = exprs(eSet[[1]])
    pdata = pData(eSet[[1]])
    

    return(eSet)
  }}


library(GEOquery)
eSet <- downGSE('GSE235210')


library(dplyr)

### Group names (constrast design matrix)
data_group = pData(eSet[[1]]) ### 141 samples, 

# colnames(data)

library(stringr)
group_list = as.character(data_group$characteristics_ch1.3)

# group_list_replace <- gsub("zika infection: ", "", str_trim(group_list))
group_list_replace <- gsub("-", "_", str_trim(group_list))
group_list_replace <- gsub(" ", "_", str_trim(group_list_replace))


group_list_replace <- gsub("treatment:_", "", str_trim(group_list_replace))

group_list_replace <- gsub("_MOI_1", "", str_trim(group_list_replace))

group_list_replace[1:3] <- paste(group_list_replace[1:3], "12H")

group_list_replace[4:6] <- paste(group_list_replace[4:6], "24H")

group_list_replace[7:9] <- paste(group_list_replace[7:9], "48H")

group_list_replace[10:12] <- paste(group_list_replace[10:12], "72H")

group_list_replace <- gsub(" ", "_", str_trim(group_list_replace))


group_list_replace





### Row count

# Specify the path to the TAR file containing .txt.gz files
tar_file_path <- "h:\\Users\\yrt05\\Desktop\\Zika virus drug\\GSE235210_RAW.tar"
# Specify the directory where you want to extract the files
extract_dir <- "h:\\Users\\yrt05\\Desktop\\Zika virus drug\\extracted_files\\"

# Create the extraction directory if it doesn't exist
dir.create(extract_dir, showWarnings = FALSE)

# Extract the TAR file
untar(tar_file_path, exdir = extract_dir)

# List all .txt.gz files in the extraction directory
txt_gz_files <- list.files(extract_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)

my_tibble <- tibble()
# combined_columns <- tibble()

n <- 25071  # Number of rows
c <- 1  # Number of columns

my_matrix <- matrix(0, nrow = n, ncol = c)
my_tibble <- as_tibble(my_matrix)

# Loop through the .txt.gz files, unzip, and combine them
for (file in txt_gz_files) {
  unzipped_file <- paste(file)
  # Read the unzipped file as a tibble (assuming it's tab-delimited)
  unzipped_data <- read.table(unzipped_file, header = FALSE, sep = "\t")
  
  col_add <- unzipped_data[2]

  result_string <- sub("^.*ZIKV_MOI", "ZIKV_MOI", unzipped_file)

  result_string <- gsub("ZIKV_MOI-1_", "", str_trim(result_string))
  result_string <- gsub("ZIKV_MOI_1_", "", str_trim(result_string))
  
  result_string <- str_extract(result_string, ".*_RNA")
  
  result_string <- gsub("_RNA", "", str_trim(result_string))
  
  colnames(col_add) <- result_string
  
  my_tibble <- cbind(my_tibble, col_add)
  
}
my_tibble <- my_tibble[2:ncol(my_tibble)]
rownames(my_tibble) <- unzipped_data$V1

library(readr)
library(stringr)
library(dplyr)
data <- my_tibble






#### Normalization (Voom) ####


# BiocManager::install("limma")
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list_replace))


colnames(design) = levels(factor(group_list_replace))

pdata = data
rownames(design)=colnames(pdata)

# save(design,file='design.Rdata')

constrast_group <- c('Mock - ZIKV_infected_12H',  'Mock - ZIKV_infected_24H', 'Mock - ZIKV_infected_48H', "Mock - ZIKV_infected_72H")

#### Constrast Matrix ####
contrast.matrix <- makeContrasts( 'Mock - ZIKV_infected_12H',  'Mock - ZIKV_infected_24H', 'Mock - ZIKV_infected_48H', 'Mock - ZIKV_infected_72H', levels = design)
head(pdata)
rownames(pdata) = rownames(data)
v <- voom(pdata,design,normalize="none")  # normalization

head(v[[1]])

#### Find DEG ####
fit <- lmFit(v,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
tempOutput = topTable(fit3, coef=1, n=Inf, adjust="fdr")
DEG_voom = na.omit(tempOutput)
head(DEG_voom)


logFC_threshold = 0.58  # 1.5 fold
# logFC_threshold = 0.3785  # 1.3 fold
# logFC_threshold = 1
adj_p_val = 0.05
p_val = 0.05
find_DEG_up = DEG_voom[which(DEG_voom[1]> logFC_threshold & DEG_voom[4]<p_val),]    # logFC = 0.58 is value 1.5
find_DEG_down = DEG_voom[which(DEG_voom[1]< (-logFC_threshold) & DEG_voom[4]<p_val),]
# find_DEG_up = DEG_voom[which(DEG_voom[1]> logFC_threshold & DEG_voom[5]<adj_p_val),]    # logFC = 0.58 is value 1.5
# find_DEG_down = DEG_voom[which(DEG_voom[1]< (-logFC_threshold) & DEG_voom[5]<adj_p_val),]


print(dim(find_DEG_up))  ### 983
print(dim(find_DEG_down)) ### 715




DEG_up = rownames(find_DEG_up)      
DEG_down = rownames(find_DEG_down)  
DEG_total = c(DEG_up,DEG_down)     
length(DEG_total) 

#### updated in "Similarity_matrix_protein.R"
# write(DEG_gene_name,"DEG_total_map.txt",sep = '')  # save TXT data

# DEG_total <- DEG_update$V1 ### updated 659 DEGs


# # ### save DEG
# save(find_DEG_up,file = 'find_DEG_up.Rdata')
# save(find_DEG_down,file = 'find_DEG_down.Rdata')






#### Filter gene probe ####

data_sample_probe <- v[[1]]
data_sample_probe <- as.data.frame(data_sample_probe)
library(dplyr)
DEG_voom %>% distinct()
DEG_voom <- DEG_voom[!duplicated(rownames(DEG_voom)), ] 


data_sample_probe <- data_sample_probe[!duplicated(rownames(data_sample_probe)), ] 



library(rentrez)

exprSet <- data_sample_probe
# dat=exprSet
# exprSet=dat
### colorful boxplot
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(exprSet, col = cols,main="expression value",las=2)
# boxplot(exprSet[,5])


### Map gene name

# Example using biomaRt for Ensembl IDs
# if (!requireNamespace("biomaRt", quietly = TRUE))
#   install.packages("biomaRt")
library(biomaRt)
# remotes::install_version("dbplyr", version = "2.3.4")
library(dbplyr)

library(remotes)



# Connect to the Ensembl database
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Define your probe ID
probe_id <- DEG_total
# probe_id <- 'ENSMNEG00000033800'
# Retrieve gene information

# gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
#                    filters = "ensembl_gene_id",
#                    values = probe_id,
#                    mart = ensembl)
gene_info <- DEG_total

print(gene_info)



DEG_gene_name <- DEG_total
dim(DEG_total)


length(DEG_total)




#### heatmap (no need) ####
library(pheatmap)
choose_gene = head(DEG_total,30)
choose_matrix = exprSet[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
pheatmap(choose_matrix, )
exprSet['Ptrh1',]





#### volcano plot (no need) #### 
#### use DEG_voom not DEG
nrDEG = DEG_voom
DEG=nrDEG

# DEG = DEG_P
DEG$logFC
# logFC_threshold = 0.58
# p_val = 0.05
# logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
DEG$change = as.factor(ifelse(DEG$P.Value < p_val & abs(DEG$logFC) > logFC_threshold,
                              ifelse(DEG$logFC > logFC_threshold ,'Up-regulated','Down-regulated'),'NOT')
)

DEG_voom$Change = as.factor(ifelse(DEG_voom$P.Value < p_val & abs(DEG_voom$logFC) > logFC_threshold,
                                   ifelse(DEG_voom$logFC > logFC_threshold ,'Up-regulated','Down-regulated'),'Not Significant')
)
# round(logFC_cutoff,3)
this_tile <- paste0('Cutoff:' , 'Fold Change > 1.5,', 'P Value < ', p_val,
                    '\nNumber of up-regulated gene: ',nrow(DEG[DEG$change =='Up-regulated',]) ,
                    '\nNumber of down-regulated gene: ',nrow(DEG[DEG$change =='Down-regulated',])
)

font_size = 15

this_tile
# install.packages("ggplot2")
library("ggplot2")
head(DEG)
DEG_voom
g = ggplot(data=DEG_voom, aes(x=logFC, y=-log10(P.Value), color=Change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=font_size)))+
  xlab("log2 fold change") + ylab("-log10 P-value") +
  theme(plot.title = element_text(size=font_size,hjust = 0))+
  scale_colour_manual(values = c('blue','black','red')) + 
  geom_hline(yintercept= -log10(p_val) -0.05, linetype="dashed", color = "black", size=2) + 
  geom_vline(xintercept = logFC_threshold - 0.05, linetype="dashed", color = "black", size=2) +
  geom_vline(xintercept = -logFC_threshold + 0.05, linetype="dashed", color = "black", size=2) + 
  theme(
    axis.title.x = element_text(size = font_size),
    axis.text.x = element_text(size = font_size),
    axis.title.y = element_text(size = font_size),
    axis.text.y = element_text(size = font_size),
    legend.title = element_text(size=font_size),
    legend.text = element_text(size=font_size)
  )
print(g)
















#### Pathway mapping ####

library("org.Hs.eg.db")

# BiocManager::install("KEGGREST")

library("KEGGREST")

hsa_path_eg  <- keggLink("pathway", "hsa") %>% 
  tibble(pathway = ., eg = sub("hsa:", "", names(.)))

hsa_kegg_anno <- hsa_path_eg %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),
    ensembl = mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID")
  )


#### Select a pathway ####
# hsa_kegg_anno_hsa05168 <- hsa_kegg_anno[hsa_kegg_anno$pathway == 'path:hsa05168',]
# 
# length(intersect(hsa_kegg_anno$symbol, DEG_total)) ### 777 DEGs are mapped to database
# 
# DEG_hsa05168 <- intersect(hsa_kegg_anno_hsa05168$symbol, DEG_total)
# # write(DEG_hsa05168, "DEG_hsa05168.txt")





########################## online pathway analyses

### STRING
# BiocManager::install("STRINGdb")
library(STRINGdb)

# data(diff_exp_example1)
head(DEG_voom)
rownames(DEG_voom)
nrDEG <- DEG_voom%>%filter(rownames(DEG_voom)%in%DEG_total)
gene_List = c(nrDEG$P.Value,nrDEG$logFC)
names(gene_List) = rownames(nrDEG)






### read filtered 659 protein names
# DEG_update <- read.table('DEG_total_filtered.txt') ### some DEGs miss match in STRING
# gene_name <- DEG_update$V1
# gene_name   # 659 proteins (updated)
# DEG_total <- gene_name

####################### Enrichment KEGG
library(clusterProfiler)
# data(geneList, package="DOSE") # Disease Ontology Semantic and Enrichment analysis
gene_getname.df <- bitr(DEG_total, fromType = "SYMBOL",
                        toType = c("ENTREZID"),
                        OrgDb = "org.Hs.eg.db")
# test <- mapIds(org.Hs.eg.db, keys = DEG_total, column = "SYMBOL", keytype = "ALIAS")
# 
# test <- mapIds(org.Hs.eg.db, 
#                    keys = DEG_total, 
#                    column = "ENTREZID",  # or "ENSEMBL" etc., depending on what ID type you need
#                    keytype = "SYMBOL",  # this is for gene symbols
#                    multiVals = "first")  # How to handle multiple matches
# test <- na.omit(test)

head(gene_getname.df)


library(org.Hs.eg.db)

kk <- enrichKEGG(gene = gene_getname.df$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.05)
# kk <- enrichKEGG(gene = test,
#                  organism = "hsa",
#                  pvalueCutoff = 0.05)

font_size = 15
# dotplot(kk, showCategory=10, x = 'Count', font.size = font_size, orderBy="x")



### Pathway plot

head(gene_getname.df)

font_size = 11
library(scales) 
library('cowplot')
library(ggplot2)
head(kk@result)
kk@result[["Description"]][1:3]

kegg <- kk
p<- dotplot(kegg, showCategory=5, x = 'Count', font.size = font_size, orderBy="x")

p1 <- p + theme(
  plot.title = element_text(color="black", size=font_size, face="bold", hjust=0.5))



gene_id <- bitr(rownames(nrDEG), fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = "org.Hs.eg.db")

fold_vector <- nrDEG$logFC
# fold_vector <- nrDEG[,1]

names(fold_vector) <- gene_id$ENTREZID


### Pathway
p2 <- cnetplot(kk, showCategory=5, foldChange=fold_vector, node_label = "category", layout = "kk",
               colorEdge = TRUE, color_category ='steelblue', node_label_size = NULL,
               cex_category = 1, label_format = 5, cex_label_category = 1)

p2 <- p2 + theme(legend.position="right",
                 legend.text = element_text(size = font_size),
                 # legend.background = element_rect(fill="lightblue", size = 1),
                 legend.title = element_text(size = font_size, face = "bold"),
                 legend.justification = c("right", "top"),
                 # legend.margin = margin(6, 6, 6, 6)
)

library(cowplot)

# p2 pathway
# p1 GO
plot_grid(p2, p1, labels = c("a","b"), label_size = 15,  nrow = 2, align = c("h"))







#### Significant DEGs ####
significant_gene_id <- kk@result[["geneID"]][1]

significant_gene_id <- gsub("/"," ", str_trim(significant_gene_id))

significant_gene_id <- str_split(significant_gene_id, " ")[[1]]

significant_gene <- bitr(significant_gene_id, fromType = "ENTREZID",
                         toType = c("SYMBOL"),
                         OrgDb = "org.Hs.eg.db")
significant_gene_KEGG <- significant_gene$SYMBOL


intersect(significant_gene_KEGG, selected_genes_name)

intersect(intersect(significant_gene_KEGG, selected_genes_name), DEG_WGCNA_common)


up_DEG <- intersect(significant_gene_KEGG, rownames(find_DEG_up))
down_DEG <- intersect(significant_gene_KEGG, rownames(find_DEG_down))
up_DEG <- bitr(up_DEG, fromType = "SYMBOL",
               toType = c("ENTREZID"),
               OrgDb = "org.Hs.eg.db")
down_DEG <- bitr(down_DEG, fromType = "SYMBOL",
                 toType = c("ENTREZID"),
                 OrgDb = "org.Hs.eg.db")
rownames(nrDEG)
P_up <- dat[rownames(dat) %in% up_DEG$SYMBOL,4]
P_down <- dat[rownames(dat) %in% down_DEG$SYMBOL,4]



### GO enrichment
# 
# Genes<-c(significant_gene_id)
# # g <- goana(Genes)
# # BiocManager::install("org.mmu.eg.db")
# # library(org.mmu.eg.db)
# g <- goana(list(Up=up_DEG$ENTREZID, Down=down_DEG$ENTREZID, P.Up = P_up, P.Down = P_down), 
#            species="Mm", covariate = TRUE, plot = TRUE)
# GO_result <- topGO(g)
# GO_result$Term
# 
# GO_result$

### GO enrichment 2
# library(clusterProfiler)
# # de <- c(up_DEG$ENTREZID, down_DEG$ENTREZID)
# de <- bitr(rownames(nrDEG), fromType = "SYMBOL",
#            toType = c("ENTREZID"),
#            OrgDb = "mouse430a2.db")
# de <- de$ENTREZID
# go <- enrichGO(de, OrgDb = "org.Mm.eg.db", ont="all")
# # ego <- enrichGO(de, OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
# # go <- enrichDAVID(de, idType="ENTREZ_GENE_ID",  listType="Gene", annotation="KEGG_PATHWAY")
# library(ggplot2)




library(ggplot2)
library(GOplot)
data(EC)
test <- EC$david
# setwd("C:/Users/yrt05/Desktop/Systems biology project/GSE137482_RAW")
# setwd('c:\\Users\\yrt05\\Desktop\\Covid_coference_presentation')
readfile <- read.csv("gProfiler_hsapiens.csv")
readfile <- subset(readfile, select = c("source","term_name", "intersections", "adjusted_p_value", "term_id") )
colnames(readfile)[1:4] <-colnames(test)[c(1,3,4,5)]
readfile <- readfile[, c(1,5,2,3,4)]
colnames(readfile)[2] <-colnames(test)[2]
test3 <- as.factor(readfile$Genes)
length(test3)
test3 <- str_split(test3, ",")
# unlist(test3)
library("org.Hs.eg.db")

nrDEG_data <- cbind(rownames(nrDEG), nrDEG)

colnames(nrDEG_data) <- c("ID", colnames(nrDEG))
readfile$Category <- gsub("GO:","", str_trim(readfile$Category))
rownames(nrDEG_data) <- c(1:nrow(nrDEG_data))

# length(test3)
for (i in 1:length(test3)){
  map <- mapIds(org.Hs.eg.db, keys = unlist(str_split(unlist(test3[i]),",")) , keytype = "SYMBOL", column = "SYMBOL", multiVals="first")
  # map <- mapIds(org.Mm.eg.db, keys = map , keytype = "SYMBOL", column = "PROBEID", multiVals="list")
  result <- data.frame(paste(toupper(map), collapse=', '))
  readfile$Genes[i] <- result
  # print(as.factor(paste(toupper(map), collapse=', ')))
}

readfile <- filter(readfile, Category %in% c("CC", "BP", "MF"))

circ <- circle_dat(readfile, nrDEG_data)


# GOBubble(circ, labels = 10)
# Reduce redundant terms with a gene overlap >= 0.75...
# reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# # ...and plot it
# GOBubble(reduced_circ, labels = 2.8)

# GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 10)  

# GOBar(subset(circ, category == 'BP'))







#### GO plot ####

### Enrichment GO
GO_enrich <- enrichGO(gene_getname.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")


gene_getname.df

go <- GO_enrich

go2 <- dotplot(GO_enrich, showCategory=5, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

### Separate analyses
go_result <- summary(go)
go_result_All <- go_result$ID[order(go_result$Count, decreasing=T)]
# head(go_result_All)
go_result_All[1:5]
IDs <- c(go_result_All[1:5])
# install.packages('GOplot')
library(GOplot)
go1 <- GOCircle(circ, nsub = IDs,  rad1 = 2, rad2 = 3,  zsc.col = c('yellow', 'black', 'cyan'))

# go_result_BP <- filter(go_result, ONTOLOGY %in% c("BP"))
# go_result_BP <- go_result_BP$ID[order(go_result_BP$Count, decreasing=T)]
# go_result_MF <- filter(go_result, ONTOLOGY %in% c("MF"))
# go_result_MF <- go_result_MF$ID[order(go_result_MF$Count, decreasing=T)]
# go_result_CC <- filter(go_result, ONTOLOGY %in% c("CC"))
# go_result_CC <- go_result_CC$ID[order(go_result_CC$Count, decreasing=T)]


go2 <- go2 + theme(legend.position="right",
                   legend.text = element_text(size = font_size),
                   # legend.background = element_rect(fill="lightblue", size = 1),
                   legend.title = element_text(size = font_size, face = "bold"),
                   legend.justification = c("right", "top"),
                   # legend.margin = margin(6, 6, 6, 6)
)
# plot(go2)
library(cowplot)
# ,   rel_heights = c(0.3, 1.2)
# par(mfrow=c(1,2))
# dotplot(go, showCategory=5, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
# GOCircle(circ, nsub = IDs, label.size = 6, rad1 = 3, rad2 = 4,  zsc.col = c('yellow', 'black', 'cyan'))

plot_grid(go2, go1, labels = c("c","d"), label_size = 22, ncol = 2, align = c("h"), greedy = TRUE, 
          rel_widths = c(1.2, 1),
          axis = "b", scale = c(1,1), cex= 2)









