### DTI (search by target proteins from DrugBank data)
###
# DTI = sif2igraph("stitch_interactions.sif", directed = FALSE)
# plot(DTI, vertex.size=6)  
library(dplyr)
library(ggplot2)
# BiocManager::install("OmnipathR")
library(OmnipathR)
library(igraph)
library(ggraph)
library(dbparser)
library(XML)
library(xml2)

getwd()
setwd('\Zika virus drug')  # \\ abs path



#### Old version "dbparser" 1.0.4 ####

# read_drugbank_xml_db("full database.xml")
# get_xml_db_rows("full database.xml")
# xml_data <- read_xml("full database.xml")
# drug_nodes <- xml_find_all(xml_data, ".//drug")
# # 
# # xml_data <- xmlParse("full database.xml")
# # 
# # root_node <- xmlRoot(xml_data)
# # 
# # # Extracting specific elements, assuming drug entries are under the tag <drug>
# # drugs <- xpathSApply(root_node, "//drug", xmlToList)
# # 
# # # Parsing each drug entry
# # drug_info <- lapply(drugs, function(drug) {
# #   name <- xmlValue(drug$name)
# #   description <- xmlValue(drug$description)
# #   return(c(Name = name, Description = description))
# # })




#### Latest version "dbparser" 2.0.1 ####
if(FALSE)
{
  dvobj <- parseDrugBank(db_path            = "full database.xml",
                         drug_options       = drug_node_options(),
                         parse_salts        = TRUE,
                         parse_products     = TRUE,
                         references_options = references_node_options(),
                         cett_options       = cett_nodes_options())
  # 
}
  
# # save(dvobj, file = "dvobj.Rdata")

load('dvobj.Rdata')

test <- dvobj$drugs

drugs <- dvobj$drugs$drug_interactions


# load All drugs data ---------------------------------------------------------------

# drugs <- dbparser::parse_drug() 
# drugs <- drugs() 
# save(drugs, file='drugs.Rdata')


# load('drugs.Rdata')
drugs <- dvobj$drugs$general_information
drugs_data <- drugs %>% as_tibble() %>% dplyr::select(primary_key, name, type, state, average_mass, monoisotopic_mass)
colnames(drugs_data)[2] <- "drug_name"

test <- dvobj$drugs


target <- dvobj$cett$targets$polypeptides$general_information ### "cett": carriers, enzymes,targets and transporters
target_gene_map <- target %>% as_tibble() %>% dplyr::select(id, gene_name, parent_id)
# colnames(target_gene_map)[3] = 'id'

drug_target_map <- dvobj$cett$targets$general_information
drug_target_map <- drug_target_map %>% as_tibble() %>% dplyr::select(id, parent_key)
drug_target_map <- drug_target_map[, c(2,1)]
colnames(drug_target_map)[2] = 'parent_id'




### join drug-target and filter NA
drug_target_joined <- inner_join(drug_target_map, target_gene_map, by = "parent_id")
drug_target_joined <- drug_target_joined %>% filter(gene_name != "" & !is.na(gene_name))

map_gene_prot <- drug_target_joined[, c(3,4)]
drug_target_joined <- drug_target_joined[, c(1,4)]




### DTI_attr edge
DTI_attr = dvobj$cett$targets$general_information %>% as_tibble() %>% dplyr::select(id, parent_key, position)
DTI_attr_map <- DTI_attr[, c(2,1,3)]
colnames(DTI_attr_map)[2] = 'parent_id'

DTI_attr_map <- inner_join(DTI_attr_map, target_gene_map, by = "parent_id")
DTI_attr_map <- DTI_attr_map %>% filter(gene_name != "" & !is.na(gene_name))

# map_gene_prot <- DTI_attr_map[, c(4,5 , 3)]
DTI_attr <- DTI_attr_map[, c(1,5,3)]

# write.csv(DTI_attr, file = "DTI_attr.csv")



### updated DEG
DEG_update <- read.table('DEG_total_filtered.txt') ### some DEGs miss match in STRING
DEG_names


### filter DTI by updated DEG (936 DTIs)
drug_target_filtered <- filter(drug_target_joined, gene_name %in% DEG_update$V1) ### DTI: 667 nodes, 936 edges


### original not changed (in use)
# save(map_gene_prot, file = "map_gene_prot.RData")

unique(drug_target_filtered$gene_name)  ### 141 prots in DTI
unique(drug_target_filtered$parent_key)  ### 796 drugs in DTI







### DEGs (filtered by STRING PPI)
# DEG_names <- read.table('DEG_total.txt') ### some DEGs miss match in STRING
library(influential)
PPI_Network_1 = sif2igraph("string_interactions_short.tsv.sif", directed = FALSE)
DEG_names <- as_ids(V(PPI_Network_1))

unique(DEG_names)  ### 667 prots in PPI


DEG_update <- read.table('DEG_total_filtered.txt') ### some DEGs miss match in STRING

DEG_names

# V(PPI_Network_1)!= setdiff(DEG_names, DEG_update$V1)


### updated graph
PPI_Network_1 <- subgraph(PPI_Network_1, V(PPI_Network_1)[name %in% DEG_update$V1])





#### Drug-protein ####
# drug_targets_1 <- parse_drug_targets()
# drug_targets_1 <- drug_targets_1[,c(1,2,3,6)]
# ## load polypeptide data
# 
# drug_peptides_1 <- parse_drug_targets_polypeptides()
# drug_peptides_1 <- drug_peptides_1[,c(1,3,6,20)] # id, name, general_function, specific_function, gene_name, parent_id
# 
# # join the 3 datasets
# drug_targets_full <- inner_join(drug_targets_1, drug_peptides_1,
#                                 by=c("id"="parent_id", "name")) # %>% inner_join(drugs, by=c("parent_key"="primary_key"))
# 
# drug_target_inter <- drug_targets_full %>% dplyr::select(parent_key, gene_name)
# 
# 
# # save(drug_target_inter, file = "drug_target_inter.Rdata")  ### save All "drug_target_inter" from database
# 
# 
# 
# 
# # Filter drugs ---------------------------------------------------------------
# 
# ### protein of interest (POI)
# load("drug_target_inter.Rdata")

drug_target_inter <- drug_target_filtered

drug_target_inter   ### 936 DTIs,  796 drugs



## read all drug property
drug_calc_prop <- dvobj$drugs$calculated_properties

# save(drug_calc_prop, file = "drug_calc_prop.Rdata")
# load("drug_calc_prop.Rdata")


# load("gene_name.Rdata")
# # write.table(gene_name, file = "gene_name.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

gene_name <- DEG_names  ### 667 DEGs
# save(gene_name, file = "gene_name.Rdata")


### All drugs
drug_target_inter_data <- drug_target_inter   ### All DTI
drug_SMILES <- (drug_calc_prop[which(drug_calc_prop[,1] == "SMILES"  ),]) [, c(4,2)]  # All Drugs SMILES


# protein = gene_name
protein = DEG_update$V1  ### updated 659 DEGs
drug_nam <- filter(drug_target_inter_data, gene_name %in% protein)  ### filter by given POI
drug_nam <- unique(drug_nam$parent_key)    # 703 drugs in DTI

# drug_nam <- unique(drug_target_inter_data$parent_key)    ### 796 drugs





### drugs filtered by SMILES
drug_common <- intersect(drug_nam, drug_SMILES$parent_key) ### Drugs (chemical, SMILE, DTI)

length(drug_common) # updated: 703 drugs with SMILES



### 825 DTI interactions (filtered by "drug_common" and "DEG_update")
drug_target_inter_data <- filter(drug_target_inter, parent_key %in% drug_common)
drug_target_inter_data <- filter(drug_target_inter_data, gene_name %in% protein)


length(unique(drug_target_inter_data$parent_key))   ### 703 drugs
length(unique(drug_target_inter_data$gene_name))   ### 135 POI

# save(drug_common, file = 'drug_common.Rdata')
# save(drug_target_inter_data, file = 'drug_target_inter_data.Rdata')





### Total DDI drugs
# save(drug_interaction, file = "drug_interaction.Rdata")




load("drug_name.Rdata") ### from DDI.R
drug_name   ## 1009 drugs

# setdiff(unique(drug_target_inter_data$parent_key), DTI_data_frame$V1)





# DTI ---------------------------------------------------------------

# load("drug_target_inter_data.Rdata")
load("drug_target_inter_data.Rdata")
drug_target_inter_data      ### 825 DTI interactions

DTI_graph = graph_from_data_frame(drug_target_inter_data, directed = FALSE)   ### 838 nodes (703 drugs + 135 prot), 825 DTI edges
# save(drug_target_inter_data, file = 'DTI_graph.Rdata')


# load("drug_name.Rdata")
# drug_name

# load("drug_common.Rdata")
drug_nam <- unique(drug_target_inter_data$parent_key)    ### 703 drugs

prot_DTI <- unique(drug_target_inter_data$gene_name)    ### 135 drugs


































library(RandomWalkRestartMH)
MultiplexObject <- create.multiplex(list(DTI = DTI_graph))  ### DTI: 278 nodes, 295 edges

AdjMatrix_DTI <- compute.adjacency.matrix(MultiplexObject)
DTI <- as.matrix(AdjMatrix_DTI)
DTI <- DTI[order(rownames(DTI)), order(colnames(DTI))]  ###  matrix 278 x 278
DTI <- DTI[startsWith(rownames(DTI),"DB"), !startsWith(colnames(DTI),"DB")] ### 214 x 64, row = drugs, col = proteins
rownames(DTI)<-gsub("_1","",rownames(DTI))
colnames(DTI)<-gsub("_1","",colnames(DTI))
# save(DTI, file = "DTI_known.Rdata")  ### Known DTI 214 x 64



add_zero_col <- matrix(0, nrow = nrow(DTI), ncol = abs( dim(DTI)[2] - length(protein)  ) )
add_zero_col <- as.data.frame(add_zero_col)

colnames(DTI)<-gsub("_1","",colnames(DTI))
rownames(DTI)<-gsub("_1","",rownames(DTI))

colnames(add_zero_col) <- setdiff(protein, colnames(DTI))

DTI <- cbind(DTI, add_zero_col)  ### DTI: 214 drug x 495 protein
rownames(DTI)<-gsub("_1","",rownames(DTI))

# save(DTI, file = "DTI_all.Rdata")



DTI_drug_order <- rownames(DTI)   # 214 drug names
# save(DTI_drug_order,file = "DTI_drug_order.Rdata")











# load("drug_SMILES_296.Rdata")
# load("drug_SMILES.Rdata")
# drug_SMILES_all
# 
# # DTI_drugs_30 <- filter(drug_SMILES, parent_key %in% rownames(DTI))
# 
# drug_SMILES$parent_key   ### 296 drugs
# drug_common   ### 214 drugs
# drug_nam <- unique(drug_target_inter_data$parent_key) ###  == drug_common
# 






###
###
### Save DTI_tibble (Python)
load("DTI_known.Rdata")  ### Known DTI 214 x 64
DTI
dim(DTI)

DTI_all.Rdata

load("DTI_all.Rdata")  ### Known DTI 214 x 495
DTI
dim(DTI)


load('drug_SMILES_468.Rdata')
drug_SMILES
drug_SMILES_468 = drug_SMILES
dim(drug_SMILES_468)




# load("jaccard_sim_DTI_drug.Rdata")
# dim(jaccard_sim_DTI_drugs)
# jaccard_sim
# 
# intersect(drug_nam, drug_SMILES$parent_key) # 214 drugs
# intersect(drug_nam, drug_common) # 214 drugs


### 296 - 30 
# add_zero_row <- matrix(0, nrow = length(drug_SMILES$parent_key) - dim(DTI)[1] , ncol = ncol(DTI))  ### 258 - 36 = 222 rows

add_zero_row <- matrix(0, nrow = length(drug_SMILES$parent_key) - dim(DTI)[1] , ncol = ncol(DTI) )  ### 468 - 214 = 254 rows
add_zero_row <- as.data.frame(add_zero_row)   ### 254 x 495



rownames(add_zero_row) <- setdiff(drug_SMILES$parent_key, rownames(DTI))  ### 254 x 495
colnames(add_zero_row) <- colnames(DTI)

# DTI <- rbind(DTI, add_zero_row)   # (228 + 36) x 9
# DTI <- DTI[order(rownames(DTI)), order(colnames(DTI))]  ###  matrix 37 x 37
DTI <- rbind(DTI, add_zero_row)   # 468 x 495
DTI <- DTI[order(rownames(DTI)), order(colnames(DTI))]  ###  468 x 495
# DTI <- subset(DTI, rownames(DTI) %in% rownames(jaccard_sim)  )
# DTI <- DTI[order(rownames(DTI)) , order(colnames(DTI))]

dim(DTI)

plot(DTI, vertex.size=1)

DTI_tibble <- as_tibble(DTI)
rownames(DTI_tibble) <- rownames(DTI)

# save(DTI, file = "DTI.Rdata")
# save(DTI_tibble, file = "DTI_tibble.Rdata")

