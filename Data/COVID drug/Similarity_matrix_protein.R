# BiocManager::install("impute")

library("protr")
setwd('\\')  # \\ abs path

# # load FASTA files


library("readxl")


## All 495 were mapped to protein based on UniProt database
protein_id <- read_excel("uniprot-protein_id.xlsx")
# protein_id <- protein_id$Entry

### protein of interest
# POI <- c("TLR2","ITGB3","CFP","PIK3CA","BAX","CCL5","IL6","IL1B","PTPN11") 
load("gene_name.Rdata")
POI = unique(gene_name)

gene_name_uniprot <- protein_id$`yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y`  ### 499 identifier, to 492 genes
# gene_name_uniprot <- protein_id$`Gene Names` 

# protein_id$`Entry name`
# name <- protein_id$`Gene names  (primary )`
# name <- gsub("_HUMAN*","",name)

diff <- setdiff(unique(gene_name_uniprot),POI)


`%!in%` <- Negate(`%in%`)
prot_diff <- filter(protein_id, `yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y` %!in% POI)
prot_diff_tibble <- add_row(prot_diff, prot_diff)
prot_diff_tibble$`yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y`[1] <- 'NXF2B'
prot_diff_tibble$`yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y`[3] <- 'NXF2'

prot_diff_tibble$`yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y`[2] <- 'ZNF224'
prot_diff_tibble$`yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y`[4] <- 'ZNF233'

protein_iden <- filter(protein_id, `yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y` %in% POI)

unique(protein_iden$`yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y`)


protein_iden <- add_row(protein_iden, prot_diff_tibble)

unique(protein_iden$`yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y`)

protein_iden <- filter(protein_iden, protein_iden$`yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y` %in% POI)
protein_iden <- protein_iden %>% distinct(`yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y`, .keep_all= TRUE)



prots <- getUniProt(protein_iden$Entry)   # 493 Proteins Sequence
names(prots) = protein_iden$`yourlist:M202203234ABAA9BC7178C81CEBC9459510EDDEA3472F29Y`
# save(prots, file = "prots.RData")

load("prots.RData")
prots




# Protein sequence feature------------------------------
# # install.packages('rentrez')
# library(rentrez)
# # shrooms <- rentrez::entrez_fetch(db = "protein",
# #                                  id = names(prots),
# #                                  rettype = "fasta",
# #                                  use_history=TRUE)
# 
# # install.packages('seqinr')
# library(seqinr)
# # shrooms_list <- read.fasta('uniprot-download_10_10.fasta')
# 
# 
# fasta_cleaner <- function(fasta_object, parse = TRUE){
# 
#   fasta_object <- sub("^(>)(.*?)(\\n)(.*)(\\n\\n)","\\4",fasta_object)
#   fasta_object <- gsub("\n", "", fasta_object)
# 
#   if(parse == TRUE){
#     fasta_object <- stringr::str_split(fasta_object,
#                                        pattern = "",
#                                        simplify = FALSE)
#   }
# 
#   return(fasta_object[[1]])
# }
# 
# i <- 1
# prots[i]
# fasta_cleaner(prots[[i]], parse = F)
# 
# 
# for(i in 1:length(shrooms_list)){
#   shrooms_list[[i]] <- fasta_cleaner(shrooms_list[[i]], parse = F)
# }
# 
# # appending shrooms_list to a vector
# shrooms_vector <- rep(NA, length(shrooms_list))
# 
# # accessing specific columns
# for(i in 1:length(shrooms_vector)){
#   shrooms_vector[i] <- shrooms_list[[i]]
# }
# 
# #  naming the vector
# names(shrooms_vector) <- names(shrooms_list)
# 
# 
# # magrittr::extract(DTI_dat_frame$gene_name)
# #
# # unique(magrittr::extract(DTI_dat_frame$gene_name))
# #
# # # prots_data = prots[1:10]
# # prots_data = prots
# # prot_PPI <- (prots_data %>%  magrittr::extract(DTI_dat_frame$gene_name)      )
# #
# # # prot_PPI <- (prots_data %>%  unique(magrittr::extract(DTI_dat_frame$gene_name))      )
# 
# 
# # BiocManager::install("msa")
# # library(msa)

# library(xlsx)
library(readxl)
load("prots.RData")
prots

prot_PPI <- prots
prot_PPI_names <- names(prot_PPI)
# test <- as.data.frame(prot_PPI)
# 
# mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
# mySequences <- readAAStringSet(mySequenceFile)
# head(mySequences)
# 
# Biostrings::AAStringSet(prot_PPI)
# 
# test <- msa(prot_PPI)



# Conjoint Triad  343  dim <<Predicting Protein-protein Interactions Based Only on Sequences Information>>
library(protr)

prot_property = matrix(, nrow = length(prot_PPI_names), ncol = length( names(extractCTriad(prot_PPI$PIK3CA))    ) )
i=1
for (seq in prot_PPI){
  encode_343d <- as.numeric(extractCTriad(seq))
  prot_property[i,] <- encode_343d
  i = i+1
}
prot_property <- as.data.frame(prot_property)
colnames(prot_property) <- names(extractCTriad(prot_PPI$PIK3CA))
row.names(prot_property) <- prot_PPI_names
prot_property <- prot_property[order(rownames(prot_property)) , order(colnames(prot_property))]



library(dplyr )

# additional properties
protein_data <- read_excel("uniprot-compressed_10_10.xlsx")
prot_data <- filter(protein_data, From %in% gene_name)
# prot_data <- prot_data %>% as_tibble() %>% dplyr::select(From, Mass, Features, Length, Sequence, `Gene Ontology (GO)`,Proteomes )
prot_data <- prot_data %>% as_tibble() %>% dplyr::select(From, Mass, Length )
colnames(prot_data)[1] <- 'Gene_name'
prot_data_no_rep <- prot_data[!duplicated(prot_data$Gene_name), ]
prot_data_no_rep <- prot_data_no_rep[order(prot_data_no_rep$Gene_name) , ]




# install.packages("BBmisc")
library(BBmisc)
prot_length <- as.numeric(prot_data_no_rep$Length)
prot_length <- normalize(prot_length, method = "range", range = c(0, 1))

prot_mass <- as.numeric(prot_data_no_rep$Mass)
prot_mass <- normalize(prot_mass, method = "range", range = c(0, 1))


# combine
prot_property <- cbind(prot_property, prot_length, prot_mass)






## Add lost protein to all with zero
add_name <- setdiff(POI, names(prots))
add_zero_row <- matrix(0, nrow = length(add_name), ncol = ncol(prot_property))
add_zero_row <- as.data.frame(add_zero_row)
rownames(add_zero_row) <- add_name
colnames(add_zero_row) <- colnames(prot_property)

prot_property <- rbind(prot_property, add_zero_row)
dim(prot_property) # 495 by 343

prot_property <- prot_property[order(rownames(prot_property)) , order(colnames(prot_property))]

# save(prot_property, file = "x_prot.RData")



load("x_prot.RData")













# Pairwise similarity ----------------------------------------
psimmat_prot_DTI <- parSeqSim(prot_PPI, cores = 6, type = "local", batches = 20, submat = "BLOSUM62", verbose = TRUE)

psimmat <- psimmat_prot_DTI
rownames(psimmat) <- prot_PPI_names
colnames(psimmat) <- prot_PPI_names

# save(psimmat, file = "psimmat_prot_DTI.RData")

# psimmat <- parSeqSimDisk(prots_data, cores = 6, type = "local", batches = 10, submat = "BLOSUM62", verbose = TRUE)
# save(psimmat, file = "psimmat.RData")


# the amino acid type sanity check and remove the non-standard sequences, To ensure that the protein sequences only have the 20 standard amino acid types which is usually required for the descriptor computation
# remove 1 protein sequences

load("psimmat.RData")
load("prots.RData")

dim(protein_seq_similarity)

protein_seq_similarity <- psimmat
rownames(protein_seq_similarity) <- names(prots) 
colnames(protein_seq_similarity) <- names(prots) 

## Add lost protein to all with zero
add_name <- setdiff(POI, names(prots))
add_zero_col <- matrix(0, nrow = nrow(protein_seq_similarity) , ncol = length(add_name))
colnames(add_zero_col) <- add_name
add_zero_row <- matrix(0, nrow = length(add_name), ncol = ncol(protein_seq_similarity)+dim(add_zero_col)[2])
rownames(add_zero_row) <- add_name

# add_zero_row <- matrix(0, nrow = length(add_name), ncol = ncol(protein_seq_similarity))
protein_seq_similarity <- cbind(protein_seq_similarity, add_zero_col)
protein_seq_similarity <- rbind(protein_seq_similarity, add_zero_row)
dim(protein_seq_similarity)

protein_seq_similarity <- protein_seq_similarity[order(rownames(protein_seq_similarity)) , order(colnames(protein_seq_similarity))]

# save(protein_seq_similarity, file = "protein_seq_similarity.RData")

load("protein_seq_similarity.RData")
protein_seq_similarity
dim(protein_seq_similarity)









# ### Drug data
# 
# Download_data = FALSE
# library("dbparser")
# if (Download_data){
#   # Load dbparser package
#   library(dbparser)
#   # Create SQLite database connection
#   database_connection <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
#   # # DrugBank database sample name
#   # biotech <- "drugbank.xml"
#   # Use DrugBank database in the library
#   read_drugbank_xml_db("drugbank.xml")
#   
#   # Parse all available drug tibbles
#   run_all_parsers(save_table = TRUE, database_connection = database_connection)
#   
#   
#   # List saved tables
#   DBI::dbListTables(database_connection)
#   # Close SQLite connection
#   DBI::dbDisconnect(database_connection)
#   ## load drugs data
#   drugs <- drugs()
#   
#   ## load drug groups data
#   drug_groups <- drug_groups()
#   
#   drug_interaction <- drug_interactions()
#   
#   # drug_sequences()  ### only for biotech drugs
#   
#   drug_salts <- drug_salts()  ### similarity SMILES
#   
#   drug_calc_prop <- drug_calc_prop()
#   
#   
#   save(drugs, file = "drugs.Rdata")
#   save(drug_interaction, file = "drug_interaction.Rdata")
#   save(drug_salts, file = "drug_salts.Rdata")
#   save(drug_calc_prop, file = "drug_calc_prop.Rdata")
# }
# 
# 
# load("drugs.Rdata")
# load("drug_interaction.Rdata")
# load("drug_salts.Rdata")
# load("drug_calc_prop.Rdata")
# druggs <- drugs
# 
# ### filter chemical drugs (not biotech)
# chemical_drugs_index <- which(druggs[["general_information"]][["type"]] == "small molecule")
# chemical_drugs_name <- (drugs[["general_information"]][["primary_key"]])[chemical_drugs_index]
# 
# save(chemical_drugs_name, file = "chemical_drugs_name.Rdata")   
















# ### Drug SMILES similarity
# load("drug_name.Rdata")
# # load("drugs.Rdata")
# # load("drug_interaction.Rdata")
# # load("drug_salts.Rdata")
# load("drug_calc_prop.Rdata")
# load("drug_common.Rdata")
# drug_common
# 
# drug_SMILES_name <- drug_name$`drugbank-id`    ### 258 chemical drugs
# # drug_SMILES_name <- drug_common    ###  36 common drugs
# drug_SMILES <- (drug_calc_prop[which(drug_calc_prop[,1] == "SMILES"),]) [, c(4,2)]
# 
# ### With drug_id,   258 out of 267 drugs have SMILES data
# # drug_SMILES<-drug_SMILES%>%filter(parent_key%in%drug_SMILES_name$`drugbank-id`) ### 258 out of 267 drugs have SMILES data
# drug_SMILES<-drug_SMILES%>%filter(parent_key%in%drug_SMILES_name) ### 
# 
# # save(drug_SMILES, file='drug_SMILES_258.Rdata')
# load('drug_SMILES_258.Rdata')
# drug_SMILES
# 
# SMILES_data <- as.matrix(drug_SMILES[,2])  ## Only SMILES vector
# 
# library(stringdist)
# jaccard_sim <- stringsimmatrix(SMILES_data,SMILES_data, method = "jaccard")  ### Jaccard without padding (no SMILES)
# rownames(jaccard_sim) <- drug_SMILES$parent_key
# colnames(jaccard_sim) <- drug_SMILES$parent_key
# 
# 
# add_name <- setdiff(drug_common, rownames(jaccard_sim))   ### common drugs that not have SMILES
# add_zero_col <- matrix(0, nrow = nrow(jaccard_sim) , ncol = length(add_name)) 
# add_zero_row <- matrix(0, nrow = length(add_name), ncol = ncol(jaccard_sim)+dim(add_zero_col)[2]) 
# 
# jaccard_sim <- cbind(jaccard_sim, add_zero_col)
# jaccard_sim <- rbind(jaccard_sim, add_zero_row)
# 
# rownames(jaccard_sim) <- c(drug_SMILES$parent_key,add_name)
# colnames(jaccard_sim) <- c(drug_SMILES$parent_key,add_name)
# 
# Jaccard_drug_without_SMILES <- jaccard_sim
# # intersect(rownames(jaccard_sim), drug_common)
# # intersect(rownames(DTI), drug_common)
# 
# jaccard_sim <- jaccard_sim[, colnames(jaccard_sim) %in% drug_common]    
# `%!in%` <- Negate(`%in%`)
# 
# ### No need for this line for 36 drugs
# jaccard_sim <- jaccard_sim[rownames(jaccard_sim) %!in% drug_common,] ### 228 drugs verse 36 drug of interest
# 
# jaccard_sim <- jaccard_sim[order(rownames(jaccard_sim)) , order(colnames(jaccard_sim))]
# # save(jaccard_sim, file = "jaccard_sim_known_30_drug.Rdata")  ### 30 drugs of interest
# # save(jaccard_sim, file = "jaccard_sim_known_36_drug.Rdata")  ### 36 drugs of interest
# # save(jaccard_sim, file = "jaccard_sim_known_drug.Rdata")    ### 228 drugs for test




### Drugs without SMILES
jaccard_sim <- Jaccard_drug_without_SMILES
drug_without_SMILES <- setdiff(as.matrix(drug_SMILES_name), rownames(jaccard_sim)) ### 267 - 258 = 9 different

jaccard_sim_all <- rbind(jaccard_sim, matrix(nrow = length(drug_without_SMILES), ncol = dim(jaccard_sim)[2]) )
jaccard_sim_all <- cbind(jaccard_sim_all, matrix(nrow = dim(jaccard_sim_all)[1]  , ncol = length(drug_without_SMILES)  ) )

# jaccard_sim_all[nrow(jaccard_sim): (nrow(jaccard_sim) + length(drug_without_SMILES)),] <- NaN
# jaccard_sim_all[,ncol(jaccard_sim): (ncol(jaccard_sim) + length(drug_without_SMILES))] <- NaN
jaccard_sim_all[nrow(jaccard_sim): (nrow(jaccard_sim) + length(drug_without_SMILES)),] <- 0
jaccard_sim_all[,ncol(jaccard_sim): (ncol(jaccard_sim) + length(drug_without_SMILES))] <- 0
diag(jaccard_sim_all) <- 1

### 267 x 267
rownames(jaccard_sim_all) <- c(drug_SMILES$parent_key,drug_without_SMILES)
colnames(jaccard_sim_all) <- c(drug_SMILES$parent_key,drug_without_SMILES)

head(jaccard_sim_all)

jaccard_sim_all <- jaccard_sim_all[order(rownames(jaccard_sim_all)) , order(colnames(jaccard_sim_all))]


jaccard_sim <- as_tibble(jaccard_sim)
rownames(jaccard_sim) <- drug_SMILES$parent_key


load("DTI_drug_order.Rdata")
DTI_drug_order


# jaccard_sim_DTI <- subset(jaccard_sim, rownames(jaccard_sim) %in% DTI_drug_order  )
# rownames(jaccard_sim_DTI) <- rownames(jaccard_sim)[rownames(jaccard_sim) %in% DTI_drug_order]


# DTI_drug_order[rownames(jaccard_sim) %in% DTI_drug_order]



load("drug_SMILES.Rdata")
load("jaccard_sim_SMILES_only.Rdata")
load("jaccard_sim_all.Rdata")
# save(drug_SMILES, file = "drug_SMILES.Rdata")     ### 258 out of 267 drugs have SMILES data
# save(jaccard_sim, file = "jaccard_sim_SMILES_only.Rdata")    ### 258 x 258 drug Jaccard similarity 
# save(jaccard_sim_all, file = "jaccard_sim_all.Rdata")    ### 267 x 267 drug Jaccard similarity 




# save(jaccard_sim_DTI, file = "jaccard_sim_DTI.Rdata")





# ### DTI (search by target proteins from DrugBank data)
# ###
# # DTI = sif2igraph("stitch_interactions.sif", directed = FALSE)
# # plot(DTI, vertex.size=6)  
# library(dplyr)
# library(ggplot2)
# library(OmnipathR)
# library(igraph)
# library(ggraph)
# library(dbparser)
# library(XML)
# 
# interactions = import_omnipath_interactions() %>% as_tibble()
# # Convert to igraph objects:
# OPI_g = interaction_graph(interactions = interactions )
# 
# get_xml_db_rows("drugbank.xml")
# 
# 
# 
# ## load drugs data
# 
# drugs <- parse_drug() 
# # test <- (drugs %>% as_tibble()) %>% dplyr::select(primary_key, name, type)
# drugs <- drugs %>% as_tibble() %>% dplyr::select(primary_key, name, type)
# # drugs <- rename(drugs, drug_name = name)
# colnames(drugs)[2] <- "drug_name"
# 
# drug_targets_1 <- parse_drug_targets()
# drug_targets_1 <- drug_targets_1[,c(1,2,3,6)]
# ## load polypeptide data
# 
# drug_peptides_1 <- parse_drug_targets_polypeptides()
# drug_peptides_1 <- drug_peptides_1[,c(1,3,6,20)] # id, name, general_function, specific_function, gene_name, parent_id
# # join the 3 datasets
# drug_targets_full <- inner_join(drug_targets_1, drug_peptides_1,
#                                 by=c("id"="parent_id", "name")) # %>% inner_join(drugs, by=c("parent_key"="primary_key"))
# 
# drug_target_inter <- drug_targets_full %>% dplyr::select(parent_key, gene_name)
# 
# save(drug_target_inter, file = "drug_target_inter.Rdata")  ### save all "drug_target_inter" from database










# ### protein of interest (POI)
# load("drug_target_inter.Rdata")
# protein = c("TLR2","ITGB3","CFP","PIK3CA","BAX","CCL5","IL6","IL1B","PTPN11")
# save(protein, file = "protein_POI.Rdata")
# 
# 
# ### DTI
# mylistt<- c()
# for (i in 1:length(protein)  ){
#   index <- which(drug_target_inter$gene_name == protein[i])
#   mylistt <- c(mylistt, index)
# }
# drug_target_inter_data <- drug_target_inter[mylistt,]   ### filter with given POI
# 
# load("DTI_tibble.Rdata")
# DTI_tibble
# rownames(DTI_tibble)
# 
# # drug_nam <- drug_target_inter_data$parent_key
# drug_nam <- rownames(DTI_tibble)
# 
# 
# load("chemical_drugs_name.Rdata")  ### Load chemical drug database
# drug_common <- intersect(drug_nam, chemical_drugs_name)
# 
# 
# # length(drug_nam)    # 53 drugs (36 chemical drugs) for 7 targets
# length(drug_common) # 30 drugs in common (manual & database)
# 
# ### 40 interactions for DTI (7 targets and 36 drugs)
# mylistt<- c()
# for (i in 1:length(drug_common)  ){
#   index <- which(drug_target_inter_data$parent_key == drug_common[i])
#   mylistt <- c(mylistt, index)
# }
# drug_target_inter_data <- drug_target_inter_data[mylistt,]  ### Drugs that have interactions with POI
# 
# 
# # test <- filter(drug_target_inter_data, parent_key %in% as_ids(vertex_DDI) )
# 
# 
# unique(drug_target_inter_data$parent_key)
# 
# 
# 
# save(drug_common, file = 'drug_common.Rdata')   
# save(drug_target_inter_data, file = 'drug_target_inter_data.Rdata')   






# ########## Drug pathway
# load("drug_common.Rdata")
# # drugnames = drug_SMILES[,1]$parent_key          ### Drug of interest
# drug_names <- drug_common          ### Drug of interest
# # drug_name$`drugbank-id`
# drug_id <- drug_common
# 
# mylis<- c()
# for(i in 1:length(drug_id)){
#   index <- which(drugs[,1]==drug_id[i])
#   mylis <- c(mylis,index)
# }
# drug_names = drugs[mylis,]
# 
# ### given drugs of interest, find DTI pathway
# drug_target_data_sample <- drug_targets_full %>% filter(organism == "Humans",parent_key %in% drug_names$primary_key)
# 
# drug_targets <-  drug_target_data_sample %>% dplyr::select(-organism) %>% mutate(in_OP = id %in% c(interactions$source))
# # not all drug-targets are in OP.
# print(all(drug_targets$in_OP))
# 
# drug_targets %>% group_by(parent_key) %>% summarise(any(in_OP))
# 
# ### protein of interest (POI)
# POI = tibble(protein = c("TLR2","ITGB3","CFP","PIK3CA","BAX","CCL5","IL6","IL1B","PTPN11") )
# 
# POI <- POI %>% mutate(in_OP = protein %in% interactions$target_genesymbol)
# # all POI is in Omnipath
# print(all(POI$in_OP))
# 
# ### Drug name
# drug_of_interest = "Abciximab"
# source_nodes <- drug_targets %>%
#   filter(in_OP, drug_name==drug_of_interest) %>%
#   pull(gene_name)
# target_nodes <- POI %>% filter(in_OP) %>% pull(protein)
# 
# collected_path_nodes = list()
# 
# for(i_source in 1:length(source_nodes)){
#   
#   paths <- shortest_paths(OPI_g, from = source_nodes[i_source],
#                           to = target_nodes,
#                           output = 'vpath')
#   path_nodes <- lapply(paths$vpath,names) %>% unlist() %>% unique()
#   collected_path_nodes[[i_source]] <- path_nodes
# }
# collected_path_nodes <- unlist(collected_path_nodes) %>% unique()
# 
# 
# 
# cisplatin_nodes <- c(source_nodes,target_nodes, collected_path_nodes) %>%
#   unique()
# cisplatin_network <- induced_subgraph(graph = OPI_g,vids = cisplatin_nodes)
# 
# 
# V(cisplatin_network)$node_type = ifelse(
#   V(cisplatin_network)$name %in% source_nodes, "direct drug target",
#   ifelse(
#     V(cisplatin_network)$name %in% target_nodes,"Protein of Interest","intermediate node"))
# 
# ggraph(
#   cisplatin_network,
#   layout = "lgl",
#   area = vcount(cisplatin_network)^2.3,
#   repulserad = vcount(cisplatin_network)^1.2,
#   coolexp = 1.1
# ) +
#   geom_edge_link(
#     aes(
#       start_cap = label_rect(node1.name),
#       end_cap = label_rect(node2.name)),
#     arrow = arrow(length = unit(4, 'mm')
#     ),
#     edge_width = .5,
#     edge_alpha = .2
#   ) +
#   geom_node_point() +
#   geom_node_label(aes(label = name, color = node_type)) +
#   scale_color_discrete(
#     guide = guide_legend(title = 'Node type')
#   ) +
#   theme_bw() +
#   xlab("") +
#   ylab("") +
#   ggtitle("Abciximab induced network")
# 






