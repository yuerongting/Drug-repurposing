# BiocManager::install("impute")

library("protr")
# setwd('c:\\Users\\yrt05\\Desktop\\Covid_coference_presentation')  # \\ abs path

# # load FASTA files
# extracell <- readFASTA("uniprot.fasta")
# 
# extracell_test <- readFASTA("uniprot-.fasta")


library("readxl")

### protein identifier (map 667 gene) (PPI) 

load('map_gene_prot.RData') ### gene identifier mapping
map_gene_prot

library(influential)
PPI_Network_1 = sif2igraph("string_interactions_short.tsv.sif", directed = FALSE)
DEG_names <- as_ids(V(PPI_Network_1)) ### PPI proteins 667


prot_id <- filter(map_gene_prot, gene_name %in% DEG_names)
prot_id <- unique(prot_id)



### protein sequence
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

test <- listAttributes(ensembl)
gene_to_protein_mapping <- getBM(attributes = c("external_gene_name", "uniprotswissprot"),
                                 filters = "external_gene_name",
                                 values = DEG_names,
                                 mart = ensembl)

prot_id_map <- gene_to_protein_mapping %>% filter(!is.na(uniprotswissprot) & uniprotswissprot != "")

prot_id_map <- filter(prot_id_map , external_gene_name %in% unique(prot_id_map$external_gene_name))
# 
# test_gene <- prot_id_map$external_gene_name
# test_gene[duplicated(prot_id_map$external_gene_name)]

unique(prot_id_map$external_gene_name)  # 661 genes mapped with prot ID
unique(prot_id_map$uniprotswissprot) # 665 prot id

prot_id_map <- prot_id_map[!duplicated(prot_id_map$external_gene_name),]

prot_id_map <- prot_id_map[!duplicated(prot_id_map$uniprotswissprot),]




nrow(unique(prot_id_map))  ### 660 prot


setdiff(DEG_names, unique(prot_id_map$external_gene_name))


# write.table(prot_id_map$id, file = "prot_id.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




# BiocManager::install("UniProt.ws")
library(UniProt.ws)

# up <- UniProt.ws(taxId = 9606)  # 9606 is the tax ID for Homo sapiens

library(httr)
library(jsonlite)
library(stringr)

# uniprot_id = uniprot_ids[1]

uniprot_ids <- prot_id_map$uniprotswissprot

length(unique(uniprot_ids))

get_protein_sequence <- function(prot_id) {
  base_url <- "https://www.uniprot.org/uniprot/"
  response <- GET(paste0(base_url, prot_id, ".fasta"))
  if (status_code(response) == 200) {
    fasta_content <- content(response, "text")
    # Extract sequence part from FASTA format
    sequence <- str_split(fasta_content, "\n")[[1]]
    sequence <- paste(sequence[-1], collapse = "")
    return(sequence)
  } else {
    return(NA)
  }
}
protein_sequences <- sapply( uniprot_ids , get_protein_sequence)



names(protein_sequences) <- prot_id_map$external_gene_name
test <- list(protein_sequences)
test <- setNames(as.list(protein_sequences), names(protein_sequences))
protein_sequences <- test
prots <- protein_sequences ## 660 prot sequence

protein_sequences[i]

# save(prots, file = "prots.RData")



load("prots.RData")
prots

unique(prots)


# unique(names(prot_PPI))



# Protein sequence feature------------------------------

prot_PPI <- prots
prot_PPI_names <- names(prot_PPI)

# Conjoint Triad  343  dim <<Predicting Protein-protein Interactions Based Only on Sequences Information>>
library(protr)


prot_PPI_names[267] ### DIO2
index_error = 267

prot_property = matrix(, nrow = length(prot_PPI_names), ncol = length(names(extractCTriad(prot_PPI$`H3-4`)))   )
i=1
for (seq in prot_PPI){
  # if (i==267){
  #   next
  # } ### i = 269: x has unrecognized amino acid type
  if(i!= index_error){
    print(i)
    
    encode_343d <- as.numeric(   extractCTriad(seq)   )
    
    prot_property[i,] <- encode_343d

  }
  i = i+1
}


prot_property <- as.data.frame(prot_property)
colnames(prot_property) <- names(extractCTriad(prot_PPI$`H3-4`))
row.names(prot_property) <- prot_PPI_names
# test <- prot_property[-267,]

prot_property <- prot_property[order(rownames(prot_property)) , order(colnames(prot_property))]

prot_property <- prot_property[complete.cases(prot_property), ]


#### updated DEG_total
DEG_updated <- rownames(prot_property)

# write(DEG_updated,"DEG_total_filtered.txt",sep = '')  # updated version from "DEG_GSE235210.R"



# # additional properties (weight, mass, length)
# protein_data <- read_excel("uniprot-compressed_10_10.xlsx")
# prot_data <- filter(protein_data, `From` %in% gene_name)
# # prot_data <- prot_data %>% as_tibble() %>% dplyr::select(From, Mass, Features, Length, Sequence, `Gene Ontology (GO)`,Proteomes )
# prot_data <- prot_data %>% as_tibble() %>% dplyr::select(From, Mass, Length )
# colnames(prot_data)[1] <- 'Gene_name'
# prot_data_no_rep <- prot_data[!duplicated(prot_data$Gene_name), ]
# prot_data_no_rep <- prot_data_no_rep[order(prot_data_no_rep$Gene_name) , ]
# 
# # install.packages("BBmisc")
# library(BBmisc)
# prot_length <- as.numeric(prot_data_no_rep$Length)
# prot_length <- normalize(prot_length, method = "range", range = c(0, 1))
# 
# prot_mass <- as.numeric(prot_data_no_rep$Mass)
# prot_mass <- normalize(prot_mass, method = "range", range = c(0, 1))
# 
# 
# # combine
# prot_property <- cbind(prot_property, prot_length, prot_mass)
# 
# 
# ## Add lost protein to all with zero
# add_name <- setdiff(POI, names(prots))
# add_zero_row <- matrix(0, nrow = length(add_name), ncol = ncol(prot_property))
# add_zero_row <- as.data.frame(add_zero_row)
# rownames(add_zero_row) <- add_name
# colnames(add_zero_row) <- colnames(prot_property)
# 
# prot_property <- rbind(prot_property, add_zero_row)
# dim(prot_property) # 495 by 343
# 
# prot_property <- prot_property[order(rownames(prot_property)) , order(colnames(prot_property))]


# save(prot_property, file = "x_prot.RData")
load("x_prot.RData") ### 659 proteins




prot_PPI_names[267] ### DIO2

if ('DIO2' %in% names(prot_PPI)) {
  # Remove the element from the list
  prot_PPI[['DIO2']] <- NULL
}

# prot_PPI <- 








# Pairwise similarity ----------------------------------------
psimmat_prot_DTI <- parSeqSim(prot_PPI, cores = 6, type = "local", batches = 20, submat = "BLOSUM62", verbose = TRUE)

psimmat <- psimmat_prot_DTI

prot_PPI_names <- prot_PPI_names[-index_error]
rownames(psimmat) <- prot_PPI_names
colnames(psimmat) <- prot_PPI_names

# save(psimmat, file = "psimmat.RData")


# the amino acid type sanity check and remove the non-standard sequences, To ensure that the protein sequences only have the 20 standard amino acid types which is usually required for the descriptor computation
# remove 1 protein sequences

load("psimmat.RData")
load("prots.RData")

protein_seq_similarity <- psimmat

prot_names <- names(prots) 
prot_names <- prot_names[-index_error]


rownames(protein_seq_similarity) <- prot_names
colnames(protein_seq_similarity) <- prot_names

protein_seq_similarity <- protein_seq_similarity[order(rownames(protein_seq_similarity)) , order(colnames(protein_seq_similarity))]


dim(protein_seq_similarity)

# save(protein_seq_similarity, file = "protein_seq_similarity.RData")



# ## Add lost protein to all with zero
# add_name <- setdiff(POI, prot_names)
# add_zero_col <- matrix(0, nrow = nrow(protein_seq_similarity) , ncol = length(add_name))
# colnames(add_zero_col) <- add_name
# add_zero_row <- matrix(0, nrow = length(add_name), ncol = ncol(protein_seq_similarity)+dim(add_zero_col)[2])
# rownames(add_zero_row) <- add_name
# # add_zero_row <- matrix(0, nrow = length(add_name), ncol = ncol(protein_seq_similarity))
# protein_seq_similarity <- cbind(protein_seq_similarity, add_zero_col)
# protein_seq_similarity <- rbind(protein_seq_similarity, add_zero_row)
# dim(protein_seq_similarity)
# 
# protein_seq_similarity <- protein_seq_similarity[order(rownames(protein_seq_similarity)) , order(colnames(protein_seq_similarity))]
# 
# # save(protein_seq_similarity, file = "protein_seq_similarity.RData")
# 
# load("protein_seq_similarity.RData")
# protein_seq_similarity
# dim(protein_seq_similarity)









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


















