library(dplyr)
library(ggplot2)
library(OmnipathR)
library(igraph)
library(ggraph)
library(dbparser)
library(XML)
# getwd()

### Drug SMILES similarity

load("drug_name.Rdata")  ### All 1099 drugs
drug_name
# load("drugs.Rdata")
# load("drug_salts.Rdata")

# load("drug_interact.Rdata")
load("drug_calc_prop.Rdata")


### All SMILES 11924
drug_SMILES <- (drug_calc_prop[which(drug_calc_prop[,1] == "SMILES"),]) [, c(4,2)]  
drug_SMILES_all <- drug_SMILES

# save(drug_SMILES_all, file = 'drug_SMILES.Rdata') # All SMILEs with DDIs




# drug_SMILES_name <- unique(drug_name$`drugbank-id`)    ### 468 chemical drugs in DDI

load('drug_name_DDI.Rdata')

drug_name_DDI
# drug_SMILES_name <- unique(drug_name$`drugbank-id`)    ### 468 chemical drugs in DDI

drug_SMILES_name <- drug_name_DDI  ### 3410 drugs in DDI

# drug_SMILES<-drug_SMILES%>%filter(parent_key%in%drug_SMILES_name$`drugbank-id`) ### 258 out of 267 drugs have SMILES data
drug_SMILES_DDI <- drug_SMILES%>%filter(parent_key%in%drug_SMILES_name)

# save(drug_SMILES_DDI, file='drug_SMILES_DDI.Rdata')  






# load('drug_SMILES_468.Rdata')
# drug_SMILES
# drug_SMILES_468 = drug_SMILES

load('drug_SMILES.Rdata')
drug_SMILES_all ### All drug SMILES
drug_SMILES <- drug_SMILES_all


# load('DDI_all.Rdata')
load("DDI.Rdata")   ### 3410 drugs, 504534 edges
DDI_node <- as_ids(V(DDI))




### SMILES similarity

load('DTI_graph.Rdata')
DTI_graph

load('drug_target_inter_data.Rdata')
drug_target_inter_data


load('drug_SMILES_DDI.Rdata')
drug_SMILES_DDI

drug_SMILES <- drug_SMILES_DDI

SMILES_data <- as.matrix(drug_SMILES[,2])  ## Only SMILES vector

library(stringdist)
jaccard_sim <- stringsimmatrix(SMILES_data,SMILES_data, method = "jaccard")  ### Jaccard without padding (no SMILES)
rownames(jaccard_sim) <- drug_SMILES$parent_key
colnames(jaccard_sim) <- drug_SMILES$parent_key


### Drugs with no SMILE data

load('drug_common.Rdata')
drug_common  ### 703 drugs from DTI






### Column only drugs in DTI
# jaccard_sim <- jaccard_sim[, colnames(jaccard_sim) %in% drug_common]  ### 3490 x 214
dim(jaccard_sim)  



# ### 3410 drugs in DDI have 353 drugs overlapped in DTI
# jaccard_sim_DTI_drugs <- jaccard_sim[rownames(jaccard_sim) %in% drug_common,] ### 214 drugs
# dim(jaccard_sim_DTI_drugs)
# 
# # order by names
# jaccard_sim_DTI_drugs <- jaccard_sim_DTI_drugs[order(rownames(jaccard_sim_DTI_drugs)) , order(colnames(jaccard_sim_DTI_drugs))]
# jaccard_sim_DTI_drugs <- as.data.frame(jaccard_sim_DTI_drugs)
# dim(jaccard_sim_DTI_drugs)  



### 3410 Drugs in  DDI
jaccard_sim_468_DDI <- jaccard_sim[rownames(jaccard_sim) %in% drug_SMILES_DDI$parent_key,] 
# jaccard_sim_468_DDI <- jaccard_sim_468_DDI[rownames(jaccard_sim_468_DDI) %!in% drug_common,] 
jaccard_sim_468_DDI <- jaccard_sim_468_DDI[order(rownames(jaccard_sim_468_DDI)) , order(colnames(jaccard_sim_DTI_drugs))]
jaccard_sim_468_DDI <- as.data.frame(jaccard_sim_468_DDI)
dim(jaccard_sim_468_DDI)   ### 3410 x 3410


# save(jaccard_sim_468_DDI, file = "jaccard_sim_468_DDI.Rdata")

# ### Other drugs
# `%!in%` <- Negate(`%in%`)
# jaccard_sim_other <- jaccard_sim[rownames(jaccard_sim) %!in% union(drug_common, drug_SMILES_468$parent_key),] ### 266 drugs verse 30 drug of interest
# jaccard_sim_other <- jaccard_sim_other[order(rownames(jaccard_sim_other)) , order(colnames(jaccard_sim_other))]
# jaccard_sim_other <- as.data.frame(jaccard_sim_other)
# dim(jaccard_sim_other)    ### 3022 x 214




# intersect(row.names(jaccard_sim) , rownames(DTI))
# intersect(drug_common, drug_SMILES_name)
# intersect(intersect(drug_common, drug_SMILES_name), rownames(DTI))

# save(jaccard_sim_DTI_drugs, file = "jaccard_sim_DTI_drug.Rdata")  ### 214 to 214 DTI drugs
# save(jaccard_sim_468_DDI, file = "jaccard_sim_label_0.Rdata")  ### 468 DDI drugs (training, with labels)
# save(jaccard_sim_other, file = "jaccard_sim_all_other_to_all.Rdata")  ### Other drugs (testing)






# 