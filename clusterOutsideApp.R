# Clustering of data  - script
library(tidyverse)

load('24Feb23_cell_level_neighborhood_data20cells.RData')
filename <- "24Feb23_7_LMA_types_20_cells.csv"

# Get summary statistics per neighborhood
community_summary <- lapply(all_neighborhoods, function(x) {
  summary_no_location <- x %>% filter(!is.na(Neighborhood)) %>% group_by(Neighborhood, Phenotype) %>% summarise(sum = n()) %>%
    pivot_wider(names_from = "Phenotype", values_from = "sum")
  
  summary_no_location[is.na(summary_no_location)] <- 0

  summary_no_location <- summary_no_location %>% filter(Neighborhood != "Free_cell")
})
  
# Convert to data frame
community_summary <-data.table::rbindlist(community_summary, fill=TRUE, idcol="Image")
community_summary[is.na(community_summary)] <- 0

# Convert to matrix
communities_matrix <- community_summary %>% select(-Neighborhood, -Image) %>%
  as.matrix()

# Convert to proportions
communities_matrix_prop <- t(apply(communities_matrix,1, function(x) x/sum(x)))

# Calculate dissimilarity and build dendrogram
hier_bray <- vegan::vegdist(communities_matrix_prop, method = "bray")
tree <- hclust(hier_bray, method="ward.D2")
tree_dend <- as.dendrogram(tree)

# Cut tree and cluster
clust <- cutree(tree_dend, k=7)
neighborhood_df_with_LMA <- community_summary %>% mutate(LMA = clust)

write_csv(neighborhood_df_with_LMA, filename)
