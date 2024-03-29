# Make nice plots of neighborhoods per image and overlay them on original RGB images
# 14Feb2023

library(tidyverse)
library(ggsci)

# load cell level data with neighborhoods
load("24Feb23_cell_level_neighborhood_data20cells.RData")

# Import Cluster to LMA key
neighToLMAKey <- vroom::vroom("clusteringResults/24Feb23_7_LMA_types_20_cells.csv")

# The X and Y positions are wrong - need to update them by matching by Cell ID with original cell level data
cell_level_data <- vroom::vroom("23Feb23_cellLevelData.csv")
cell_level_data <- cell_level_data %>% group_by(Image) %>%
  mutate(Cell.ID = paste0("Cell_", row_number()))

# Loop over all images in the list to generate plot
for(image in 1:length(all_neighborhoods)) {
roi <- all_neighborhoods[[image]]
roi$Image <- names(all_neighborhoods)[image] 

# Update location by Cell ID
roi <- roi %>% select(-Cell.X.Position, -Cell.Y.Position) %>%
  left_join(., cell_level_data)

# Invert Y-axis
roi$Centroid_y <- max(roi$Centroid_y) - roi$Centroid_y

#Update to LMA ID
LMAKey <- neighToLMAKey %>% filter(Image == roi$Image[1]) %>%
  group_by(Neighborhood) %>%
  summarise(LMA = unique(LMA))

roi <- left_join(roi, LMAKey)

# Get points that make the convex-hull of each cluster
roi_hull <- roi %>% filter(!is.na(Neighborhood) & Neighborhood != "Free_cell") %>%
  group_by(Neighborhood) %>%
  slice(chull(Centroid_x, Centroid_y))

# Get subset of cells that are in neighborhood
roi_not_negative <- roi %>% filter(Phenotype != "Negative")

# Plot
roi %>% ggplot()+
  geom_point(data=roi, aes(x=Centroid_x, y=Centroid_y), size=0.1, alpha=0.1, color="grey")+
  geom_point(data=roi_not_negative, aes(x=Centroid_x, y=Centroid_y, color=Phenotype), size=0.3)+
  geom_polygon(data=roi_hull, aes(x=Centroid_x, y=Centroid_y, group=Neighborhood, fill=factor(LMA)), alpha=0.5, show.legend = TRUE)+
  scale_color_manual(values = c("black", "#afcbe3","#377eb8","#f4a3a5",
                                          "#e41a1d", "#4DAF4A", "#FFFF33", "#800080", "#c080c0"))+
  scale_fill_npg(name="LMA")+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme_void()+
  theme(panel.background = element_rect(fill='white', colour='white'),
        plot.background = element_rect(fill = "white"))

ggsave(paste0(names(all_neighborhoods)[image], ".png"),
       width=14.70,height=8.38,units='in')
}
             