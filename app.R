library(shiny)
library(SPIAT)
library(tidyverse)
library(ggsci)

options(repos = BiocManager::repositories())

images <- tibble(n = c(1:4),
                 filename = c("HGS_B_4_B+T_Scan1_Core[2,3,J]_[19600,46803]_component_data.tif",
                              "HGS_B_4_B+T_Scan1_Core[4,2,E]_[14543,54060]_component_data.tif",
                              "HGS_D_4_B+T_Scan1_Core[2,2,B]_[13421,44528]_component_data.tif",
                              "HGS_D_4_B+T_Scan1_Core[4,4,C]_[14067,54839]_component_data.tif"))
ihc_data <- read.csv("ihc_data.csv.gz")

ui <- fluidPage(titlePanel("SPIAT package demo"),
                sidebarLayout(
  sidebarPanel(selectInput("image", "Select image", choices = setNames(images$filename, images$n)),
               sliderInput("n_size", "Minimum neighborhood size", value = 12, min = 0, max = 30),
               sliderInput("radius", "Maximum cell-cell distance", value = 15, min = 0, max = 30),
               submitButton("Update View", icon("refresh"))),
  mainPanel(imageOutput("selectedImage"),
            plotOutput("neighborhood_map", width = 400, height=400))))

server <- function(input, output) {

  single_image_data <- reactive({
    ihc_data %>% filter(Image == input$image) %>%
    mutate(ID = row_number())
    })


  spe_object <- reactive({
    coord_x = single_image_data()$Centroid_x
    coord_y = single_image_data()$Centroid_y
    phenotypes = single_image_data()$Class
    dummy_intensity = rep(0, nrow(single_image_data()))
    intensity_matrix = matrix(dummy_intensity, nrow=1, ncol=nrow(single_image_data()))
    colnames(intensity_matrix) = single_image_data()$ID
    spe_object = format_image_to_spe(intensity_matrix=intensity_matrix,
                                   phenotypes = phenotypes, coord_x = coord_x,coord_y = coord_y)
  })


output$selectedImage <- renderImage({
  list(
    src = file.path("images", paste0(input$image, ".png")),
    contentType = "image/jpeg",
    width = 400,
    height = 400
    )
  }, deleteFile = FALSE)

output$neighborhood_map <- renderPlot({
  plot1 <- identify_neighborhoods(spe_object(),
                                  cell_types_of_interest = c("CD20p", "CD3pCD8p", "CD3pCD8n", "CD79ApCD20n", "CD3pCD8nFoxP3p",
                                                             "CD3pCD8pFoxP3p"),
                                  method="hierarchical",
                                  radius=input$radius,
                                  min_neighborhood_size = input$n_size,
                                  feature_colname = "Phenotype")

  coordinates <- data.frame(spatialCoords(plot1)) %>%
    mutate(Cell.ID = row_number())

  coordinates$Cell.Y.Position <- max(coordinates$Cell.X.Position) - coordinates$Cell.Y.Position

  phenotypes <- colData(plot1) %>%
    as.data.frame() %>%
    mutate(Cell.ID = as.numeric(Cell.ID)) %>%
    arrange(Cell.ID)

  coordinates <- left_join(coordinates, phenotypes)
  in_neighborhoods <- subset(coordinates, Neighborhood != "Free_cell" & !is.na(Neighborhood))

  label_location <- vector()
  for (Clusternumber in unique(in_neighborhoods$Neighborhood)) {
    cells_in_Cluster <- in_neighborhoods[in_neighborhoods$Neighborhood ==
                                           Clusternumber, ]
    minX <- min(cells_in_Cluster$Cell.X.Position)
    maxX <- max(cells_in_Cluster$Cell.X.Position)
    minY <- min(cells_in_Cluster$Cell.Y.Position)
    maxY <- max(cells_in_Cluster$Cell.Y.Position)
    averageX <- (minX + maxX)/2
    averageY <- (minY + maxY)/2
    label_location <- rbind(label_location, c(Clusternumber,
                                              averageX, averageY))
  }

  label_location <- as.data.frame(label_location)
  
  label_location[,2:3] <- sapply(label_location[,2:3], as.numeric)
  colnames(label_location) <- c("Cluster", "Cell.X.Position", "Cell.Y.Position")
  label_location$Cluster <- gsub("Cluster_", "", label_location$Cluster)


  ggplot(coordinates, aes(Cell.X.Position, Cell.Y.Position))+
    geom_point(alpha = 0.2)+
    geom_point(data=in_neighborhoods, aes(color=Neighborhood))+
    geom_text(data=label_location, aes(x=Cell.X.Position, y=Cell.Y.Position,
                                       label=Cluster), size=14)+
    scale_color_npg()+
    theme_void()+
    theme(legend.position="none")
})

output$neighborhood_composition <- renderPlot({plot(1,1)})

}


shinyApp(ui = ui, server = server)
