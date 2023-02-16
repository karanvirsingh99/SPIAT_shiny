library(shiny)
library(SPIAT)
library(tidyverse)
library(ggsci)
library(vegan)
library(dendextend)
library(colourpicker)

options(repos = BiocManager::repositories())
options(shiny.maxRequestSize = 1000 * 1024^2)
options(dplyr.summarise.inform = FALSE)


#################################################################
##                          Functions                          ##
#################################################################

find_communities_list <- function(spe_object, radius, neighborhood_size, cell_types){
  
  error_catch <- tryCatch({
    identify_neighborhoods(spe_object = spe_object,
                           cell_types_of_interest = cell_types,
                           method="hierarchical",
                           radius=radius,
                           min_neighborhood_size = neighborhood_size,
                           feature_colname = "Phenotype")
    1
  }, error=function(e) 2)
  
  # error_catch == 1 when no error 
  if(error_catch == 1){
    
    communities <- identify_neighborhoods(spe_object = spe_object,
                                          cell_types_of_interest = cell_types,
                                          method="hierarchical",
                                          radius=radius,
                                          min_neighborhood_size = neighborhood_size,
                                          feature_colname = "Phenotype")
    
    communities_as_df <- colData(communities) %>% as.data.frame()
    communities_as_df <- cbind(communities_as_df, spatialCoords(communities))
    
    return(communities_as_df)
    
  } else(
    return(NULL) # return empty dataframe if no error
  )
} 

#################################################################
##                         Sample data                         ##
#################################################################

images <- tibble(n = c(1:4),
                 filename = c("HGS_B_4_B+T_Scan1_Core[2,3,J]_[19600,46803]_component_data.tif",
                              "HGS_B_4_B+T_Scan1_Core[4,2,E]_[14543,54060]_component_data.tif",
                              "HGS_D_4_B+T_Scan1_Core[2,2,B]_[13421,44528]_component_data.tif",
                              "HGS_D_4_B+T_Scan1_Core[4,4,C]_[14067,54839]_component_data.tif"))
ihc_data <- read.csv("ihc_data.csv.gz")

#################################################################
##                        UI function                          ##
#################################################################

##------------------------------------------------------------------------------------------------------
##  The three features in my app are                                                                                                
##  1) a slider that lets you choose the minimimum neighborhood size 
##    (useful in visualizing output with different thresholds)                                                                
##  2) a slider that lets you chose the maximum distance between cells to be considered interacting 
##    (useful in visualizing output with different thresholds)                            
##  3) a dropdown image selector, that allows you to run the same settings on different sample 
##    images and their corresponding data  
##------------------------------------------------------------------------------------------------------

ui <- fluidPage(
  tabsetPanel(
  tabPanel("Demo",
           titlePanel("SPIAT package demo"), #Title 
           sidebarLayout(
             #Sidebar Panel
             sidebarPanel(selectInput("image", "Select image", choices = setNames(images$filename, images$n)), #Image selection
                          sliderInput("n_size", "Minimum neighborhood size", value = 12, min = 0, max = 30), #Neighborhood size input
                          sliderInput("radius", "Maximum cell-cell distance", value = 15, min = 0, max = 30) #Radius input
                          # submitButton("Update View", icon("refresh"))
             ), #Update view button
             #Main Panel
             mainPanel(imageOutput("selectedImage"), #Selected image output
                       plotOutput("neighborhood_map", width = 400, height=400) # Image output
             )
           )),
  tabPanel("Start here",
           sidebarLayout(sidebarPanel(
             tabsetPanel(
               tabPanel("Analyze",
               fileInput("upload", "Data", buttonLabel = "Browse"),
               checkboxGroupInput("cell_types", "Phenotypes", ""),
               numericInput("final_radius", "Maximimum cell-cell distance", 15),
               numericInput("min_neighborhood_size", "Minimum neighborhood size", 15),
               actionButton("run_detection", "Run neighborhood detection")),
               tabPanel("Plot",
               numericInput("number_of_clusters","Number of clusters", 2),
               textInput("colors", "Enter colors for the barplot", value = pal_npg()(10)),
               downloadButton("download", "Download .csv")
             ))),
             mainPanel(tabsetPanel(
               tabPanel("Instructions",
                       h1("SPIAT workflow"),
                       h4("Step 1. Upload a csv file using the Browse button"),
                       ("Every row should be a cell. It should contain columns Image, Class, Centroid_x, Centroid_y"),
                       h4("Step 2. After uploading the data, the app will read it into R, and convert each image into a spatial object. This takes some time!"),
                       ("Two lines of text will appear here to notify you when that is done"),
                       br(),
                       br(),
                       textOutput("data_read_in"),
                       textOutput("list_loop"),
                       h4("Step 3. Once the data is ready, the phenotypes (from the Class column of the file) will pop up on the left"),
                       ("Choose which cell types you want to use for the neighborhood detection"),
                       h4("Step 4. Choose the max cell-cell distance and neighborhood size and press 'Run neighborhood detection'"),
                       ("The app will start detecting neighborhoods in every image. This may take a while depending on the size of the data"),
                       h4("A dendrogram and barplots will pop up in the Plots")),
               tabPanel("Plots",
                       plotOutput("dendrogram"),
                       plotOutput("bar_plot")),
               tabPanel("Summary",
                        tableOutput("summary_stats"),
                        ("This table shows the average proportion of each cell type in each neighborhood type")
                        
               )))
           ))))



#################################################################
##                       Server function                       ##
#################################################################

server <- function(input, output, session) {
  
  #################################################################
  ##                            TAB 1                            ##
  #################################################################
  
  # Filter sample data based on which image is chosen
  single_image_data <- reactive({
    ihc_data %>% filter(Image == input$image) %>%
      mutate(ID = row_number())
  })
  
  
  # Create a spatial object based on the sample data from chosen image
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
  
  # Render image based on which image is chosen
  output$selectedImage <- renderImage({
    list(
      src = file.path("images", paste0(input$image, ".png")),
      contentType = "image/jpeg",
      width = 400,
      height = 400
    )
  }, deleteFile = FALSE)
  
  # Output neighborhood map based on image chosen, neighborhood radius and minimum neighborhood size
  
  # Step 1 - use the SPIAT::identify_neighborhoods function to create a dataframe that tells me which cells are in which neighborhood
  output$neighborhood_map <- renderPlot({
    plot1 <- identify_neighborhoods(spe_object(),
                                    cell_types_of_interest = c("CD20p", "CD3pCD8p", "CD3pCD8n", "CD79ApCD20n", "CD3pCD8nFoxP3p",
                                                               "CD3pCD8pFoxP3p"),
                                    method="hierarchical",
                                    radius=input$radius,
                                    min_neighborhood_size = input$n_size,
                                    feature_colname = "Phenotype")
    
    # Step 2 - extract the cell coordinates from the neighborhood output
    coordinates <- data.frame(spatialCoords(plot1)) %>%
      mutate(Cell.ID = row_number())
    
    # Step 3 - invert y-axis so that it matches the image (just a QuPath coordinates quirk)
    coordinates$Cell.Y.Position <- max(coordinates$Cell.X.Position) - coordinates$Cell.Y.Position
    
    # Step 4 - extract the phenotypes from the neighborhood output (this also has the neighborhood ID)
    phenotypes <- colData(plot1) %>%
      as.data.frame() %>%
      mutate(Cell.ID = as.numeric(Cell.ID)) %>%
      arrange(Cell.ID)
    
    # Join coordinates and neighborhood id
    coordinates <- left_join(coordinates, phenotypes)
    
    # I want to colour only the cells that are in a neighborhood, so I subset the previous df
    in_neighborhoods <- subset(coordinates, Neighborhood != "Free_cell" & !is.na(Neighborhood))
    
    # Calculate where to plot the label - this is part of the source code for identify_neighborhood function
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
    
    # Clean up the label locations (convert to numeric, add column names)
    label_location[,2:3] <- sapply(label_location[,2:3], as.numeric)
    colnames(label_location) <- c("Cluster", "Cell.X.Position", "Cell.Y.Position")
    label_location$Cluster <- gsub("Cluster_", "", label_location$Cluster)
    
    #Make a colour vector (some colours are repeated but that is A-okay!)
    colors = rep(ggsci::pal_npg()(10), 100)[1:length(label_location$Cluster)]
    
    # Plot the cells, colored by neighborhood with text label with neighborhood ID
    ggplot(coordinates, aes(Cell.X.Position, Cell.Y.Position))+
      geom_point(alpha = 0.2)+
      geom_point(data=in_neighborhoods, aes(color=Neighborhood))+
      geom_text(data=label_location, aes(x=Cell.X.Position, y=Cell.Y.Position,
                                         label=Cluster), size=14)+
      scale_color_manual(values=colors)+
      theme_void()+
      theme(legend.position="none")
  })
  
  
  #################################################################
  ##                            TAB 2                            ##
  #################################################################
  
  cell_level <- reactive({
    req(input$upload)
    # From Mastering Shiny (David Heidl)
    table_in <- read_csv(input$upload$datapath)
    updateCheckboxGroupInput(session = session,inputId = "cell_types", label = "Phenotypes", choices = sort(unique(table_in$Class)),
                             selected = sort(unique(table_in$Class)))
    return(table_in)
  })
  
  output$data_read_in <- renderText({
    req(input$upload)
    paste("File", input$upload$name,
          "successfully read in. There are ", nrow(cell_level()), "rows and", ncol(cell_level()), "columns!")
  })
  
  observeEvent(input$cell_types,
               {print(input$cell_types)})
  
  ##----------------------------------------------------------------
  ##                      Convert data to list                     -
  ##----------------------------------------------------------------
  
  images_as_list <- reactive({
    req(input$upload)
    imgDataList <- list() #empty list
    image_names <- unique(cell_level()$Image) #get all image names
    
    for (img in image_names) { 
      print(paste("Processing image", img))
      # Get all cells corresponding to image i and add them to imgDataList
      imgData <- cell_level() %>% filter(Image == img) %>% select(Image, Class, Parent, Centroid_x, Centroid_y)
      # Add cell ID column
      imgData <- imgData %>% mutate(ID = paste0("Cell_", row_number()))
      # Extract spatial info
      coord_x = imgData$Centroid_x
      coord_y = imgData$Centroid_y
      phenotypes = imgData$Class
      dummy_intensity = rep(0, nrow(imgData))
      intensity_matrix = matrix(dummy_intensity, nrow=1, ncol=nrow(imgData))
      colnames(intensity_matrix) = imgData$ID
      
      imgDataList[[img]] = SPIAT::format_image_to_spe(intensity_matrix=intensity_matrix,
                                              phenotypes = phenotypes, coord_x = coord_x,coord_y = coord_y)
    }
    
    return(imgDataList)
    })
  
  output$list_loop <- renderText({
    req(images_as_list())
    paste("Finished processing", length(images_as_list()), "images into a list")
  })
  
  global = reactiveValues(neighborhoods_df = NULL, cell_level_neighborhoods = NULL)
  
  observeEvent(input$run_detection,
               {
                 req(images_as_list()) # Need to put images into lists first
                 
                 # Run neighborhood detection, returns a list with cell level data per image
                 all_neighborhoods <- lapply(images_as_list(), find_communities_list, radius = input$final_radius, neighborhood_size = input$min_neighborhood_size,
                                             cell_types = input$cell_types)
                 
                 print("!!!MESSAGE!!! FINISHED DETECTING NEIGHBORHOODS")
                 
                 # Delete entries for which no neighborhoods were detected
                 all_neighborhoods <- all_neighborhoods %>% purrr::discard(is.null)
                 
                 save(all_neighborhoods, file="cell_level_neighborhood_data.RData")
                 
                 community_summary <- lapply(all_neighborhoods, function(x) {
                   summary_no_location <- x %>% filter(!is.na(Neighborhood)) %>% group_by(Neighborhood, Phenotype) %>% summarise(sum = n()) %>%
                     pivot_wider(names_from = "Phenotype", values_from = "sum")
                   
                   summary_no_location[is.na(summary_no_location)] <- 0
                   
                   summary_no_location <- summary_no_location %>% filter(Neighborhood != "Free_cell")
                 })
                 
                 community_summary <-data.table::rbindlist(community_summary, fill=TRUE, idcol="Image")
                 community_summary[is.na(community_summary)] <- 0
                 print(paste("!!!MESSAGE!!!", "Detected", nrow(community_summary), "neighborhoods in", length(unique(community_summary$Image)), "images"))
                 global$neighborhoods_df <- community_summary
               })
  
  observe({print(head(global$cell_level_neighborhoods))})
  
  tree <- reactive({
    req(global$neighborhoods_df)
    print("Building the dendrogram")
    print(names(global$neighborhoods_df))
    communities_matrix <- global$neighborhoods_df %>% select(-Neighborhood, -Image) %>%
      as.matrix()
    
    categories_of_interest  = input$cell_types
    
    # Convert to proportions
    communities_matrix_prop <- t(apply(communities_matrix,1, function(x) x/sum(x)))
    
    # Calculate dissimilarity and build dendrogram
    hier_bray <- vegan::vegdist(communities_matrix_prop, method = "bray")
    tree <- hclust(hier_bray, method="ward.D2")
    tree_dend <- as.dendrogram(tree)
    return(tree_dend)
  })
  
  output$dendrogram <- renderPlot({
    req(tree())
    print("Plotting the dendrogram")
    clust <- cutree(tree(), k=input$number_of_clusters)
    clust.cutree <- dendextend:::cutree(tree(), k=input$number_of_clusters, order_clusters_as_data = FALSE)
    idx <- order(as.numeric(names(clust.cutree)))
    clust.cutree <- clust.cutree[idx]
    tbl <- table(clust, clust.cutree)
    lbls <- apply(tbl,2,which.max)
    dend1 <- dendextend::color_branches(tree(), k = input$number_of_clusters, groupLabels = lbls)
    dend1 %>% dendextend::set("labels", NULL) %>% plot()
  })
  
  output$bar_plot <- renderPlot({
    req(tree())
    print("Plotting the barplot")
    print(input$colors)
    clust <- cutree(tree(), k=input$number_of_clusters)
    neighborhood_df_with_LMA <- global$neighborhoods_df %>% mutate(New_Neighborhood = clust)

    bar_data <- neighborhood_df_with_LMA %>%
      mutate(COMMUNITY_ID = row_number()) %>%
      select(-c("Image", "Neighborhood")) %>% pivot_longer(cols = -c("COMMUNITY_ID","New_Neighborhood")) %>%
      group_by(COMMUNITY_ID) %>%
      mutate(prop = value/sum(value))
    
      barplot <- bar_data %>%
      ggplot(aes(x=factor(COMMUNITY_ID), y=prop, fill=name))+
      geom_col()+
      facet_wrap(~New_Neighborhood, scales="free")+
      scale_fill_manual(name = "Cell Type",
                        values = (str_split(input$colors, " "))[[1]])+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      xlab(NULL)+
      ylab("# of cells")
      
    print(barplot)
  })
  
  output$summary_stats <- renderTable({
    req(tree())
    print("Calculating summary stats")
    clust <- cutree(tree(), k=input$number_of_clusters)
    neighborhood_df_with_LMA <- global$neighborhoods_df %>% mutate(New_Neighborhood = clust)
    
    bar_data <- neighborhood_df_with_LMA %>%
      select(-c("Image", "Neighborhood")) %>%
      mutate(across(-"New_Neighborhood")/rowSums(across(-"New_Neighborhood"))) %>%
      group_by(New_Neighborhood) %>%
      summarise(across(everything(), mean))
    
  })
  
  output$download <- downloadHandler(
    filename = function() {
      "SPIAT_output.csv"
    },
    content = function(file) {
      clust <- cutree(tree(), k=input$number_of_clusters)
      neighborhood_df_with_LMA <- global$neighborhoods_df %>% mutate(New_Neighborhood = clust)
      vroom::vroom_write(neighborhood_df_with_LMA, file)
    }
  )
}





shinyApp(ui = ui, server = server)
