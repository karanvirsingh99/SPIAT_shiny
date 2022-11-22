# Shiny app with SPIAT demo

The [SPIAT package](https://trigosteam.github.io/SPIAT/articles/introduction.html) is an R package with tools to analyze multiplex immunofluorescence data (mcIF). One of the tools embedded in this package is neighborhood detection, which allows users to identify cell neighborhoods in mcIF neighborhoods. The two main parameters to the SPIAT::identify_neighborhoods functions are the minimum size of the neighborhood, and the maximum distance to consider two cells as being part of the same neighborhood.

This shiny app provides an interactive visualization of how changing these two parameters affects neighborhood detection, using 4 sample mcIF images.

## Prerequisites

This app requires tidyverse, ggsci and SPIAT to be installed. You can run the following code block in your R session to install the required packages.

```{r}
install_packages("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")}
# The following initialises usage of Bioc devel
BiocManager::install(version="devel")
BiocManager::install("SPIAT")

install_packages("ggsci")
```

## How to run

To run this shiny app, run the following command in R:

```{r}
shiny::runGitHub("SPIAT_shiny", "karanvirsingh99")
```
