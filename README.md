# ST.mtx
R package for creating Visium Seurat object without HDF5 file

## Install:
devtools::install_github("https://github.com/ElyasMo/ST.mtx")

## Vignette:

```{r}
library(ST.mtx)

Directory <- "Path" # The directory to the Spaceranger output which includes the
  #"filtered_feature_bc_matrix/" and "Spatial" folders.
Name <- "Project name" # A character which explaines the project name. It will be 
  #place in the orig.ident column of your metadata and also the image key. 

Seurat_object <- ST.mtx(Directory,Name)

```

## Cite:
... (Inpress)

