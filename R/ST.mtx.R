
#' Creating Seurat object from mtx
#'
#' Normally the Visium Seurat objects are made directly from HDF5 feature-barcade matrix format.
#' However, somethimes the HDF5 (.h5) format is not available when using public data.
#' Using ST.mtx, one can create the Visium Seurat objetcs directly from the gene
#' count matrix and Spatial transcriptomics images.
#'
#' @param Path Path to the SpaceRnager output directory which contains the filtered_feature_bc_matrix folder and spatial folders
#' @param Name Name of the project
#' @return A Visium Seurat object
#' @export
ST.mtx <- function(path # The Spaceranger output directory which includes filtered_feature_bc_matrix and spatial folders
                   ,name # The project name
){
  x1 <- Read10X(paste0(path,"filtered_feature_bc_matrix/"),gene.column = 1) # Reading in the gene count matrix
  x2 <- Read10X_Image(paste0(path,"spatial/")) # Reading in the ST image
  x2@coordinates <- x2@coordinates[intersect(rownames(x2@coordinates),x1@Dimnames[[2]]),]
  x1 <- x1[,rownames(x2@coordinates)]
  x1 <- CreateSeuratObject(counts = x1,assay = "Spatial",project = name) # Creating the Seurat object
  x1@images$x1 <- Read10X_Image(paste0(path,"spatial/"),filter.matrix = T) # Adding the image to the Seurat object
  x1@images$x1@assay <- "Spatial" # Adding the metadata required for the image
  x1@images$x1@key <- tolower(paste0(name,"_")) # Adding the metadata required for the image
  x1@images$x1@coordinates <- x1@images$x1@coordinates[na.omit(match(names(as.data.frame(x1@assays$Spatial@counts)),
                                                                     row.names(x1@images$x1@coordinates))),] # Making consistent the ST barcodes in the image and gene counts
  return(x1)
}
