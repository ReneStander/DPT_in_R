# In this script, an example will be shown of the implementation of the RMPA on
# an irregular lattice. 
# We will make use of the sf package

# Import the packages
library(sf)
library(dplyr)
library(ggplot2)
library(leaflet)
library(scales)
library(ggmap)
library(stringr)
library(Matrix)

# For illustration purposes, a built-in irregular lattice data will be used
nc <- st_read(system.file("shape/nc.shp", package="sf"))

# select only the column to be used (this is not completely necessary)
nc.use <- nc %>% select(BIR74)

# Create the connectivity matrix
(neigh <- st_intersects(nc.use, nc.use))

# convert the neighbourhood to a sparse matrix to be used in the RMPA
d <- length(neigh)
neigh_sparse <- Matrix(matrix(rep(0, d*d),nrow = d), sparse = TRUE)

for(i in 1:d){
  d_neigh <- unlist(neigh[i])
  for(j in 1:length(d_neigh)){
    if(i != d_neigh[j] & neigh_sparse[i,d_neigh[j]] != 1){
      neigh_sparse[i,d_neigh[j]] <- 1
    }
  }
}

# Initialise the RMPA
RMPA_init(data = nc.use, attr = "BIR74", neighbours = neigh_sparse)

# Implement the RMPA
RMPA_implement()


# Extract certain scales and plot them
extracted <- scale_selection_lattice(nc.use, 1)
plot(extracted)

extracted2 <- scale_selection_lattice(nc.use, c(1,4), interval = TRUE)
plot(extracted2)

extracted3 <- scale_selection_lattice(nc.use, c(1,4), interval = TRUE, binary = TRUE)
plot(extracted3)
