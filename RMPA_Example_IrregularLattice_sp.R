# In this script, an example will be shown of the implementation of the RMPA on
# an irregular lattice. 
# We will make use of the sp package

library(sp)
library(raster)
library(maptools)
library(rgdal)
library(spdep)
library(Matrix)

# For illustration purposes, a built-in irregular lattice data will be used
eire <- readOGR(system.file("shapes/eire.shp", package="spData")[1])

# Create the connectivity matrix
eire.nb <- poly2nb(eire)
# convert the neighbourhood to a sparse matrix to be used in the RMPA
d <- length(eire.nb)
neigh_sparse <- Matrix(matrix(rep(0, d*d),nrow = d), sparse = TRUE)

for(i in 1:d){
  d_neigh <- unlist(eire.nb[i])
  for(j in 1:length(d_neigh)){
    if(i != d_neigh[j] & neigh_sparse[i,d_neigh[j]] != 1){
      neigh_sparse[i,d_neigh[j]] <- 1
    }
  }
}

# Initialise the RMPA
RMPA_init(eire, attr = "INCOME", neigh_sparse)

# Implement the RMPA
RMPA_implement()


# Extract certain scales and plot them
extracted <- scale_selection_lattice(eire, 1)
spplot(extracted, "pulses")

extracted2 <- scale_selection_lattice(eire, c(1,4), interval = TRUE)
spplot(extracted2, "pulses")

extracted3 <- scale_selection_lattice(eire, c(1,4), interval = TRUE, binary = TRUE)
spplot(extracted3, "pulses")








