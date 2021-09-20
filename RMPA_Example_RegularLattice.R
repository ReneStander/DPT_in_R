# In this script, an example will be shown of the implementation of the RMPA on
# a regular lattice. 

# Import packages
library(tidyverse)
library(spatstat)
library(Matrix) # sparse matrices
library(reticulate) # to reshape arrays
library(igraph) # for the graph plots
library(ggraph) # creating graph plots

# Import the data/image
# for illustration purposes, we will use a small image
org_img <- as.im(matrix(c(2,1,3,5,4,3), ncol =3, byrow = T))
plot(org_img) # this is just to show the image

# Neighbourhood structure
# Disclaimer: This is not the most efficient way to get the neighbourhood structure
neigh <- matrix(c(0,1,0,1,0,0,1,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,1,0), ncol =6, 
                byrow = T)
neigh_sparse <- Matrix(neigh, sparse = TRUE) # This line converts the neighbourhood matrix to a sparse matrix
neigh_sparse # This line of code prints the neighbourhood matrix

# Initialise the RMPA
RMPA_init(org_img)

# Implement the RMPA
RMPA_implement()

#Extract certain scales and plot them

extracted <- scale_selection(org_img, 1, interval = FALSE, binary = FALSE)
plot(extracted)

extracted2 <- scale_selection(org_img, c(1,2), interval = TRUE, binary = FALSE)
plot(extracted2)

extracted3 <- scale_selection(org_img, 1, interval = FALSE, binary = TRUE)
plot(extracted3)
