# DPT_in_R
The implementation of the Discrete Pulse Transform in R is included in this repository. This code implements the DPT through the Roadmaker's Pavage algorithm. 

**Disclaimer:** The code is not optimised in anyway. It may be slow when working with data sets with many regions.
The code will be updated as it is improved.

The code within this repository was mainly written for the implementation on an irregular lattice. 

## 1. RMPA_functions.R 
This file contains all the functions needed to perform the algorithm. Run everything in this file when starting. 

* RMPA_init: This function initialises the RMPA for regular lattice.
* merge_nodes: This function merges nodes that are of equal value and are neighbours.
* create_feature_table: This function creates the feature table containing local minimums and maximums in the working graph
* apply_Un: This function applies the Un operator to the working graph.
* apply_Ln: This function applies the Ln operator to the working graph.
* clear_feature_table: This function deletes the rows currently in the feature table.
* RMPA_implement: This function implements the RMPA and iterates throught the whole process automatically for regular lattice.
* scale_selection: This function extracts certain pulses from the pulse graph for partial reconstruction for regular lattice.
* scale_selection_lattice: This function extracts certain pulses from the pulse graph for irregular lattice.

## 2. RMPA_Example_RegularLattice.R
This file contains an example of the implementation of the RMPA on a regular lattice.

Ensure that the regular lattice is of class "im". This is the image type from the spatstat package. 

```r
library(spatstat)
#change the image type to "im"
as.im(image) 
```

## 3. RMPA_Example_IrregularLattice_sf.R 
This file contains an example of the implementation of the RMPA on an irregular lattice when using the sf package.

Ensure that you are handling the irregular lattice with the sf package when making use of the method in this file.

## 4. RMPA_Example_IrregularLattice_sp.R 
This file contains an example of the implementation of the RMPA on an irregular lattice when using the sp package.

Ensure that you are handling the irregular lattice with the sp package when making use of the method in this file.