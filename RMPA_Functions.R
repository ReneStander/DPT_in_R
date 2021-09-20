# In this script we will define the functions to be used for the RMPA
# Disclaimer: The code shown here is not optimised in any way yet. It will most likely be slow

# This function initialises the RMPA
# works on regular lattices
RMPA_init <- function(im){
  org_img_array <- as.array(im)
  org_img_flatten <<- 
    array_reshape(org_img_array, dim = cbind(1,dim(org_img_array)[1]*dim(org_img_array)[2]))
  
  WG_nodes <<- tibble(id = 1:length(org_img_flatten), Value = t(org_img_flatten), Scale = 1)
  
  edgelist <- get.edgelist(graph.adjacency(neigh_sparse))
  WG_edges <<- tibble(from = edgelist[,1], to = edgelist[,2])
  
  PG_nodes <<- tibble(id = WG_nodes$id, scale = WG_nodes$id, value = 0)
  
  PG_edges <<- tibble(from = numeric(0), to = numeric(0))
  
  WG_PG_edge <<- tibble(from = WG_nodes$id, to = PG_nodes$scale)
  
  feature_table <<- tibble(id = numeric(0), size = numeric(0), value = numeric(0), type = character(0))
}

# This function merges the nodes that are of equal value and are neighbours
# Will work on any type of data
merge_nodes_new <- function(){
  # find nodes of equivalent value in WG
  dups <- duplicated(WG_nodes$Value)
  dups_vals <- unique(WG_nodes$Value[dups])
  
  m <- 1
  
  while(m <= length(dups_vals)){
    (dup_ind <- which(WG_nodes$Value == dups_vals[m]))
    dim(neigh_sparse)
    
    # see if they are neighbours
    if(!is.null(dim(neigh_sparse))){
        for(j in 1:length(dup_ind)){
        for(k in 1:length(dup_ind)){
          dup_neigh <- (neigh_sparse[dup_ind[j], dup_ind[k]] == 1)
          if(dup_neigh == 1) break
        }
        if(dup_neigh == 1) break
      }
    } else{
      dup_neigh <- FALSE
    }
    
    
    dup_neigh
    dup_ind <- dup_ind[c(j,k)]
    
    # if they are neighbours - merge the nodes
    # merge to the smaller index
    if(dup_neigh == TRUE){
      
      # change edges between the pulse and working graph too
      WG_PG_edge[(WG_PG_edge$from == as.numeric(WG_nodes[dup_ind[2], "id"])), "from"] <<-  as.numeric(WG_nodes[dup_ind[1], "id"])
      
      # merge nodes in WG
      WG_nodes[dup_ind[1],"Scale"] <<- WG_nodes[dup_ind[1], "Scale"] + WG_nodes[dup_ind[2], "Scale"]
      # remove row in WG
      WG_nodes <<- WG_nodes[-dup_ind[2],]
      
      # update neighbourhood matrix
      neigh_sparse[dup_ind[1],] <<- neigh_sparse[dup_ind[1],] + neigh_sparse[dup_ind[2],]
      neigh_sparse[,dup_ind[1]] <<- neigh_sparse[,dup_ind[1]] + neigh_sparse[,dup_ind[2]]
      neigh_sparse <<- neigh_sparse[-dup_ind[2],]
      if(!is.null(ncol(neigh_sparse))){
        neigh_sparse <<- neigh_sparse[,-dup_ind[2]]
        neigh_sparse[dup_ind[1], dup_ind[1]] <<- 0
      }else{
        neigh_sparse <<- neigh_sparse[1]
        neigh_sparse <<- 0
      }
      
      neigh_sparse[which(neigh_sparse>1)] <<- 1
      
      # update edge list
      edgelist <- get.edgelist(graph.adjacency(neigh_sparse, mode = 'undirected'))
      WG_edges <<- tibble(from = edgelist[,1], to = edgelist[,2])
    } else{
      m = m + 1
      
    }
  }
}

# This function creates the feature table
# Will work on any data type
create_feature_table <- function(){
  for(i in 1:nrow(WG_nodes)){
    diff <- WG_nodes$Value[i] - WG_nodes$Value
    
    d_n <- diff[neigh_sparse[i,]==1]
    if(sum(d_n>0) == length(d_n)){
      feature_table <<- feature_table %>% add_row(id = WG_nodes$id[i], 
                                                 size = WG_nodes$Scale[i],
                                                 value = WG_nodes$Value[i],
                                                 type = "b")
    }
    
    if(sum(d_n<0) == length(d_n)){
      feature_table <<- feature_table %>% add_row(id = WG_nodes$id[i], 
                                                 size = WG_nodes$Scale[i],
                                                 value = WG_nodes$Value[i],
                                                 type = "p")
    }
  }
}

#This function applies the Un operator on the data
# Will work on any data type
apply_Un <- function(n){
  to_smooth <- feature_table %>% filter(type == "p" & size == n)
  for(i in 1:nrow(to_smooth)){
    d <- as.numeric(to_smooth$id[i])
    d1 <- which(d == WG_nodes$id)[1]
    # find neighbours
    nn <- neigh_sparse[d1,]
    # find values at neighbours
    vals <- WG_nodes$Value[(nn == 1)]
    # value to smooth to
    new_val <- min(vals)
    
    # add node to pulse graph
    node_id <- nrow(PG_nodes) + 1
    
    PG_nodes <<- PG_nodes %>% add_row(id = node_id, 
                                     scale = n,
                                     value = (WG_nodes$Value[d1] - new_val))
    
    d_edges <- WG_PG_edge %>% filter(from == d)
    for(j in 1:nrow(d_edges)){
      PG_edges <<- PG_edges %>% add_row(from = d_edges$to[j], to = node_id)
    }
    
    # replace value
    WG_nodes$Value[d1] <<- new_val
  }
}

# This function applies the Ln operator to the data
# Will work on any data type
apply_Ln <- function(n){
  to_smooth <- feature_table %>% filter(type == "b" & size == n)
  for(i in 1:nrow(to_smooth)){
    d <- as.numeric(to_smooth$id[i])
    d1 <- which(d == WG_nodes$id)[1]
    # find neighbours
    nn <- neigh_sparse[d1,]
    # find values at neighbours
    vals <- WG_nodes$Value[(nn == 1)]
    # value to smooth to
    new_val <- max(vals)
    
    # add node to pulse graph
    node_id <- nrow(PG_nodes) + 1
    
    PG_nodes <<- PG_nodes %>% add_row(id = node_id,
                                     scale = n,
                                     value = (WG_nodes$Value[d1] - new_val))
    
    d_edges <- WG_PG_edge %>% filter(from == d)
    for(j in 1:nrow(d_edges)){
      PG_edges <<- PG_edges %>% add_row(from = d_edges$to[j], to = node_id)
    }
    
    # replace value
    WG_nodes$Value[d1] <<- new_val
  }
}

# This function deletes all the rows in the feature table
# Will work on any data type
clear_feat_table <- function(){
  feature_table <<- feature_table[-seq(1,nrow(feature_table),1),]
}

# This function implements the RMPA and iterates through the whole process automatically
# Will only work on regular lattice
RMPA_implement <- function(){
  # Step 1: Merge nodes
  merge_nodes_new()
  
  # Step 2: Apply LULU operators from n = 1 up to n = N
  N <- length(org_img_flatten)
  for(n in 1:N){
    nr <- nrow(WG_nodes)
    if(nr > 1){
      create_feature_table()
      
      if(sum(feature_table$size == n & feature_table$type == "b") > 0)  apply_Ln(n)
      if(nrow(feature_table) > 0) clear_feat_table()
      create_feature_table()
      
      if(sum(feature_table$size == n & feature_table$type == 'p') > 0)  apply_Un(n)
      merge_nodes_new()
      if(nrow(feature_table) > 0) clear_feat_table()
      
    }
  }
}

# This function extracts certain pulses from the pulse graph
# Will only work on grid images
scale_selection <- function(im, scales, interval = FALSE, binary = FALSE){
  
  org_img_flatten_copy <- Matrix(org_img_flatten, sparse = TRUE)
  org_img_flatten_copy[which(org_img_flatten_copy != 0)] <- 0
  
  if(interval == FALSE){
    PG_nodes_extracted <- PG_nodes %>% filter(scale %in% scales & value !=0)
    
    for(i in 1:nrow(PG_nodes_extracted)){
      d <- as.numeric(PG_nodes_extracted$id[i])
      if(binary == FALSE){
        org_img_flatten_copy[PG_edges$from[which(PG_edges$to == d)]] <- PG_nodes_extracted$value[i]
      } else{
        org_img_flatten_copy[PG_edges$from[which(PG_edges$to == d)]] <- 1
      }
    }
    
  } else{
    PG_nodes_extracted <- PG_nodes %>% filter(scale >= scales[1] & scale <= scales[2]  & value !=0)
    d_scales <- seq(scales[1],scales[2],1)
    for(i in 1:length(d_scales)){
      
      pulse <- PG_nodes_extracted %>% filter(scale == d_scales[i])
      
      for(j in 1:nrow(pulse)){
        d <- as.numeric(pulse$id[j])
        dd <- PG_edges$from[which(PG_edges$to == d)]
        
        if(binary == FALSE){
          org_img_flatten_copy[dd] <- org_img_flatten_copy[dd] + pulse$value[j]
        } else{
          if(length(dd) > 1){
            for(k in 1:length(dd)){
              if(org_img_flatten_copy[dd[k]] != 1){
              org_img_flatten_copy[dd[k]] <- 1
              }
            }
          } else{
            if(org_img_flatten_copy[dd] != 1){
              org_img_flatten_copy[dd] <- 1
            }
          }
          
          }
        }
        
      }
      
    }
  
  return(as.im(array_reshape(as.array(org_img_flatten_copy), dim = dim(im))))
}

# This function initialises the RMPA when you have lattice data and are using the 
# sf package
RMPA_init_lattice <- function(data, attr = NULL, neighbours){
  WG_nodes <<- tibble(id = 1:nrow(data), Value = data[,attr], Scale = 1)
  
  edgelist <- get.edgelist(graph.adjacency(neighbours))
  WG_edges <<- tibble(from = edgelist[,1], to = edgelist[,2])
  
  PG_nodes <<- tibble(id = WG_nodes$id, scale = WG_nodes$id, value = 0)
  
  PG_edges <<- tibble(from = numeric(0), to = numeric(0))
  
  WG_PG_edge <<- tibble(from = WG_nodes$id, to = PG_nodes$scale)
  
  feature_table <<- tibble(id = numeric(0), size = numeric(0), value = numeric(0), type = character(0))
}

# This function implements the RMPA and iterates through the whole process automatically
# Will only work on irregular lattice
RMPA_implement_lattice <- function(){
  # Step 1: Merge nodes
  merge_nodes_new()
  
  # Step 2: Apply LULU operators from n = 1 up to n = N
  N <- 27
  for(n in 1:N){
    nr <- nrow(WG_nodes)
    if(nr > 1){
      create_feature_table()
      
      if(sum(feature_table$size == n & feature_table$type == "b") > 0)  apply_Ln(n)
      if(nrow(feature_table) > 0) clear_feat_table()
      create_feature_table()
      
      if(sum(feature_table$size == n & feature_table$type == 'p') > 0)  apply_Un(n)
      merge_nodes_new()
      if(nrow(feature_table) > 0) clear_feat_table()
      
    }
  }
}

# This function initialises the RMPA when you have lattice data from the sp package
RMPA_init_lattice_sp <- function(data, attr = NULL, neighbours){
  
  WG_nodes <<- tibble(id = 1:nrow(data), Value = data[[attr]], Scale = 1)
  
  edgelist <- get.edgelist(graph.adjacency(neighbours))
  WG_edges <<- tibble(from = edgelist[,1], to = edgelist[,2])
  
  PG_nodes <<- tibble(id = WG_nodes$id, scale = WG_nodes$id, value = 0)
  
  PG_edges <<- tibble(from = numeric(0), to = numeric(0))
  
  WG_PG_edge <<- tibble(from = WG_nodes$id, to = PG_nodes$scale)
  
  feature_table <<- tibble(id = numeric(0), size = numeric(0), value = numeric(0), type = character(0))
}

# This function extracts certain pulses from the pulse graph
# Will only work on lattice data
scale_selection_lattice <- function(lat_dat, scales, interval = FALSE, binary = FALSE){
  
  lat_dat$pulses <- rep(0, nrow(lat_dat))
  
  if(interval == FALSE){
    PG_nodes_extracted <- PG_nodes %>% filter(scale %in% scales & value !=0)
    
    for(i in 1:nrow(PG_nodes_extracted)){
      d <- as.numeric(PG_nodes_extracted$id[i])
      if(binary == FALSE){
        lat_dat$pulses[PG_edges$from[which(PG_edges$to == d)]] <- PG_nodes_extracted$value[i]
      } else{
        lat_dat$pulses[PG_edges$from[which(PG_edges$to == d)]] <- 1
      }
    }
    
  } else{
    PG_nodes_extracted <- PG_nodes %>% filter(scale >= scales[1] & scale <= scales[2]  & value !=0)
    d_scales <- seq(scales[1],scales[2],1)
    for(i in 1:length(d_scales)){
      
      pulse <- PG_nodes_extracted %>% filter(scale == d_scales[i])
      
      for(j in 1:nrow(pulse)){
        d <- as.numeric(pulse$id[j])
        dd <- PG_edges$from[which(PG_edges$to == d)]
        
        if(binary == FALSE){
          lat_dat$pulses[dd] <- lat_dat$pulses[dd] + pulse$value[j]
        } else{
          if(length(dd) > 1){
            for(k in 1:length(dd)){
              if(lat_dat$pulses[dd[k]] != 1){
                lat_dat$pulses[dd[k]] <- 1
              }
            }
          } else{
            if(lat_dat$pulses[dd] != 1){
              lat_dat$pulses[dd] <- 1
            }
          }
          
        }
      }
      
    }
    
  }
  
  return(lat_dat)
}

