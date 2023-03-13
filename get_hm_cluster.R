# dat : matrix used to make heatmap
# hm : ComplexHeatmap object
# dimension : Either "row" or "col" to get clusters from

# Example
# cluster.row <- get_hm_cluster(dat = corr.R, hm = hm.up, dimension = "row")

get_hm_cluster <- function(dat, hm, dimension){
  cluster.result <- data.frame()
  
  # Rows of heatmap
  if(dimension == "row"){
    #Deal with single cluster results
    if(is.list(row_order(hm))){
      clust.tot <- length(row_order(hm))
    } else {
      clust.tot <- 1
    }
    
    for (i in 1:clust.tot){
      #Get row indices
      if(is.list(row_order(hm))){
        cluster.index <- row_order(hm)[[as.character(i)]]
      } else {
        cluster.index <- row_order(hm)
      }
      
      #Clusters with >1 element
      if(length(cluster.index) > 1){
        #Pull clusters
        cluster.result <- t(t(row.names(dat[cluster.index,]))) %>% 
          #Convert to data frame
          as.data.frame() %>% 
          #row order within cluster
          mutate(row_within_cluster = 1:length(cluster.index)) %>% 
          #add cluster name
          mutate(cluster = paste0("cluster", i)) %>% 
          #Rename default column
          rename(row=V1) %>% 
          #concatenate results
          bind_rows(cluster.result)
      } else {
        #Clusters with only 1 element
        cluster.result <- data.frame(row = rownames(dat)[cluster.index],
                                     row_within_cluster = 1,
                                     cluster = paste0("cluster", i)) %>% 
          bind_rows(cluster.result)
      }
    }
  } else if(dimension == "col"){
    # Columns of heatmap
    
    #Deal with single cluster results
    if(is.list(column_order(hm))){
      clust.tot <- length(column_order(hm))
    } else {
      clust.tot <- 1
    }
    
    for (i in 1:clust.tot){
      #Get column indices
      if(is.list(column_order(hm))){
        cluster.index <- column_order(hm)[[as.character(i)]]
      } else {
        cluster.index <- column_order(hm)
      }
      
      #Clusters with >1 element
      if(length(cluster.index) > 1){
        cluster.result <- t(t(colnames(dat[,cluster.index]))) %>% 
          as.data.frame() %>% 
          mutate(col_within_cluster = 1:length(cluster.index)) %>% 
          mutate(cluster = paste0("cluster", i)) %>% 
          rename(col=V1) %>% 
          bind_rows(cluster.result)
      } else {
        #Clusters with only 1 element
        cluster.result <- data.frame(col = colnames(dat)[cluster.index],
                                     col_within_cluster = 1,
                                     cluster = paste0("cluster", i)) %>% 
          bind_rows(cluster.result)
      }
    }
  } else{ stop("dimension must be one of row or col.") }
  
  return(cluster.result)
}
