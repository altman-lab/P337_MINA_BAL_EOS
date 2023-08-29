# dat : matrix used to make heatmap
# hm : ComplexHeatmap object
# dimension : Either "row" or "col" to get clusters from

# Example
# cluster.row <- get_hm_cluster(dat = corr.R, hm = hm.up, dimension = "row")

get_hm_cluster <- function(dat, hm, dimension){
  cluster.result <- data.frame()
  
  if(dimension == "row"){
    #get clusters, rows of heatmap
    for (i in 1:length(row_order(hm))){
      #Pull clusters
      cluster.result <- t(t(row.names(dat[row_order(hm)[[as.character(i)]],]))) %>% 
        #Convert to data frame
        as.data.frame() %>% 
        #add cluster name
        mutate(cluster = paste0("cluster", i)) %>% 
        #Rename default column
        rename(row=V1) %>% 
        #concatenate results
        bind_rows(cluster.result)
    }
  } else if(dimension == "col"){
    # get clusters, columns of heatmap
    for (i in 1:length(column_order(hm))){
      cluster.result <- t(t(colnames(dat[,column_order(hm)[[as.character(i)]]]))) %>% 
        as.data.frame() %>% 
        mutate(cluster = paste0("cluster", i)) %>% 
        rename(col=V1) %>% 
        bind_rows(cluster.result)
    }
  } else{ stop("dimension must be one of row or col.") }
  
  return(cluster.result)
}
