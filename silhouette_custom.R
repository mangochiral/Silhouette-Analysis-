library(amap)
library(factoextra)
library(readr)
library(tidyverse)
library(cluster)
#Silhoutte Analysis for pearson, spearman distance matrix not available in other 
# package

# For calculation of Kmeans the df can be a range of selected sample
silhoutteCustom <- function(df,k,methodtype){
  
  # Make correlation matrix based on peaks
  diss_cor <- cor(t(df), method = as.character(methodtype))
  
  # Make dissimarility matrix
  spearman_dis <- 1-diss_cor
  
  # Make a empty dataframe to record cluster number for each k means clustering
  AllClusters <- data.frame(matrix(NA, nrow = length(rownames(df))))
  
  #Loop for various Kmeans
  for(i in 2:k){
    #Random seed for reproducibility
    set.seed(123)
    
    #Calculate kmeans for each k
    KmeansAMAP <- Kmeans(as.matrix(df), centers = i, 
                         iter.max = 1000, nstart = 25, method = "spearman")
    ClusterNos <- as.data.frame(KmeansAMAP$cluster)
    colnames(ClusterNos) <- paste0("kmeans_", i)
    
    #Append each kmeans clustering output
    AllClusters <- cbind(AllClusters, ClusterNos)
  }
  # Remove the first NA column
  AllClusters <- AllClusters[, -c(1)]
  
  #Make a empty dataframe to record avg. silhouette score
  silscore <- data.frame(matrix(NA, nrow = length(colnames(AllClusters))))
  colnames(silscore) <- "silscores"
  
  #Loop through each kmeans data by parts for memory issues of computing
  for (t in 1:length(colnames(AllClusters))){
    
    # Split into smaller chunks if needed
    chunk_size <- 1000
    n <- nrow(df)
    
    # Empty vector
    sil_scores <- vector()
    
    #Loop through each chunk
    for(i in seq(1, n, chunk_size)) {
      end <- min(i + chunk_size - 1, n)
      chunk_dist <- as.dist(spearman_dis[i:end, i:end])
      chunk_clusters <- AllClusters[[t]][i:end]
      chunk_sil <- silhouette(chunk_clusters, chunk_dist)
      sil_scores <- c(sil_scores, chunk_sil[, "sil_width"])
    }
    # Calculate overall silhouette score
    mean_sil <- mean(sil_scores)
    silscore[["silscores"]][t] <- mean_sil
  }
  silscore <- silscore |> rownames_to_column("kmeans_val")
  silscore$kmeans_val <- 1+ as.numeric(silscore$kmeans_val)
  # Return silhouette score
  return(as.data.frame(silscore))
  
}

plotSilhoutte <- function(df){
  ggplot(df, 
         aes(x=kmeans_val, y=silscores)) +
    geom_point()+ geom_line()+
    geom_vline(xintercept = 5, linetype="dashed", color="red") +
    geom_vline(xintercept = 2, linetype="dashed", color="blue") +
    scale_x_continuous(breaks = 1:12) +
    theme_minimal() +
    labs(x="K-means Value", y="Silhouette Score")
}

silhOptmial <- silhoutteCustom(dfNormSelected, 12, "spearman")
# write_tsv(silhOptmial, paste0(path, "silscore.tsv"))
plotSilhoutte(silhOptmial)
# ggsave(paste0(path, "Silhouette_Analysis_for_kmeans.png"), width = 8, height = 6)










