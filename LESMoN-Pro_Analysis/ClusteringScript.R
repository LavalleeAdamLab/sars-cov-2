require(graphics)
library("gplots")

setwd("~/Documents/Lavallée-Adam_Lab/SARS-CoV-2/EnumerateMotifs/")

distance_matrix <- read.csv("SignificantMotifs_DistanceMatrix_0.05FDR_0.25identity_0.75accepted.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
m <- data.matrix(distance_matrix)

# convert infinity distance to a value one greater than the max non-infinity distance in the matrix
m[is.infinite(m)] <- max(m[!is.infinite(m)]) + 1


# get clusters
row_clusters <- hclust(dist(m), method="average") # using the average method
row_clusters <- hclust(dist(m), method="ward.D") # using the ward method
row_clusters <- hclust(dist(m), method="complete") # using the ward method

# break into clusters
clusters <- cutree(row_clusters, 9)

#Function that convert cluster tree to a list
cluster_tree_to_list <- function(cluster_tree){
  cluster_names=unique(cluster_tree)

  list_of_motif_clusters=list()

  for (i in 1:length(cluster_names)) {
    list_of_motif_clusters[[i]]<-names(which(clusters==i))

  }

  return(list_of_motif_clusters)
}

#Generate the list (R version of list not the Python version) of clusters with their corresponding motifs
list_of_motifs_in_clusters=cluster_tree_to_list(clusters)

#Name of the file with p_value of significant motifs
p_value_table_name="krogan_significantMotifs_details_FDR0.05_proteinsMapped_filteredHomology_0.25identity_0.75accepted.txt"

#Read the p_value motif table
motif_pvalue=read.table(p_value_table_name,header = TRUE)

#Make a data frame that stores the clusters and their assocauited proteins
storage <- data.frame(matrix(data = NA, nrow = length(list_of_motifs_in_clusters),ncol=2))
colnames(storage) <- c("cluster_number", "motif(s)")

#For each cluster, find the motifs with lowest p-values
#Then add it to the dataframe storage
for (i in 1:length(list_of_motifs_in_clusters)) {
  p_value_info <- motif_pvalue[motif_pvalue$Motif %in% list_of_motifs_in_clusters[[i]],][,c(1,3)]
  min_p_value <- min(p_value_info[,2]) #Get the mimimum p_value
  most_sig_motifs <-  p_value_info[p_value_info$P.value == min_p_value,]$Motif #Get the motif(s) that have the minimum p-value
  storage[i,1] <- i
  storage[i,2] <- paste(most_sig_motifs,collapse = "|")
}

#Write the storage into a table
write.table(storage,file = "clusters_sig_motifs_details_FDR0.05_proteinsMapped_filteredHomology_0.25identity_0.75accepted_9clusters_completelinkage.tsv",sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)


#######################################################################################################
#
#This is where the python program should be executed.
#
#
######################################################################################################

require(ggplot2)
require("ggseqlogo")

#Read the cluster and sequences generated from the python
cluster_motif_data <- read.table("clusters_corresponding_original_motifs_details_FDR0.05_proteinsMapped_filteredHomology_0.25identity_0.75accepted_9clusters_wardlinkage.tsv",header = TRUE)

#Function that converts the dataframe (containing cluster and original sequences) into a list
df_to_list <- function(dataframe){
  temp_list=list()

  for (i in 1:nrow(dataframe)) {
    temp_list[paste(as.character(dataframe[i,1]),as.character(dataframe[i,2]),sep = " ")] <-c(strsplit(as.character(dataframe[i,3]),split = "|",fixed = TRUE))
  }

  return(temp_list)
}

#Get the list of sequence logos
list_of_sequence_clusters<-df_to_list(cluster_motif_data)

#Generate the seqlogo with 4 columns
pdf("ggseqlogo_details_FDR0.05_proteinsMapped_filteredHomology_0.25identity_0.75accepted_9clusters_wardlinkage_bits.pdf")
ggseqlogo(list_of_sequence_clusters, ncol = 3)
dev.off()


# OTHER ANALYSES

# compare the proteins mapping to the chosen representative motif to those mapping to all other proteins in the cluster
cluster_proteins_data <- data.frame(1:nrow(cluster_motif_data))
colnames(cluster_proteins_data) <- "Cluster.id"
number_of_motifs <- numeric(nrow(cluster_motif_data))
proteins_associated_to_representative_motifs <- character(nrow(cluster_motif_data))
proteins_associated_to_clusters <- character(nrow(cluster_motif_data))
for (cluster_id in 1:nrow(cluster_motif_data)) {
  number_of_motifs[cluster_id] <- length(which(clusters == cluster_id))
  representative_motif <- cluster_motif_data$representative_motif[cluster_id]
  proteins_associated_to_representative_motifs[cluster_id] <- motif_pvalue[which(motif_pvalue$Motif == representative_motif), "ProtAccessions"]
  protein_list <- character()
  for (motif in names(clusters[which(clusters == cluster_id)])) {
    proteins_to_add <- strsplit(motif_pvalue[which(motif_pvalue$Motif == motif), "ProtAccessions"], "\\|")[[1]]
    for (protein in proteins_to_add) {
      if (!(protein %in% protein_list)) {
        protein_list <- append(protein_list, protein)
      }
    }
  }
  proteins_associated_to_clusters[cluster_id] <- paste(protein_list, collapse = "|")
}
cluster_proteins_data$Number.of.motifs <- number_of_motifs
cluster_proteins_data$Proteins.associated.to.representative.motif <- proteins_associated_to_representative_motifs
cluster_proteins_data$Proteins.associated.to.clusters <- proteins_associated_to_clusters
cluster_proteins_data$Number.of.proteins.representative.motif <-
  sapply(cluster_proteins_data$Proteins.associated.to.representative.motif,
         function(x) length(strsplit(x, "\\|")[[1]]))
cluster_proteins_data$Number.of.proteins.cluster <-
  sapply(cluster_proteins_data$Proteins.associated.to.clusters,
         function(x) length(strsplit(x, "\\|")[[1]]))


write.csv(cluster_proteins_data, "motif_associated_proteins_details_FDR0.05_proteinsMapped_filteredHomology_0.25identity_0.75accepted_9clusters_completelinkage.csv")








# store cluster information in a data frame
clusters_df <- data.frame(Motif = names(clusters), Cluster=clusters)
rownames(clusters_df) <- 1:nrow(clusters_df)

# add p-value and associated protein data to the motifs
additional_motif_data_df <- read.table("krogan_significantMotifs_details_FDR0.05_proteinsMapped_filteredHomology_0.25identity_0.75accepted.txt",
                                     header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep="\t")

# sort by motif and cbind
clusters_df <- clusters_df[order(clusters_df$Motif),]
additional_motif_data_df <- additional_motif_data_df[order(rownames(additional_motif_data_df)),]
full_motif_data <- cbind(clusters_df, additional_motif_data_df)

# sort by pvalue
full_motif_data <- full_motif_data[order(full_motif_data$P.value),]

# add FDRs
p_value_to_fdr <- read.csv("covFDRsMono.tsv", stringsAsFactors = FALSE, sep="\t")[order(p_value_to_fdr$threshold),]
fdrs <- sapply(full_motif_data$P.value, function(x) p_value_to_fdr[which(p_value_to_fdr$threshold == ceiling(x*10^6)/10^6), "fdr"])
full_motif_data$FDR <- fdrs

# sort columns
rownames(full_motif_data) <- 1:nrow(full_motif_data)
full_motif_data <- full_motif_data[, c("Motif", "Cluster", "TPD", "P.value", "FDR", "NumOfProteins", "ProtAccessions", "ProteinNames")]

# output the full motif information data frame
write.csv(full_motif_data, "lesmon_pro_motif_summary_table_FDR0.05.tsv", row.names = FALSE)



# heatmaps
my_ward <- function(d) hclust(d, method="ward.D")
my_average <- function(d) hclust(d, method="average")
my_complete <- function(d) hclust(d, method="complete")


heatmap(m,symm=TRUE,labCol=NA,col=topo.colors(100),scale="none",hclustfun=my_average)

heatmap(m,symm=TRUE,labRow=NA,labCol=NA,col=topo.colors(100),scale="none",hclustfun=my_average)

heatmap(m,symm=TRUE,col=topo.colors(100),scale="none",hclustfun=my_average)	# with structure index labels




# average linkage with scale
pdf("heatmapOfClustering_Average.pdf")
heatmap.2(m, symm=TRUE, labRow=NA, labCol=NA, col=topo.colors(100), scale="none", hclustfun=my_average, cexRow=0.25, cexCol=0.25, key=TRUE, density.info="none", trace="none", Colv=TRUE)
dev.off()

# with names
heatmap.2(m, symm=TRUE, col=topo.colors(100), scale="none", hclustfun=my_average, cexRow=0.25, cexCol=0.25, key=TRUE, density.info="none", trace="none", Colv=TRUE)


# ward linkage with scale
pdf("heatmapOfClustering_Ward.pdf")
heatmap.2(m, symm=TRUE, labRow=NA, labCol=NA, col=topo.colors(100), scale="none", hclustfun=my_ward, cexRow=0.25, cexCol=0.25, key=TRUE, density.info="none", trace="none", Colv=TRUE)
dev.off()

# complete linkage
heatmap.2(m, symm=TRUE, labRow=NA, labCol=NA, col=topo.colors(100), scale="none", hclustfun=my_complete, cexRow=0.25, cexCol=0.25, key=TRUE, density.info="none", trace="none", Colv=TRUE)
