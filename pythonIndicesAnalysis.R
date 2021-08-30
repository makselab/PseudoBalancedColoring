require(ggplot2)
require(foreach)
require(fibrationSymmetries)
require(igraph)
require(keras)
require(MASS)
require(tidyr)
require(dplyr)

getQuotientEigenSpectrum <- function(adjacency_matrix, balanced_coloring) {
  indicator_matrix = to_categorical(balanced_coloring)[, -1]
  
  Q = ginv(indicator_matrix) %*% adjacency_matrix %*% indicator_matrix
  setQ = eigen(Q)$values
  
  to_remove = NULL
  eigen_data = eigen(adjacency_matrix)
  for(j in 1:length(setQ)) {
    to_remove = c(to_remove, which.min(abs(eigen_data$values - setQ[j])))
  }
  
  return(eigen_data$values[-to_remove])
}


setwd("/home/ian/Dropbox (City College)/Research/PhD work/shared folders/PROTEIN-FOLDING/COLLABORATION-ANALYSIS/DAVID-FRANCESCO/FINAL-ANALYSIS-MAY-2021/")

directoryName = "INTEGER-PROGRAM-COST-EDGES"
# directoryName = "HUBS-ONLY-NETWORKS"
# directoryName = "INTEGER-PROGRAM-COST-ONE-OVER-MAX-DEGREE"

directories = list.dirs(paste0("DATA/", directoryName, "/OUTPUT/"))
# directories = list.dirs(paste0("DATA/", directoryName))
directories = directories[!grepl("OUTPUT/$", directories)]

dirId = 1
files = list.files(directories[dirId])
files = files[grepl("_indices.txt", files)]

indices = foreach(i = 1:length(files), .combine = rbind) %do% {
  ncolors = as.integer(gsub(".*_([0-9]+)_indices.txt", "\\1", files[i]))
  edgeList = read.table(paste0(directories[dirId], "/", gsub("_indices.txt", "_edges.txt", files[i])), stringsAsFactors = F, sep = ",", header = T)
  balancedColoring = get.balanced.coloring.Kamei(raw_edges = edgeList[, 1:2])
  graph = graph_from_edgelist(as.matrix(edgeList[, 1:2]), directed = F)
  adjacency_matrix = as.matrix(get.adjacency(graph))
  kappa = max(eigen(adjacency_matrix)$value) + 1
  shifted_adjacency_matrix = adjacency_matrix - kappa * diag(ncol(adjacency_matrix))
  quotient_eigenspectrum = getQuotientEigenSpectrum(adjacency_matrix = shifted_adjacency_matrix, balanced_coloring = balancedColoring$Color)
  quotient_eigenspectrum = Re(quotient_eigenspectrum)
  
  index = read.table(paste0(directories[dirId], "/", files[i]), stringsAsFactors = F, sep = ",")
  index[grepl("epsilon", index$V1), 2] = index[grepl("epsilon", index$V1), 2] / nrow(edgeList) / 4
  index = rbind(index, c("Number of Trivial Colors", sum(table(balancedColoring$Color) == 1)))
  index = rbind(index, c("Number of Non-Trivial Colors", max(balancedColoring$Color) - sum(table(balancedColoring$Color) == 1)))
  index = rbind(index, c("Number of Nodes in Non-Trivial Colors", nrow(balancedColoring) - sum(table(balancedColoring$Color) == 1)))
  index = rbind(index, c("Critical Eigenvalue Shifted Adjacency", max(quotient_eigenspectrum)))
  index$V3 = ncolors
  
  return(index)
}
indices$V2 = as.double(indices$V2)
indices$V3 = as.double(indices$V3)

# variables = c("Number of Trivial Colors", "Number of Non-Trivial Colors", "Number of Nodes in Non-Trivial Colors", "Normalized Fiedler Value", "Mean epsilon", "Critical Eigenvalue Shifted Adjacency")
variables = c("Number of Non-Trivial Colors", "Number of Trivial Colors", "Number of Nodes in Non-Trivial Colors", "Normalized Fiedler Value", "Max epsilon")
indices = filter(indices, V1 %in% variables)
indices$V1 = factor(indices$V1, levels = variables)
# indices[indices$V1 == "Max epsilon", 2] = indices[indices$V1 == "Max epsilon", 2] / 53
# indices[indices$V1 == "Max epsilon", 2] = indices[indices$V1 == "Max epsilon", 2] / 59
# # indices = filter(indices, V1 != "Mean epsilon all")
# 
manualIndices = data.frame(V1 = read.table("manualRepairIndices", sep = "\n", stringsAsFactors = F)[c(2:4, 9), ], stringsAsFactors = F)
# # manualIndices = data.frame(V1 = read.table("manualRepairIndices", sep = "\n", stringsAsFactors = F)[c(7), ], stringsAsFactors = F)
# manualIndices = data.frame(V1 = read.table("manualRepairIndices", sep = "\n", stringsAsFactors = F)[c(17:19, 24), ], stringsAsFactors = F)
# # manualIndices = data.frame(V1 = read.table("manualRepairIndices", sep = "\n", stringsAsFactors = F)[c(21), ], stringsAsFactors = F)
manualIndices = separate(manualIndices, V1, c("V1", "V2"), sep = "\t")
manualIndices$V1 = factor(manualIndices$V1, levels = variables)
manualIndices$V2 = as.double(manualIndices$V2)
# # c("Mean epsilon", nrow(edgelist) * 0.25 * 4)

# require(egg)

plotData = 
  ggplot(data = indices) +
    geom_hline(data = manualIndices, aes(yintercept = V2), color = "red", size = 2) +
    geom_vline(xintercept = 9, color = "red", size = 2) +
    geom_vline(xintercept = 9, color = "green", size = 2, linetype = "dashed") +
    # geom_vline(xintercept = 8, color = "red", size = 2) +
    # geom_vline(xintercept = 11, color = "green", size = 2) +
    # geom_vline(xintercept = 12, color = "green", size = 2) +
    # geom_rect(aes(xmin = 11, xmax = 12, ymin = -Inf, ymax = Inf), fill = "palegreen", alpha = 0.11) +
    geom_point(mapping = aes(x = V3, y = V2), size = 8) +
    geom_line(mapping = aes(x = V3, y = V2), size = 2) +
    facet_wrap(~V1, scales = "free", ncol = 1) +
    theme_bw() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 40, face = "bold"),
          strip.text = element_text(size = 35, face = "bold"),
          panel.spacing = unit(1.5, "lines"),
          axis.title.x = element_text(vjust = 0.5, hjust = 1),
          panel.grid.major = element_line(colour = "black")) +
    scale_x_continuous(breaks = 1:17, name = "Number of colors in a repaired graph") +
    scale_y_continuous(name = element_blank())
  # plotData =
#   tag_facet_outside(plotData,
#                     tag_fun_top = function(i) c("a", "b", "c", "d", "e")[i])
ggsave("IndicesSummary.png", plot = plotData, width = 1, height = 2, scale = 10)
# ggsave("IndicesSummary.png", plot = plotData, width = 5, height = 2.5, scale = 8)
