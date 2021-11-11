require(ggplot2)
require(foreach)
require(fibrationSymmetries)
require(igraph)
require(keras)
require(MASS)
require(tidyr)
require(dplyr)

getIndexFiles <- function(directoryName, prefix) {
  directories = list.dirs(paste0(directoryName, "/OUTPUT/"))
  #directories = directories[!grepl("OUTPUT/$", directories)]
  directories = directories[grepl(prefix, directories)]
  
  files = list.files(directories, full.names = T)
  files = files[grepl("_indices.txt", files)]
  
  return(files)
}

makePlot <- function(nodes, edges, coordinates, pdfname) {
  nodes = merge(nodes, coordinates)
  
  edges$ColorCode = gsub("Original", "#000000", edges$Color)
  edges$ColorCode = gsub("Repaired", "#FF0000", edges$ColorCode)
  
  graph = graph_from_edgelist(as.matrix(edges[, 1:2]), directed = F)
  nodes$Id = factor(nodes$Id, levels = as.character(V(graph)$name))
  nodes = arrange(nodes, Id)
  
  pdf(pdfname, width = 10, height = 10)
  plot(graph,
       layout = as.matrix(nodes[, c("x", "y")]),
       vertex.size = 15,
       vertex.color = nodes$IPColor,
       edge.width = 2,
       edge.color = edges$ColorCode,
       vertex.label.cex = 1,
       vertex.label.color = "#000000"
  )
  dev.off()
}
  
getIndices <- function(files) {
  indices = foreach(i = 1:length(files), .combine = rbind) %do% {
    ncolors = as.integer(gsub(".*_([0-9]+)_indices.txt", "\\1", files[i]))
    edgeList = read.table(paste0(gsub("_indices.txt", "_edges.csv", files[i])), stringsAsFactors = F, sep = ",", header = T)
    balancedColoring = get.balanced.coloring.Kamei(raw_edges = edgeList[, 1:2])
    graph = graph_from_edgelist(as.matrix(edgeList[, 1:2]), directed = F)
    adjacency_matrix = as.matrix(get.adjacency(graph))
    kappa = max(eigen(adjacency_matrix)$value) + 1
    shifted_adjacency_matrix = adjacency_matrix - kappa * diag(ncol(adjacency_matrix))
    quotient_eigenspectrum = getQuotientEigenSpectrum(adjacency_matrix = shifted_adjacency_matrix, balanced_coloring = balancedColoring$Color)
    quotient_eigenspectrum = Re(quotient_eigenspectrum)
    
    index = read.table(paste0(files[i]), stringsAsFactors = F, sep = ",")
    index[grepl("epsilon", index$V1), 2] = index[grepl("epsilon", index$V1), 2] / nrow(edgeList) / 4
    index = rbind(index, c("Number of Trivial Colors", sum(table(balancedColoring$Color) == 1)))
    index = rbind(index, c("Number of Non-Trivial Colors", max(balancedColoring$Color) - sum(table(balancedColoring$Color) == 1)))
    index = rbind(index, c("Number of Nodes in Non-Trivial Colors", nrow(balancedColoring) - sum(table(balancedColoring$Color) == 1)))
    index = rbind(index, c("Critical Eigenvalue Shifted Adjacency", max(quotient_eigenspectrum)))
    index$V3 = ncolors
    index$V4 = length(unique(balancedColoring$Color))
    
    # add colors to node file
    nodes = read.table(paste0(gsub("_indices.txt", "_nodes.csv", files[i])), stringsAsFactors = F, sep = ",", header = T)
    colnames(balancedColoring)[1] = "Id"
    nodes = merge(nodes, balancedColoring)
    nodes = arrange(nodes, Sector, Orbit)
    nodes = nodes[, c("Id", "Label", "Sector", "Orbit", "Color", "IPColor", "Fixed")]
    write.table(nodes, paste0(gsub("_indices.txt", "_nodes.csv", files[i])), sep = ",", row.names = F, quote = F)
    
    makePlot(nodes, edgeList, coordinates, gsub("_indices.txt", ".pdf", files[i]))
    
    return(index)
  }
  indices$V2 = as.double(indices$V2)
  indices$V3 = as.double(indices$V3)
  
  return(indices)
}

getManualIndices <- function(file, prefix) {
  manualIndices = read.table(file, sep = "\n", stringsAsFactors = F)
  manualIndices = separate(manualIndices, V1, c("V1", "V2"), sep = "\t")
  if(prefix == "bw") {
    manualIndices = manualIndices[2:15,]
  } else {
    manualIndices = manualIndices[17:30,]
  }
  manualIndices$V2 = as.double(manualIndices$V2)
  return(manualIndices)
}

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

filterIndices <- function(indices, variables) {
  indices = filter(indices, V1 %in% variables)
  indices$V1 = factor(indices$V1, levels = variables)
  return(indices)
}

directoryName = "/home/ian/Dropbox (City College)/Research/PhD work/shared folders/PROTEIN-FOLDING/COLLABORATION-ANALYSIS/DAVID-FRANCESCO/FINAL-ANALYSIS-MAY-2021/DATA/INTEGER-PROGRAM-COST-EDGES"
# directoryName = "HUBS-ONLY-NETWORKS"
# directoryName = "/home/ian/Dropbox (City College)/Research/PhD work/shared folders/PROTEIN-FOLDING/COLLABORATION-ANALYSIS/DAVID-FRANCESCO/FINAL-ANALYSIS-MAY-2021/DATA/INTEGER-PROGRAM-COST-ONE-OVER-MAX-DEGREE"
# directoryName = "/home/ian/Dropbox (City College)/Research/PhD work/shared folders/PROTEIN-FOLDING/COLLABORATION-ANALYSIS/DAVID-FRANCESCO/FINAL-ANALYSIS-MAY-2021/DATA/INTEGER-PROGRAM-COST-ONE-OVER-SUM-DEGREES"
# directoryName = "/home/ian/Dropbox (City College)/Research/PhD work/shared folders/PROTEIN-FOLDING/COLLABORATION-ANALYSIS/DAVID-FRANCESCO/FINAL-ANALYSIS-MAY-2021/DATA/INTEGER-PROGRAM-COST-ONE-OVER-MULTI-DEGREE"
prefix = "c_fw_no"
manualIndexFile = "/home/ian/Dropbox (City College)/Research/PhD work/shared folders/PROTEIN-FOLDING/COLLABORATION-ANALYSIS/DAVID-FRANCESCO/FINAL-ANALYSIS-MAY-2021/CODE/manualRepairIndices"
if(prefix == "c_fw_no") {
  coordinates = read.table("/home/ian/Dropbox (City College)/Research/PhD work/shared folders/PROTEIN-FOLDING/COLLABORATION-ANALYSIS/DAVID-FRANCESCO/FINAL-ANALYSIS-MAY-2021/DATA/Forward_gap_coordinates.txt", header = T, stringsAsFactors = F)
} else {
  coordinates = read.table("/home/ian/Dropbox (City College)/Research/PhD work/shared folders/PROTEIN-FOLDING/COLLABORATION-ANALYSIS/DAVID-FRANCESCO/FINAL-ANALYSIS-MAY-2021/DATA/Backward_gap_coordinates.txt", header = T, stringsAsFactors = F)
}
colnames(coordinates)[1] = "Id"

# variables = c("Number of Trivial Colors", "Number of Non-Trivial Colors", "Number of Nodes in Non-Trivial Colors", "Normalized Fiedler Value", "Mean epsilon", "Critical Eigenvalue Shifted Adjacency")
variables = c("Number of Non-Trivial Colors", "Number of Trivial Colors", "Number of Nodes in Non-Trivial Colors", "Normalized Fiedler Value", "Max epsilon")

files = getIndexFiles(directoryName, prefix)
indices = getIndices(files)
manualIndices = getManualIndices(manualIndexFile, prefix)
indices = filterIndices(indices, variables)
manualIndices = filterIndices(manualIndices, variables)

# plotData =
#   ggplot(data = indices) +
#     geom_hline(data = manualIndices, aes(yintercept = V2), color = "red", size = 2) +
#     {if(prefix == "c_fw_no") geom_vline(xintercept = 9, color = "red", size = 2)} +
#     {if(prefix == "c_fw_no") geom_vline(xintercept = 9, color = "green", size = 2, linetype = "dashed")} +
#     {if(prefix == "c_bw_no") geom_vline(xintercept = 8, color = "red", size = 2)} +
#     {if(prefix == "c_bw_no") geom_vline(xintercept = 12, color = "green", size = 2)} +
#     # {if(prefix == "bw") geom_vline(xintercept = 12, color = "green", size = 2)} +
#     # {if(prefix == "bw") geom_rect(aes(xmin = 11, xmax = 12, ymin = -Inf, ymax = Inf), fill = "palegreen", alpha = 0.11)} +
#     geom_point(mapping = aes(x = V4, y = V2), size = 8) +
#     geom_line(mapping = aes(x = V4, y = V2), size = 2) +
#     facet_wrap(~V1, scales = "free", ncol = 1) +
#     theme_bw() +
#     theme(axis.text = element_text(size = 25),
#           axis.title = element_text(size = 40, face = "bold"),
#           strip.text = element_text(size = 35, face = "bold"),
#           panel.spacing = unit(1.5, "lines"),
#           axis.title.x = element_text(vjust = 0.5, hjust = 1),
#           panel.grid.major = element_line(colour = "black")) +
#     scale_x_continuous(breaks = 1:17, name = "Number of colors in a repaired graph") +
#     scale_y_continuous(name = element_blank())
# plotData
# ggsave("/home/ian/Dropbox (City College)/Research/PhD work/shared folders/PROTEIN-FOLDING/COLLABORATION-ANALYSIS/DAVID-FRANCESCO/FINAL-ANALYSIS-MAY-2021/IndicesSummary.png", plot = plotData, width = 1, height = 2, scale = 10)
# # # ggsave("IndicesSummary.png", plot = plotData, width = 5, height = 2.5, scale = 8)

ggplot(data = indices) +
  geom_hline(data = manualIndices, aes(yintercept = V2), color = "red", size = 2) +
  geom_point(mapping = aes(x = V4, y = V2), size = 8) +
  geom_line(mapping = aes(x = V4, y = V2), size = 2) +
  facet_wrap(~V1, scales = "free", ncol = 1) +
  scale_x_continuous(breaks = 1:17, name = "Number of colors in a repaired graph") +
  theme_bw()
      