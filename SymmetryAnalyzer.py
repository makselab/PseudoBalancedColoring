import igraph as ig
import pynauty as pynauty
import pandas as pd
import sys as sys
from collections import defaultdict

from NormalSubgroup import *

class SymmetryAnalyzer:
	def __init__(self, edges, directed = False):
		self.edges = edges
		
		self.ig_graph = ig.Graph.TupleList(self.edges[["Source", "Target", "Color"]].itertuples(index = False),
										   directed = directed, edge_attrs = "color")
		self.pynauty_graph = self.get_pynauty_graph_from_ig_graph(self.ig_graph)
		
		self.generators, grpsize1, grpsize2, orbits, numorbits = pynauty.autgrp(self.pynauty_graph)
		if (len(self.generators) == 0):
				print("Graph has no automorphisms")
				sys.exit()

		self.permutations = [combinatorics.Permutation(generator) for generator in self.generators]
		
		self.nodes = pd.DataFrame({"Name": self.ig_graph.vs["name"],
								  "Orbit": orbits})
	
	def decompose(self):
		self.sectors = self.get_sectors_from_permutations(self.permutations)
		self.nodes["Sector"] = self.get_assigned_sectors(self.sectors, self.ig_graph)
		self.nodes["Orbit"] = self.nodes["Sector"].astype(str) + "_" + self.nodes["Orbit"].astype(str)
		
		self.NormalSubgroups = [NormalSubgroup(self.ig_graph, self.sectors[i],
											   self.get_permutations_on_sector_id(self.permutations,
																				  self.sectors, i))
								for i in range(len(self.sectors))]
		
		self.calculate_epsilons()
	
	def calculate_epsilons(self):
		for NormalSubgroup in self.NormalSubgroups:
			NormalSubgroup.calculate_epsilons()
		
		meanEpsilons = [NormalSubgroup.meanEpsilon for NormalSubgroup in self.NormalSubgroups]
		maxEpsilons = [NormalSubgroup.maxEpsilon for NormalSubgroup in self.NormalSubgroups]
		self.meanEpsilon = mean(meanEpsilons)
		self.maxEpsilon = max(maxEpsilons)
	
	# List of indices (update it every time you change anything as it's used twice in the print info section):
	# Fiedler value
	# Normalized Fiedler value
	# Eigen ratio
	###### experimental indices:
	# Random Walk Fiedler Value
	# Minimal degree
	# Maximal degree
	# Sum of degrees
	# Vertex connectivity
	# Edge connectivity
	def calculate_indices(self):
		self.calculate_fiedler_values()
		self.calculate_experimental_values()
	
	#############################################################################
	########################### Helping functions ###############################
	#############################################################################
	def get_pynauty_graph_from_ig_graph(self, ig_graph):
		# Represent igraph as adjacency list-of-dictionaries
		adjacency_dict = {}
		for source in ig_graph.vs.indices:
			adjacent_edges = ig_graph.es.select(_source = source)
			adjacent_nodes = []
			for edge in adjacent_edges:
				if(edge.source == source):
					adjacent_nodes.append(edge.target)
				else:
					adjacent_nodes.append(edge.source)
			adjacency_dict[source] = adjacent_nodes

		pynauty_graph = pynauty.Graph(number_of_vertices = ig_graph.vcount(),
									  directed = ig_graph.is_directed(),
									  adjacency_dict = adjacency_dict)
		
		return(pynauty_graph)

	def get_sectors_from_permutations(self, permutations):
		permutation_domains = []
		for permutation in permutations:
			permutation_domains.append([value for sublist in permutation.cyclic_form for value in sublist])
			
		return(list(self.merge_common(permutation_domains)))

	# this is the function taken from the internet, it would be better to re-write it
	def merge_common(self, lists):
		neigh = defaultdict(set)
		visited = set()
		for each in lists:
			for item in each:
				neigh[item].update(each)
		def comp(node, neigh = neigh, visited = visited, vis = visited.add):
			nodes = set([node])
			next_node = nodes.pop
			while nodes:
				node = next_node()
				vis(node)
				nodes |= neigh[node] - visited
				yield node
		for node in neigh:
			if node not in visited:
				yield sorted(comp(node))

	def get_permutations_on_sector_id(self, permutations, sectors, sector_id):
		is_permutation_on_sector = lambda permutation, sector_id : all([node in set(sectors[sector_id]) for node in permutation])

		permutations_to_return = []
		for permutation in permutations:
			permutation_cyclic = [value for sublist in permutation.cyclic_form for value in sublist]
			if(is_permutation_on_sector(permutation_cyclic, sector_id)):
				permutations_to_return.append(permutation)
		return(permutations_to_return)
		
	def get_assigned_sectors(self, sectors, ig_graph):
		sectorId = []
		for i in ig_graph.vs.indices:
			search_result = [i in sector for sector in sectors]
			if (any(search_result) == False):
				sectorId.append(-1)
			else:
				sectorId.append([idx + 1 for idx in range(len(search_result)) if search_result[idx] == True][0])
		
		trivialIdx = max(sectorId) + 1
		for idx in range(len(sectorId)):
			if sectorId[idx] == -1:
				sectorId[idx] = trivialIdx
				trivialIdx = trivialIdx + 1
		return(sectorId)
	#############################################################################
	########################## Calculating indices ##############################
	#############################################################################
	def calculate_fiedler_values(self):
		graph_laplacian = self.ig_graph.laplacian(normalized = False)
		graph_laplacian = np.array(graph_laplacian)
		eigenval, eigenvect = np.linalg.eig(graph_laplacian)
		# we only take the real part, because imaginary part is very small and
		# it seems like it's here due to a computational artifact
		eigenval = eigenval.real
		eigenval = np.sort(eigenval)
		self.fiedlerValue = eigenval[1]
		self.eigenRatio = eigenval[1] / eigenval[len(eigenval) - 1]
		
		graph_laplacian = self.ig_graph.laplacian(normalized = True)
		graph_laplacian = np.array(graph_laplacian)
		eigenval, eigenvect = np.linalg.eig(graph_laplacian)
		# we only take the real part, because imaginary part is very small and
		# it seems like it's here due to a computational artifact
		eigenval = eigenval.real
		eigenval = np.sort(eigenval)
		self.normalizedFiedlerValue = eigenval[1]
	
	def calculate_experimental_values(self):
		adjacency_matrix = np.array(self.ig_graph.get_adjacency().data)
		eigenval, eigenvect = np.linalg.eig(adjacency_matrix)
		eigenval = np.sort(eigenval)
		kappa = eigenval[len(eigenval) - 1] + 5
		
		graph_laplacian = self.ig_graph.laplacian(normalized = False)
		graph_laplacian = np.array(graph_laplacian)
		degreeSeq = np.array(self.ig_graph.degree())
		randomWalkLaplacian = graph_laplacian / degreeSeq[:, None]
		eigenval, eigenvect = np.linalg.eig(randomWalkLaplacian)
		# we only take the real part, because imaginary part is very small and
		# it seems like it's here due to a computational artifact
		eigenval = eigenval.real
		eigenval = np.sort(eigenval)
		self.randomWalkSecondEigenvalue = eigenval[1]
		
		degreeSeq = np.sort(np.array(self.ig_graph.degree()))
		self.minimalDegree = np.ndarray.min(degreeSeq)
		self.maximalDegree = np.ndarray.max(degreeSeq)
		self.sumDegrees = np.sum(degreeSeq)
		
		self.vertex_connectivity = self.ig_graph.vertex_connectivity()
		self.edge_connectivity = self.ig_graph.edge_connectivity()

	#############################################################################
	################################### Print ###################################
	#############################################################################
	def print_all_info(self):
		nodes = self.nodes.sort_values(by = ["Sector", "Orbit"])
		print(self.nodes)
		
		self.print_indices()
		print()
		self.print_normal_subgroups()
	
	def print_all_info_to_files(self, prefix):
		self.print_indices_to_file(prefix + "_indices.txt")
		self.print_normal_subgroups_to_file(prefix + "_normal_subgroups.txt")
		self.print_permutations_to_file(prefix + "_permutations.txt")
		self.print_nodes_to_file(prefix + "_nodes.csv")
		self.print_edges_to_file(prefix + "_edges.csv")
	
	def print_normal_subgroups(self):
		for NormalSubgroup in self.NormalSubgroups:
			NormalSubgroup.print()
			print()
	
	def print_normal_subgroups_to_file(self, fileName):
		with open(fileName, 'w') as filehandle:
			pass
		for NormalSubgroup in self.NormalSubgroups:
			NormalSubgroup.print_info_to_file(fileName)

	def print_permutations_to_file(self, fileName):
		with open(fileName, 'w') as filehandle:
			filehandle.writelines("%s\n" % self.ig_graph.vs["name"])
		for NormalSubgroup in self.NormalSubgroups:
			NormalSubgroup.print_raw_permutations_to_file(fileName)
	
	def print_nodes_to_file(self, fileName):
		nodes = self.nodes.sort_values(by = ["Sector", "Orbit"])
		nodes.to_csv(fileName, index = False)
	
	def print_edges_to_file(self, fileName):
		edges = self.edges
		edges["Type"] = "Undirected"
		edges.to_csv(fileName, index = False)
	
	def print_indices(self):
		print("Fiedler value = ", self.fiedlerValue)
		print("Normalized Fiedler value = ", self.normalizedFiedlerValue)
		print("Eigen ratio = ", self.eigenRatio)
		print("Mean epsilon = ", self.meanEpsilon)
		print("Max epsilon = ", self.maxEpsilon)
	
	def print_indices_to_file(self, fileName):
		with open(fileName, 'w') as filehandle:
			filehandle.writelines("Fiedler value,%s\n" % self.fiedlerValue)
			filehandle.writelines("Normalized Fiedler Value,%s\n" % self.normalizedFiedlerValue)
			filehandle.writelines("Eigen ratio,%s\n" % self.eigenRatio)
			filehandle.writelines("Mean epsilon,%s\n" % self.meanEpsilon)
			filehandle.writelines("Max epsilon,%s\n" % self.maxEpsilon)
			# experimental indices from here down
			filehandle.writelines("Random Walk Fiedler Value,%s\n" % self.randomWalkSecondEigenvalue)
			filehandle.writelines("Minimal degree,%s\n" % self.minimalDegree)
			filehandle.writelines("Maximal degree,%s\n" % self.maximalDegree)
			filehandle.writelines("Sum of degrees,%s\n" % self.sumDegrees)
			filehandle.writelines("Vertex connectivity,%s\n" % self.vertex_connectivity)
			filehandle.writelines("Edge connectivity,%s\n" % self.edge_connectivity)
