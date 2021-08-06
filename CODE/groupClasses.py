import sympy.combinatorics as combinatorics
from functions import *
from statistics import mean

get_names_from_list_of_ids = lambda ig_graph, ids : [ig_graph.vs["name"][j] for j in ids]

class NormalSubgroup:
	def __init__(self, ig_graph, sector, permutations):
		self.ig_graph = ig_graph.copy()
		self.sector = sector
		self.sector_names = get_names_from_list_of_ids(self.ig_graph, sector)
		self.sympy_permutations = permutations
		self.sympy_permutation_subgroup = combinatorics.PermutationGroup(self.sympy_permutations)
		
		self.orbits = [orbit for orbit in self.sympy_permutation_subgroup.orbits() if len(orbit) > 1]
		self.orbits = [get_names_from_list_of_ids(ig_graph, list(orbit)) for orbit in self.orbits]
		
		self.classification = classify_sympy_perm_group(self.sympy_permutations, len(self.sector))

		repairedEdgesBool = np.array(self.ig_graph.es["color"]) == "Repaired"
		repairedEdgesBool = list(np.flatnonzero(repairedEdgesBool))
		originalGraph = ig_graph.copy()
		originalGraph.delete_edges(repairedEdgesBool)
		self.original_adjacency_matrix = np.array(originalGraph.get_adjacency().data)
	
	def print(self):
		print("Class:\n\t", self.classification)
		print("Mean epsilon = \t", self.meanEpsilon)
		print("Sector:\n\t", self.sector_names)
		print("Orbits (", len(self.orbits), "):")
		for i in range(len(self.orbits)):
			print("\t", i ,":\t", self.orbits[i])
		
		print("Permutations:")
		for i in range(len(self.sympy_permutations)):
			permutation = self.sympy_permutations[i]
			print("\tEpsilon = ", self.epsilons[i], ":\t", [get_names_from_list_of_ids(self.ig_graph, nodes) for nodes in permutation.cyclic_form])
		#print(self.sympy_permutation_subgroup)
		
	def get_epsilon(self, permutation):
		permutation_array = permutation.array_form
		permutation_matrix = np.identity(len(permutation_array))[permutation_array,:]

		commutator = np.matmul(permutation_matrix, self.original_adjacency_matrix) - np.matmul(self.original_adjacency_matrix, permutation_matrix)

		epsilon = np.sum(np.absolute(commutator))
		return epsilon
		
	def calculate_epsilons(self):
		self.epsilons = []
		for permutation in self.sympy_permutations:
			self.epsilons.append(self.get_epsilon(permutation))

		self.meanEpsilon = mean(self.epsilons)

	def print_info_to_file(self, filename):
		with open(filename, 'a') as filehandle:
			filehandle.writelines("Class: %s\n" % self.classification)
			filehandle.writelines("Mean epsilon = %s\n" % self.meanEpsilon)
			filehandle.writelines("Sector:\n\t%s\n" % self.sector_names)
			
			filehandle.writelines("Orbits (%i):\n" % len(self.orbits))
			for i in range(len(self.orbits)):
				filehandle.writelines("\t%i:\t%s\n" % (i, self.orbits[i]))
			
			filehandle.writelines("Permutations:\n")
			for i in range(len(self.sympy_permutations)):
				permutation = self.sympy_permutations[i]
				filehandle.writelines("\tEpsilon = %s:\t%s\n" % (self.epsilons[i], [get_names_from_list_of_ids(self.ig_graph, nodes) for nodes in permutation.cyclic_form]))
			filehandle.writelines("\n")
	
	def print_raw_permutations_to_file(self, filename):
		with open(filename, 'a') as filehandle:
			for j in range(len(self.sympy_permutations)):
				filehandle.writelines("%s\n" % self.sympy_permutations[j].array_form)
			filehandle.writelines("\n")
