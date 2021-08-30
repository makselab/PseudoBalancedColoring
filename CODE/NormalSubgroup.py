import igraph as ig
import numpy as np
from statistics import mean

import sympy.combinatorics as combinatorics

import sage.groups.perm_gps.permgroup as sage_permgroup
import sage.groups.perm_gps.permgroup_named as sage_named_groups

#from functions import *

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
		
		self.classification = self.classify_sympy_perm_group(self.sympy_permutations, len(self.sector))

		repairedEdgesBool = np.array(self.ig_graph.es["color"]) == "Repaired"
		repairedEdgesBool = list(np.flatnonzero(repairedEdgesBool))
		originalGraph = ig_graph.copy()
		originalGraph.delete_edges(repairedEdgesBool)
		self.original_adjacency_matrix = np.array(originalGraph.get_adjacency().data)
	
	def calculate_epsilon(self, permutation):
		permutation_array = permutation.array_form
		permutation_matrix = np.identity(len(permutation_array))[permutation_array,:]

		commutator = np.matmul(permutation_matrix, self.original_adjacency_matrix) - np.matmul(self.original_adjacency_matrix, permutation_matrix)

		epsilon = np.sum(np.absolute(commutator))
		return epsilon
		
	def calculate_epsilons(self):
		self.epsilons = []
		for permutation in self.sympy_permutations:
			self.epsilons.append(self.calculate_epsilon(permutation))

		self.meanEpsilon = mean(self.epsilons)
		self.maxEpsilon = max(self.epsilons)

	def sympy_permutations_to_sage_group(self, sympy_permutations):
		subgroup_perms_sage_form = []
		
		for i in range(len(sympy_permutations)):
			permutation = sympy_permutations[i].cyclic_form
			permutation_sage = []
			for j in range(len(permutation)):
				permutation_sage.append(tuple(permutation[j]))
			subgroup_perms_sage_form.append(permutation_sage)

		sage_group = sage_permgroup.PermutationGroup(subgroup_perms_sage_form)
		return(sage_group)

	def classify_sympy_perm_group(self, sympy_permutations, domain_size):
		sage_group = self.sympy_permutations_to_sage_group(sympy_permutations)

		group_name = "Unknown"
		for i in range(1, domain_size + 1):
			if(sage_group.is_isomorphic(sage_named_groups.SymmetricGroup(i))):
				group_name = "S" + str(i)
				return(group_name)
		
		for i in range(1, domain_size + 1):
			if(sage_group.is_isomorphic(sage_named_groups.DihedralGroup(i))):
				group_name = "D" + str(i)
				return(group_name)

		for i in range(1, domain_size + 1):
			if(sage_group.is_isomorphic(sage_named_groups.CyclicPermutationGroup(i))):
				group_name = "C" + str(i)
				return(group_name)

		for i in range(1, domain_size + 1):
			if(sage_group.is_isomorphic(sage_named_groups.AlternatingGroup(i))):
				group_name = "A" + str(i)
				return(group_name)
			
		return(group_name)

	def print(self):
		print("Class:\n\t", self.classification)
		print("Mean epsilon = \t", self.meanEpsilon)
		print("Max epsilon = \t", self.maxEpsilon)
		print("Sector:\n\t", self.sector_names)
		print("Orbits (", len(self.orbits), "):")
		for i in range(len(self.orbits)):
			print("\t", i ,":\t", self.orbits[i])
		
		print("Permutations:")
		for i in range(len(self.sympy_permutations)):
			permutation = self.sympy_permutations[i]
			print("\tEpsilon = ", self.epsilons[i], ":\t", [get_names_from_list_of_ids(self.ig_graph, nodes) for nodes in permutation.cyclic_form])
		#print(self.sympy_permutation_subgroup)

	def print_info_to_file(self, filename):
		with open(filename, 'a') as filehandle:
			filehandle.writelines("Class: %s\n" % self.classification)
			filehandle.writelines("Mean epsilon = %s\n" % self.meanEpsilon)
			filehandle.writelines("Max epsilon = %s\n" % self.maxEpsilon)
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
