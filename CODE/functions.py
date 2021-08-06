from os import listdir
import re

import pandas as pd
import numpy as np
import igraph as ig
from collections import defaultdict

import pynauty as pynauty
import sympy.combinatorics as combinatorics

import sage.groups.perm_gps.permgroup as sage_permgroup
import sage.groups.perm_gps.permgroup_named as sage_named_groups

def get_pynauty_graph(ig_graph):
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
								  directed = False,
								  adjacency_dict = adjacency_dict)
	
	return(pynauty_graph)

def get_sectors(permutations):
	permutation_domains = []
	for permutation in permutations:
		permutation_domains.append([value for sublist in permutation.cyclic_form for value in sublist])
		
	return(list(merge_common(permutation_domains)))

# this is the function taken from the internet, it would be better to re-write it
def merge_common(lists):
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

def get_permutations_on_sector_id(permutations, sectors, sector_id):
	is_permutation_on_sector = lambda permutation, sector_id : all([node in set(sectors[sector_id]) for node in permutation])

	permutations_to_return = []
	for permutation in permutations:
		permutation_cyclic = [value for sublist in permutation.cyclic_form for value in sublist]
		if(is_permutation_on_sector(permutation_cyclic, sector_id)):
			permutations_to_return.append(permutation)
	return(permutations_to_return)

def sympy_permutations_to_sage_group(sympy_permutations):
	subgroup_perms_sage_form = []
	
	for i in range(len(sympy_permutations)):
		permutation = sympy_permutations[i].cyclic_form
		permutation_sage = []
		for j in range(len(permutation)):
			permutation_sage.append(tuple(permutation[j]))
		subgroup_perms_sage_form.append(permutation_sage)

	sage_group = sage_permgroup.PermutationGroup(subgroup_perms_sage_form)
	return(sage_group)

def classify_sympy_perm_group(sympy_permutations, domain_size):
	sage_group = sympy_permutations_to_sage_group(sympy_permutations)

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

def get_assigned_sectors(sectors, ig_graph):
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
