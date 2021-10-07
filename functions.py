#from os import listdir
import os
import re
import pandas as pd
import numpy as np

from SymmetryAnalyzer import SymmetryAnalyzer

def parseMIPOutputFile(inputFile):
	inputData = pd.read_table(inputFile, sep = "\n", header = None, names = ["Data"])
	
	#nodes = inputData.loc[(np.flatnonzero(inputData['Data'] == "Node,Color")[0] + 1):
	#					  (np.flatnonzero(inputData['Data'] == "Edges added")[0] - 1)]
	#nodes = nodes['Data'].str.split(",", expand = True)
	#self.nodes = nodes.rename(columns={0: "Name", 1: "Color"})
	
	originalEdges = inputData.loc[(np.flatnonzero(inputData['Data'] == "Input File")[0] + 2):
								  (np.flatnonzero(inputData['Data'] == "Fixed File")[0] - 1)]
	originalEdges = originalEdges['Data'].str.split(" ", expand = True)
	originalEdges.columns = ["Source", "Target", "Weight"]
	originalEdges[["Source", "Target"]] = [sorted([row[0], row[1]]) for row in originalEdges.values]
	originalEdges = originalEdges.drop_duplicates()
	originalEdges["Color"] = "Original"
	
	repairedEdges = inputData.loc[(np.flatnonzero(inputData['Data'] == "Edges added")[0] + 1):
								  (np.flatnonzero(inputData['Data'] == "Input File")[0] - 1)]
	if not(repairedEdges.empty):
		repairedEdges = repairedEdges['Data'].str.split(",", expand = True)
		repairedEdges.columns = ["Source", "Target", "Id", "Weight"]
		repairedEdges = repairedEdges[["Source", "Target", "Weight"]]
		repairedEdges[["Source", "Target"]] = [sorted([row[0], row[1]]) for row in repairedEdges.values]
		repairedEdges = repairedEdges.drop_duplicates(subset = ["Source", "Target"])
		repairedEdges["Color"] = "Repaired"
		edges = pd.concat([originalEdges, repairedEdges])
	else:
		edges = originalEdges
	return(edges)

def readEdgeFile(inputFile, sep = "\t", header = None):
	edges = pd.read_table(inputFile, sep = sep, header = header)
	if len(edges.columns) == 2:
		edges["Weight"] = 1
	
	edges.columns = ["Source", "Target", "Weight"]
	edges["Color"] = "Original"
	return(edges)

def runFile(inputFile):
	print("Running ", inputFile)
	edges = parseMIPOutputFile(inputFile)
	symmetryAnalyzer = SymmetryAnalyzer(edges)
	symmetryAnalyzer.decompose()
	symmetryAnalyzer.calculate_indices()
	symmetryAnalyzer.print_all_info()

def runFile_to_output(inputFile, outputPrefix):
	print("Running ", inputFile)
	edges = parseMIPOutputFile(inputFile)
	symmetryAnalyzer = SymmetryAnalyzer(edges)
	symmetryAnalyzer.decompose()
	symmetryAnalyzer.calculate_indices()
	symmetryAnalyzer.print_all_info_to_files(outputPrefix)

def runFolder(folderToRun, outputFolder):
	files = os.listdir(folderToRun)
	files = [file for file in files if re.search(r"\.[a-z]{3}$", file)]
	
	os.makedirs(outputFolder, exist_ok = True)
	for file in files:
		prefix = re.sub(r"\.[a-z]{3}$", "", file)
		runFile_to_output(folderToRun + "/" + file, outputFolder + "/" + prefix)
