#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 20:05:23 2021

@author: phillips

Version 4:
    conversion to gurobipy
    
Version 5:
    allows for an initial solution
    
Version 6:
    allows for a weighted objective
"""

import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import itertools as itools

##precision parameter
epsilon = .001


def readdata(fname,NumColors,FixedFlag,fixedfile):

    Colors = range(NumColors)
    
    GraphData = pd.read_csv(fname,sep=' ',index_col=[0,1],header=None,skiprows=1)
    
    
    EdgeDictRaw = GraphData.to_dict()[2]
    
    #remove duplicates
    tempedgelist = []
    for (i,j) in EdgeDictRaw.keys():
        if (i,j) not in tempedgelist and (j,i) not in tempedgelist:
            tempedgelist.append((i,j))
            
    EdgeDict = {(i,j):EdgeDictRaw[(i,j)] for (i,j) in tempedgelist}
    
    Edges,EdgeWeights = gp.multidict(EdgeDict)

    

    Nodes = []
    for tup in Edges:
        if tup[0] not in Nodes:
            Nodes.append(tup[0])
        if tup[1] not in Nodes:
            Nodes.append(tup[1])
            

    #this assumes an undirected graph so edges are in the network
    #only in one orientation.
    #WARNING: we are not data validating this assumption!
    NotEdges = []
    NodePairs = []
    for p in Nodes:
        for q in Nodes:
            if p != q:
                #not edges appear twice because of w_ijc and w_jic
                if (p,q) not in Edges and (q,p) not in Edges:
                    NotEdges.append((p,q))
                #node pairs only once
                if (p,q) not in NodePairs and (q,p) not in NodePairs:
                    NodePairs.append((p,q))

    NE_unique = list(set(NotEdges).intersection(set(NodePairs)))
    NotEdges = gp.tuplelist(NotEdges)
    NodePairs = gp.tuplelist(NodePairs)
    #big M parameter
    M = len(Nodes)
    
    if FixedFlag:
        ftable=pd.read_csv(fixedfile,index_col=0,sep=' ',header=None)        
        fdict = ftable.to_dict()[1]
    else:
        fdict = {}
        
    return Nodes,NodePairs,Colors,Edges,EdgeWeights,NotEdges,fdict,NE_unique,M


    ##create a MIP
def CreatePMIP(Nodes,Colors,NodePairs,FixedFlag,fdict,NotEdges,M,EdgeWeights,\
               env, InitialFlag=False,init_nc=[],init_sc=[],init_re=[],\
                   Wt_ObjFlag=False,ObjWeights=[]):
    
    
    pmip = gp.Model(name='Pseudo',env=env)
    
    #initialize node variables 
    nodecolor=pmip.addVars(Nodes,Colors,vtype=GRB.BINARY,name='nodecolor')
    samecolor=pmip.addVars(NodePairs,vtype=GRB.BINARY,name='samecolor')

    
    #constraints (6): every color used
    allcolors = pmip.addConstrs((gp.quicksum(nodecolor.select('*', c)) >= 1
                                 for c in Colors), name="allcolors")

    #constraints (7): all nodes assigned exactly one color
    allnodes = pmip.addConstrs((gp.quicksum(nodecolor.select(i, '*')) == 1
                                 for i in Nodes), name="allnodes")
    
    #constraints (8): node color assignments must respect the samecolor variable
    colorassignments = pmip.addConstrs((nodecolor[p,c]+nodecolor[q,c] <= \
                                         samecolor[p,q]+1 for (p,q) in NodePairs \
                                        for c in Colors),name='colorassignments')
    #constraints (10): node color assignments must respect the samecolor variable
    colorsame1 = pmip.addConstrs((nodecolor[p,c]-nodecolor[q,c] <= \
                                         1-samecolor[p,q] for (p,q) in NodePairs \
                                        for c in Colors),name='colorsame1')
    colorsame2 = pmip.addConstrs((nodecolor[q,c]-nodecolor[p,c] <= \
                                         1-samecolor[p,q] for (p,q) in NodePairs \
                                        for c in Colors),name='colorsame2')

    antisymmetry=[]
    #symmetry breaking constraints
    ColorMinus = range(max(Colors))
    antisymmetry = pmip.addConstrs((gp.quicksum(nodecolor.select('*',c)) >= \
                                         gp.quicksum(nodecolor.select('*',c+1)) \
                                         for c in ColorMinus), name='antisymmetry')


    #initialize pseudo edges        
    repair_edges = pmip.addVars(NotEdges,Colors,lb=0.0,ub=1.0,vtype=GRB.BINARY,\
                                    name='repair_edges')


    #repair_edges are directional   
    #big M enforced constraints              
    colors_agree_pq = pmip.addConstrs((gp.quicksum(repair_edges.select(p,'*',c)) + \
                                      sum(nodecolor[j,c]*EdgeWeights[i,j] for (i,j) in \
                                      EdgeWeights.keys() if i==p) + \
                                      sum(nodecolor[j,c]*EdgeWeights[j,i] for (j,i) in \
                                      EdgeWeights.keys() if i==p) - \
                                      (gp.quicksum(repair_edges.select(q,'*',c)) + \
                                      sum(nodecolor[j,c]*EdgeWeights[i,j] for (i,j) in \
                                      EdgeWeights.keys() if i==q) + \
                                      sum(nodecolor[j,c]*EdgeWeights[j,i] for (j,i) in \
                                      EdgeWeights.keys() if i==q))  \
                                      <= M*(1-samecolor[p,q]) \
                                      for (p,q) in NodePairs \
                                      for c in Colors),\
                                      name='colors_agree_pq')

    colors_agree_qp = pmip.addConstrs((-(gp.quicksum(repair_edges.select(p,'*',c)) + \
                                      sum(nodecolor[j,c]*EdgeWeights[i,j] for (i,j) in \
                                      EdgeWeights.keys() if i==p) + \
                                      sum(nodecolor[j,c]*EdgeWeights[j,i] for (j,i) in \
                                      EdgeWeights.keys() if i==p) ) + \
                                      (gp.quicksum(repair_edges.select(q,'*',c)) + \
                                      sum(nodecolor[j,c]*EdgeWeights[i,j] for (i,j) in \
                                      EdgeWeights.keys() if i==q) + \
                                      sum(nodecolor[j,c]*EdgeWeights[j,i] for (j,i) in \
                                      EdgeWeights.keys() if i==q))  \
                                      <= M*(1-samecolor[p,q]) \
                                      for (p,q) in NodePairs \
                                      for c in Colors),\
                                      name='colors_agree_qp')

        

    #only allow repair edges of the correct color
    only_correct_repairs = pmip.addConstrs((repair_edges[i,j,c]<=nodecolor[j,c]\
                                            for (i,j) in NotEdges for c in Colors), \
                                           name='only_correct_repairs')


    equal_repair_edges = pmip.addConstrs((gp.quicksum(repair_edges.select(p,q,'*')) \
                        - gp.quicksum(repair_edges.select(q,p,'*')) == 0\
                        for (p,q) in NodePairs),
                        name='equal_repair_edges')
    
    if FixedFlag:
        fixedcolors = set(fdict.values())
        fnL = []
        for i in fixedcolors:
            fnL.append({v for v in fdict.keys() if fdict[v]==i})
        
        fixedcons=[]                
        #ensure nodes in the same set are the same color
        for L in fnL:
            for p,q in itools.combinations(L,2):
                if (p,q) in samecolor.keys():
                    fixedcons.append(pmip.addConstr(samecolor[(p,q)]==1,name='fixed_eq'+str(p)+'_'+str(q)))
                else:
                    fixedcons.append(pmip.addConstr(samecolor[(q,p)]==1,name='fixed_eq'+str(p)+'_'+str(q)))
            
        #ensure nodes in different sets are different colors
        for K,L in itools.combinations(fnL,2):
            for p in K:
                for q in L:
                    if (p,q) in samecolor.keys():
                        fixedcons.append(pmip.addConstr(samecolor[(p,q)]==0,name='fixed_ne_'+str(p)+'_'+str(q)))
                    else:
                        fixedcons.append(pmip.addConstr(samecolor[(q,p)]==0,name='fixed_ne_'+str(p)+'_'+str(q)))
                    
    else:
        fixedcons=[]

    if Wt_ObjFlag:
        obj = gp.quicksum(ObjWeights[p,q,c]*repair_edges[p,q,c] for (p,q) in NotEdges for c in Colors)
    else:   
        obj = gp.quicksum(repair_edges[p,q,c] for (p,q) in NotEdges for c in Colors) 
    
    pmip.setObjective(obj,GRB.MINIMIZE)

    if InitialFlag:
        for i in Nodes:
            for c in Colors:
                nodecolor[i,c].Start = init_nc[i,c]
        
        for (i,j) in NodePairs:
            samecolor[i,j].Start = init_sc[i,j]
            
        for (i,j) in NotEdges:
            for c in Colors:
                repair_edges[i,j,c].Start = init_re[i,j,c]
        
    return pmip,nodecolor,samecolor,allcolors,allnodes,colorassignments,\
        fixedcons,antisymmetry,colorsame1,colorsame2,repair_edges,\
        colors_agree_pq,colors_agree_qp,only_correct_repairs,\
        equal_repair_edges

def UpdateNodeColor(nodecolor,Nodes,Colors):
    badsol = False
    sol = {(i,c) for i in Nodes for c in Colors if (abs(nodecolor[i,c].x-1)<epsilon)}
    if (len(sol) != len(Nodes)):
        print('Bad solution!')
        f = open('debugfile.txt','w')
        for i in Nodes:
            for c in Colors:
                print(f"nodecolor[{i},{c}].x = {nodecolor[i,c].x}",file=f)
        f.close()
        badsol = True
    return sol,badsol
    #return {(i,c) for i in Nodes for c in Colors if nodecolor[i,c].x==1}
    

def UpdateRepairEdges(repair_edges,NotEdges,Colors):
    return {(i,j,c,repair_edges[i,j,c].x) for (i,j) in NotEdges for c in Colors if \
                repair_edges[i,j,c].x > epsilon}


def run_one_pmip(fname,NumColors,FixedFlag,fixedfile,timelimit,InitFlag,\
                 init_nc,init_sc,init_re,Wt_ObjFlag=False,ObjWeights=[]):

    #create the inputs
    Nodes,NodePairs,Colors,Edges,EdgeWeights,NotEdges,fdict,NE_unique,M = \
        readdata(fname,NumColors,FixedFlag,fixedfile)
    
    #initialize an environment
    env = gp.Env()
    
    #create the model
    pmip,nodecolor,samecolor,allcolors,allnodes,colorassignments,\
        fixedcons,antisymmetry,colorsame1,colorsame2,repair_edges,\
        colors_agree_pq,colors_agree_qp,only_correct_repairs,\
        equal_repair_edges = \
        CreatePMIP(Nodes,Colors,NodePairs,FixedFlag,fdict,NotEdges,M,\
                   EdgeWeights,env,InitFlag,init_nc,init_sc,init_re,
                   Wt_ObjFlag,ObjWeights)
        
    #set the time limit
    pmip.setParam("TimeLimit", timelimit)

    #optimize
    pmip.optimize()

    #find the solution
    nodecolor_sol,bad_sol = UpdateNodeColor(nodecolor,Nodes,Colors)
    repair_edges_sol = UpdateRepairEdges(repair_edges,NotEdges,Colors)
    
    return pmip,nodecolor_sol,repair_edges_sol,bad_sol


def write_one_solution(fname,nc_sol,re_sol,obj,lb):
    f = open(fname,"w")
    
    print(f"Edge weight sum: {obj}",file=f)    
    print(f"Edge weight lower bound: {lb}",file=f)
    print("Node,Color",file=f)    
    
    
    for sol in nc_sol:
        node = sol[0]
        color = sol[1]
        print(f'{node},{color}',file=f)

    print("\n",end="",file=f)
    print("Edges added",file=f)

    for repair in re_sol:
        i = repair[0]
        j = repair[1]
        c = repair[2]
        v = repair[3]
        print(f"{i},{j},{c},{v}",file=f)

    f.close()
    
#runs one MILP and then writes the result to a file that incorporates
#the original input filename, the NumColors, and fixedfile
def RWmilp(fname,NumColors,FixedFlag,fixedfile,outpath,timelimit,\
           InitFlag=False,init_nc={},init_sc={},init_re={}):
    pmip,nc_sol,re_sol,bad_sol = run_one_pmip(fname,NumColors,FixedFlag,\
                                              fixedfile,timelimit,InitFlag,\
                                              init_nc,init_sc,init_re)
    

    write_one_solution(outpath,nc_sol,re_sol,pmip.ObjVal,pmip.ObjBound)
    
    return pmip,nc_sol,re_sol,bad_sol

##main calls

#pmip_test,nc_sol_test,re_sol_test,bad_sol_test = \
#RWmilp('test_pseudo.graph.txt',3,True,\
#       'test_pseudo_fixed.txt','test_pseudo_out.txt',1000)

#RWmilp('test_pseudo_comp1.graph.txt',1,False,\
#       'no.txt','test_comp1_c1_out.txt',1000)


pmip_test2,nc_sol_test2,re_sol_test2,bad_sol_test2 = RWmilp('bw_fw_data/fw_gap_jn.txt',\
                                                      9,True,\
                                                      'bw_fw_data/old_fw_gap_jn_balanced_nodes.txt',
                                                      'test_fw_9.out.txt',3600)


#for i in range(1,10):
#    ian_name="ian"+str(i)+".out.txt"
#    pmip_test2,nc_sol_test2,re_sol_test2,bad_sol_test2 = RWmilp('iantest.txt',\
#                                                      i,False,\
#                                                      'bw_fw_data/old_fw_gap_jn_balanced_nodes.txt',
#                                                      ian_name,3600)

 