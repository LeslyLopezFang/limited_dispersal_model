import sys
import random as rd
import numpy as np
import pandas as pd
import os.path
import argparse
import itertools
from collections import defaultdict
import limited_dispersal_model_functions as db

parser = argparse.ArgumentParser(description="Run distance-based model.")

parser.add_argument('-N',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('-r',action="store",dest="simulation_number",type=int,help="What simulation number is this?")
parser.add_argument('-L',action="store",dest="L",type=int,help="L is the number of loci simulated by one run of the limited dispersal model.")
parser.add_argument('-s',action="store",dest="total_sample_size_of_tij",type=int,help="What is the total sample size of all tij's (ie (N*(N-1))/2 * sample size of one tij?")
parser.add_argument('-alpha',action="store",nargs='?',type=float,dest="alpha",help="Alpha is strength of isolation.")
parser.add_argument('--standard',action = "store_true",default=False,help="Run standard Wright Fisher model with the db-model.")
parser.add_argument('--save_matrix',action="store_true",default=False,help="Save all coalescent times (G matrix) from run.")
parser.add_argument('--verbose',action="store_true",default=False,help="Print statements as simulation runs in 'outlog.txt'")
param = parser.parse_args()

#Set alpha if standard true
if param.standard:
  param.alpha=1.0

if param.alpha:
  pass
else:
  print("Parameter '-alpha' required.\n")
  sys.exit(1)

#Outfile paths and names
matrix_path="../data/arrays/";tree_path="../data/trees/";tmrca_path="../data/tmrca/"
G_matrix_file="{}coalescent_times_matrix__db-model_N-{}_alpha-{}_L-{}_runs-{}".format(matrix_path,param.N
             ,param.alpha,param.L,param.total_sample_size_of_tij)
matrix_file="{}coalescent_times_matrix_sum_db-model_N-{}_alpha-{}_L-{}_runs-{}".format(matrix_path,param.N
             ,param.alpha,param.L,param.total_sample_size_of_tij)
newick_file="{}newick_tree_db-model_N-{}_alpha-{}_L-{}_runs-{}.tre".format(tree_path,param.N
             ,param.alpha,param.L,param.total_sample_size_of_tij)
tmrca_file="{}tmrca_list_db-model_N-{}_alpha-{}_L-{}_runs-{}.txt".format(tmrca_path,param.N
             ,param.alpha,param.L,param.total_sample_size_of_tij)
if param.verbose:
  outlog=open("outlog.txt","w")

t = 0
parents = {}
#Keep track of branch lengths
node_time = defaultdict(int)
#Keep track of children of the parents
children = defaultdict(list)
#Keep track of all coalescent events
coalescent_events=defaultdict(lambda:None) #coalescent_events=defaultdict(list)
number_of_coalescent_events=0
#Keep track of nodes or tips
nodes=defaultdict(lambda:None, {k:k for k in range(param.N)})
string_dicc=defaultdict(lambda:None)
#Get the tips for every internal node
tips_dict=defaultdict(list, {k:[k] for k in range(param.N)})
#Keep list of a list of every individual at each time
individuals_list = [range(param.N)]
#Keep track of mrca for starting sample
mrca_list = individuals_list[t]
#Keep track of the times a coalescent event occurs
coalescent_time = []
#Keep track of triple coalescent events
#triple_coalescence=defaultdict(None)
#Multiply branch lengths for triple coalescence
scalar_branch_length=4
triple_coalescence_branch_length=0.25*scalar_branch_length

#Get position of all starting individuals
pos = range(1,param.N+1)
#Get dij matrix: distance between individuals i and individuals j
dij = db.get_dij_matrix(N=param.N,pos=pos)
#Get probability matrix for an individual in each position to get a parent based on position
pij = db.probability_picking_parent(param.alpha,dij,param.N)

if param.verbose:
  outlog.write("Starting simulation\n")
  outlog.flush()
  #To ensure the outlog is written
  os.fsync(outlog.fileno())

while(len(mrca_list)!= 1 ):
  #Get individuals we are finding parents for in the time generation t
  individuals = individuals_list[t]
  #Get possible parents
  possible_parents = range(individuals[-1]+1,individuals[-1]+1+param.N)
  individuals_list.append(possible_parents)
  t+=1
  #Pick parent for each individual at that generation
  for ind in mrca_list:
    i = individuals.index(ind)
    #Get probability density for position
    parent_index = np.random.choice(param.N,p=pij[i])
    children[possible_parents[parent_index]].append(ind)
    #Update nodes dictionary
    node = ind
    while not node in coalescent_events.values() and not node in range(param.N):
      node = nodes[node]
    nodes[possible_parents[parent_index]]=node
  mrca_list = [parent for parent in possible_parents if parent in children.keys()]
  #Was there a coalescent event? Only look at the current time step
  for coalesced_node in [key for key,values in children.items() if (len(values) > 1) & (key in mrca_list)]:
    #Update nodes dictionary
    nodes[coalesced_node] = coalesced_node
    coalescent_time.append(t)
    #Keep track of coalescent events
    number_of_coalescent_events+=1
    coalescent_events[number_of_coalescent_events]=coalesced_node
    node_time[coalesced_node] = t*scalar_branch_length
    #Keep track of tips
    tips=np.asarray([nodes[child] for child in children[coalesced_node]])
    #All offspring are tips
    if (tips < param.N).all():
      tips_dict[coalesced_node]+=list(tips)
    else:
      #Offspring include inner node
      updated_tips=list(tips[tips<param.N])
      for inner_node in tips[tips>=param.N]:
        updated_tips+=tips_dict[inner_node]
      tips_dict[coalesced_node]+=updated_tips
#    #Was there a triple coalescent event?
#    if(len(children[coalesced_node]) > 2):
#      triple_coalescence[number_of_coalescent_events]=[t,coalesced_node,children[coalesced_node]]

##If there is a triple coalescent, list the information
##if(param.multiple & len(triple_coalescence.keys())>0 ):
#if triple_coalescence:
#  with open(triple_coalescence_file,"a+") as fout:
#    for key, values in triple_coalescence.items():
##    fout.write("\t".join([str(N),str(param.alpha),str(param.run),str(key),str(values[0]),str(values[1]),str(values[2][0]),str(values[2][1]),str(values[2][2])])+"\n")
#      fout.write("\t".join([str(N),str(param.alpha),str(param.simulation_number),str(param.run),str(key),str(values[0]),str(values[1]),str(values[2])])+"\n")
##else:
##  with open(triple_coalescence_file,"w") as fout:
##    fout.write("\t".join(["N","alpha","run","coalesced event","time","coalesced node","child1","child2","child3"])+"\n")
##    for key, values in triple_coalescence.items():
##      fout.write("\t".join([str(N),str(param.alpha),str(param.run),str(key),str(values[0]),str(values[1]),str(values[2][0]),str(values[2][1]),str(values[2][2])])+"\n")

###############################################################################################################
#Ouput newick tree
#Prepare Newick Tree
for node_key in sorted(coalescent_events.keys()):
  inner_node=coalescent_events[node_key]
  children_list=[nodes[key] for key in children[inner_node]] #Get the tip or the next inner node (coalesced node)
  string_dicc[inner_node] = db.string_node(children_list,inner_node,string_dicc,node_time
                      ,scalar_branch_length,triple_coalescence_branch_length)

#Output newick tree
if os.path.exists(newick_file):
  fout_tree = open(newick_file,"a+")
  fout_tree.write("\n")
else:
  fout_tree = open(newick_file,"w")
outstring_tree= string_dicc[sorted(coalescent_events.values())[-1]]+";"
fout_tree.write(outstring_tree)
fout_tree.close()


###########################################################################################################
#				G matrix: matrix of coalescent times tij
#Check if a tij is being sampled for this simulation/coalescent tree
if (param.simulation_number > param.total_sample_size_of_tij):
  pass
else:
  #Get the tij to sample from this coalescent tree using the simulation number
  run=1;num=int(np.floor(param.L/((param.N*(param.N-1))/2)))
  sampling_dict={}
  for pair in itertools.combinations(range(param.N),2):
    for l in range(run,run+num):
      sampling_dict[l]=pair
    run+=num
  i=sampling_dict[param.simulation_number][0];j=sampling_dict[param.simulation_number][1]
  #Check if matrix_file exists to get G_sum or initiate G_sum
  if os.path.exists("{}.npy".format(matrix_file)):
    G_sum=np.load("{}.npy".format(matrix_file))
  else:
    G_sum=np.zeros((param.N,param.N))
  G=np.zeros((param.N,param.N))
  #Go through coalescent events to get tij's for this simulation
  for coalescent_event in range(1,len(coalescent_events.keys())+1):
    node_combinations=[]
    coalesced_nodes=np.asarray([nodes[child] for child in children[coalescent_events[coalescent_event]]])
    if (coalesced_nodes < param.N).all():
      #Get all combinations of tips
      node_combinations=list(itertools.combinations(coalesced_nodes,2))
    else:
      #Get tips for inner nodes
      updated_coalesced_nodes=[tips_dict[node] for node in coalesced_nodes]
      #Was there a multiple coalescence?
      if(len(coalesced_nodes)>2):
        for combo in itertools.combinations(updated_coalesced_nodes,2):
          node_combinations+=list(itertools.product(combo[0],combo[1]))
      else:
        #Only two nodes coalesced
        node_combinations+=list(itertools.product(updated_coalesced_nodes[0],updated_coalesced_nodes[1]))
    #Update matrix with the time of coalescence for every tip pair wise combination
    for node in node_combinations:
      G[min(node)][max(node)]+=node_time[coalescent_events[coalescent_event]]/scalar_branch_length
  #Update G matrix with the sample tij i.e. coalescent time for i and j
  G_sum[i][j]+=G[i][j]
  np.save(matrix_file,G_sum)
  #Save G matrix for simulation if turned on
  if param.save_matrix:
    np.save(G_matrix_file,G)

#######################################################################################################################
#					TMRCA
#Update tmrca_file with tmrca from this simulation run
tmrca=node_time[mrca_list[0]]/scalar_branch_length
if(os.path.exists(tmrca_file)):
  with open(tmrca_file,"a+") as fout:
    fout.write("\n"+str(tmrca))
else:
  with open(tmrca_file,"w") as fout:
    fout.write(str(tmrca))

################################################################################################################
if param.verbose:
  print("Finished simulation {} for N={} with alpha={} and L={}".format(
        param.simulation_number,param.N, param.alpha, param.L))
