#Will add more parameters: subpopulation size, size of the chasm between populations.
import sys
import random as rd
import numpy as np
import pandas as pd
import os.path
import argparse
from collections import defaultdict
import db_model_functions as db

parser = argparse.ArgumentParser(description="Run distance-based model.")
subparsers=parser.add_subparsers(help="sub-command help")

parser.add_argument('-N',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('-r',action="store",dest="run",help="What run number is this?")
parser.add_argument('--db_model',action="store_true",default=False,help="Run db-model.")
parser.add_argument('--multiple',action="store_true",default=False,help="Run multiple runs of the model to get an output file with multiple trees.")
parser.add_argument('-alpha',action="store",dest="alpha",type=int,help="Alpha is strength of isolation.")
parser.add_argument('-type',action="store",dest="type",choices=['linear-db','exponential-db'],help="Probability of choosing parent: 'linear-db' or 'exponential-db'.")
#per generation outfile
parser.add_argument('--per_generation',action="store_true",default=False,help="Output file with number of individuals per generation.")

#Creat subparser for linear db-model
parser_db = subparsers.add_parser('db',help="Linear and exponential db-model help.")
parser_db.add_argument('--standard',action = "store_true",default=False,help="Run standard Wright Fisher model with the db-model.")
parser_db.set_defaults(db_model=True)

#Creat subparser for two populations for db-model
parser_two=subparsers.add_parser("mp",help="Two populations within db-model help.")
parser_two.add_argument("barrier_size",action="store",dest="barrier_size",type=int,help="Length of barrier in x-axis between two populations.")
parser_two.add_argument("--multiple_populations",action="store_true",default=True,help="Run db-model with two populations.")
param = parser.parse_args()
N = param.N

if(param.multiple_populations):
  outfile_path = "/home/llopez84/qsb/model-collab/Data/two_populations/"
  tree_path = "/home/llopez84/qsb/model-collab/results/trees/two_populations/"
  triple_coalescence_file=outfile_path+"two_populations_multiple_coalescence_db-model.csv"
  outfile=outfile_path+"two_populations_tmrca_db-model.csv"
  if(param.multiple):
    newick_file = tree_path + "newick_tree_"+param.type+"_db-model_N-"+str(N)+"_alpha-"+str(param.alpha)+"_number-of-runs-"+str(param.run)+".tre"
  else:
    newick_file = tree_path + "newick_tree_"+param.type+"_db-model_N-"+str(N)+"_alpha-"+str(param.alpha)+"_run-"+str(param.run)+".tre"
else:
  outfile_path = "/home/llopez84/qsb/model-collab/Data/"
  tree_path = "/home/llopez84/qsb/model-collab/results/trees/"
  triple_coalescence_file=outfile_path+"multiple_coalescence_db-model.csv"
#  triple_coalescence_file=outfile_path+"triple_coalescence_"+param.type+"_db-model_"+str(N)+"-N_"+str(param.alpha)+"-alpha.csv"
#  outfile = outfile_path +param.type+"_db-model_"+str(N)+"-N_"+str(param.alpha)+"-alpha.csv"
  outfile = outfile_path + "tmrca_db-model.csv"
  if(param.multiple):
    newick_file = tree_path + "newick_tree_"+param.type+"_db-model_N-"+str(N)+"_alpha-"+str(param.alpha)+"_number-of-runs-"+str(param.run)+".tre"
  else:
    newick_file = tree_path + "newick_tree_"+param.type+"_db-model_N-"+str(N)+"_alpha-"+str(param.alpha)+"_run-"+str(param.run)+".tre"

#per generation
if(param.per_generation):
  per_generation_outfile = outfile_path + param.type+"_db-model_"+str(N)+"-N_"
  per_generation_outfile += str(param.alpha)+"-alpha_ind-per-generation.csv"
outlog=open("outlog.txt","w")

###The number of populations N is constant for both populations
##For linear and exponential db-model: pmf
if param.multiple_populations:
  n=param.N
  N = (2*param.N)+param.barrier_size
  #Get position of all starting individuals
  pos = range(1,N+1)
  #Get dij matrix: distance between individuals i and individuals j
  dij = db.get_dij_matrix(N=N,pos=pos)
  #Get probability matrix for an individual in each position to get a parent based on position
  pij=db.probability_picking_parent_two_populations(alpha=param.alpha,dij=dij,N=n,barrier_size=param.barrier_size
                   ,model_type=param.type)
  #Only use position indices outside of the barrier
  use_position_indices=np.append(np.arange(0,n),np.arange(n+barrier_size,N))
else:
  #Get position of all starting individuals
  pos = range(1,N+1)
  #Get dij matrix: distance between individuals i and individuals j
  dij = db.get_dij_matrix(N=N,pos=pos)
  #Get probability matrix for an individual in each position to get a parent based on position
  pij = db.probability_picking_parent(alpha=param.alpha,dij=dij,N=N,model_type=param.type)
  #Use all indices within the population size
  use_position_indices=np.arange(0,N)




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
nodes=defaultdict(lambda:None, {k:k for k in range(N)})
string_dicc=defaultdict(lambda:None)
##Last node or tip corresponding to each position
#nodes_dicc = defaultdict(lambda:None, {k:k for k in range(N)})
#Keep list of a list of every individual at each time
individuals_list = [range(N)]
#Keep track of mrca for starting sample
mrca = [individuals_list[t][i] for i in use_position_index]
#Keep track of the times a coalescent event occurs
coalescent_time = []
#Keep track of triple coalescent events
triple_coalescence=defaultdict(None)
#Multiply branch lengths for triple coalescence
scalar_branch_length=4
triple_coalescence_branch_length=0.25*scalar_branch_length

data=[]

outlog.write("Starting simulation\n")
outlog.flush()
# typically the above line would do. however this is used to ensure that the file is written
os.fsync(outlog.fileno())
while(len(mrca)!= 1 ):
  #Get individuals we are finding parents for in the time generation t
  #Make sure the position for the individuals chosen is not in the space between populations (X'ed out)
  individuals=[individuals_list[t][i] for i in use_position_index]
  #Get possible parents
  possible_parents=range(individuals[-1]+1,individuals[-1]+1+N)
  individuals_list.append(possible_parents)
  t+=1
#  outlog.write("Time t = "+str(t)+" mrca number: +" + str(len(mrca))+"\n")
#  outlog.flush()
  # typically the above line would do. however this is used to ensure that the file is written
#  os.fsync(outlog.fileno())

  for ind in mrca:
    i = individuals.index(ind)
    parent_index = np.random.choice(N,p=pij[i])
    children[possible_parents[parent_index]].append(ind)
    #Update nodes dictionary
    node = ind
    while not node in coalescent_events.values() and not node in range(N):
      node = nodes[node]
    nodes[possible_parents[parent_index]]=node
  mrca = [parent for parent in possible_parents if parent in children.keys()]
  #Was there a coalescent event? Only look at the current time step
  for coalesced_node in [key for key,values in children.items() if (len(values) > 1) & (key in mrca)]:
    #Update nodes dictionary
    nodes[coalesced_node] = coalesced_node
    coalescent_time.append(t)
    #Keep track of coalescent events
    number_of_coalescent_events+=1
    coalescent_events[number_of_coalescent_events]=coalesced_node
    node_time[coalesced_node] = t*scalar_branch_length
    #Was there a triple coalescent event?
    if(len(children[coalesced_node]) > 2):
      triple_coalescence[number_of_coalescent_events]=[t,coalesced_node,children[coalesced_node]]
  data.append((param.run,t,len(mrca)))

outlog.write(str(len(triple_coalescence.keys()))+"\n")
outlog.close()

#If there is a triple coalescent, list the information
#if(param.multiple & len(triple_coalescence.keys())>0 ):
if triple_coalescence:
  with open(triple_coalescence_file,"a+") as fout:
    for key, values in triple_coalescence.items():
#    fout.write("\t".join([str(N),str(param.alpha),str(param.run),str(key),str(values[0]),str(values[1]),str(values[2][0]),str(values[2][1]),str(values[2][2])])+"\n")
      fout.write("\t".join([str(N),str(param.alpha),str(param.run),str(key),str(values[0]),str(values[1]),str(values[2])])+"\n")
else:
  with open(triple_coalescence_file,"w") as fout:
    fout.write("\t".join(["N","alpha","run","coalesced event","time","coalesced node","child1","child2","child3"])+"\n")
    for key, values in triple_coalescence.items():
      fout.write("\t".join([str(N),str(param.alpha),str(param.run),str(key),str(values[0]),str(values[1]),str(values[2][0]),str(values[2][1]),str(values[2][2])])+"\n")

#####UNDOOOOOOOOOOO
#Output newick tree
#How do you know you will get the first coalescent event?
#np.sort(list(coalescent_events.keys()))
for node_key in sorted(coalescent_events.keys()):
  inner_node=coalescent_events[node_key]
  children_list=[nodes[key] for key in children[inner_node]] #Get the tip or the next inner node (coalesced node)
  string_dicc[inner_node] = db.string_node(children_list,inner_node,string_dicc,node_time
                      ,scalar_branch_length,triple_coalescence_branch_length)


##Output stuff
##fout_tree = open(newick_file,"w")
#if(param.multiple):
#  fout_tree = open(newick_file,"a+")
#  fout_tree.write("\n")
#else:
#  fout_tree = open(newick_file,"w")
#outstring_tree= string_dicc[sorted(coalescent_events.values())[-1]]+";"
#fout_tree.write(outstring_tree)
#fout_tree.close()
#######################################################################################################################



#Output tmrca or t
#db-model
if(param.multiple_populations):
  header="model\talpha\truns\tN\tbarrier_size\ttmrca\n"
  outline = param.type+"\t"+str(param.alpha)+"\t"+str(param.run)+"\t"+str(N)+"\t"+str(param.barrier_size)+"\t"+str(t)
else:
  header = "model\truns\tN\talpha\ttmrca\n"
  outline = param.type+"\t"+str(param.run)+"\t"+str(N)+"\t"+str(param.alpha)+"\t"+str(t)
#elif(param.binomial):
#  header = "model\t\runs\t\N\tmrca\n"
#  outline = "binomial\t"+str(param.run)+"\t"+str(N)+"\t"+str(t)
#elif(param.poisson):
#  header  = "model\t\runs\t\N\tmrca\n"
#  outline = "poisson\t"+str(param.run)+"\t"+str(N)+"\t"+str(t)


if(os.path.exists(outfile)):
  fout = open(outfile,"a+")
  fout.write(outline+"\n")
  fout.close()
else:
  fout = open(outfile,"w")
  fout.write(header)
  fout.write(outline+"\n")
  fout.close()

print(param.N,param.run,param.alpha)




##second outfile to store both number of individuals and time in generations
#if(param.per_generation):
#  columns=["run","time","number of individuals"]
#  df = pd.DataFrame(data,columns=columns)
#  if(os.path.exists(per_generation_outfile)):
#    df.to_csv(per_generation_outfile,index=False,sep="\t",mode="a",header=False)
#  else:
#    df.to_csv(per_generation_outfile,index=False,sep="\t")

##To go through every individual one by one
#range(t*N)

#for ind in individuals_list[0]:
#  

#f#In case I need individuals and corresponding parents
#parent_data = []
#for key in parents:
#  parent_data.append((key,parents[key]))
#pf = pd.DataFrame(parent_data,columns=["individual","parent"])
#print("Ran "+param.type+" model for "+str(param.alpha)+" with N = "+str(N))
