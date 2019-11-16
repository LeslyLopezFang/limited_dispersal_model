import sys
import random as rd
import numpy as np
import pandas as pd
import os.path
import argparse
import db_model_functions as db

parser = argparse.ArgumentParser(description="Run distance-based model.")

parser.add_argument('-alpha',action="store",dest="alpha",type=int,help="Alpha is strength of isolation.")
parser.add_argument('-type',action="store",dest="type",help="Probability of choosing parent: 'linear' or 'exponential'.")
parser.add_argument('-n',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('--standard',action = "store_true",default=False,help="Run standard Wright Fisher model.")
parser.add_argument('-r',action="store",dest="run",help="What run number is this?")
#type = "exponential" #"linear" # or "exponential"
param = parser.parse_args()
N = param.N

if(param.standard):
  if(param.type=="linear"):
    param.alpha = 0
  elif(param.type=="exponential"):
    param.alpha = 1
  else:
    print("Specify model type 'linear' or 'exponential'")
    sys.exit()

#Outfile
outfile_path = "/home/llopez84/qsb/model-collab/Data/"
outfile = outfile_path +param.type+"_db-model_"+str(N)+"-N_"+str(param.alpha)+"-alpha.csv"
per_generation_outfile = outfile_path + param.type+"_db-model_"+str(N)+"-N_"
per_generation_outfile += str(param.alpha)+"-alpha_ind-per-generation.csv"
ind_name_outfile = outfile_path + "positions_individuals_mrca_"+param.type+"_db-model_"+str(param.alpha)+"-alpha_"+str(N)+"-N.csv"

t = 0
parents = {}
#Keep list of a list of every individual at each time
individuals_list = [range(N)]
#Keep track of mrca for starting sample
mrca = individuals_list[t]
#Keep track of the times a coalescent event occurs
coalescent_time = []

#Get position of all starting individuals
pos = range(1,N+1)
#Get dij matrix: distance between individuals i and individuals j
dij = db.get_dij_matrix(N=N,pos=pos)
#Get probability matrix for an individual in each position to get a parent based on position
pij = db.probability_picking_parent(param.alpha,dij,N,param.type)

data=[]

#Open individuals outfile
if(os.path.exists(ind_name_outfile)):
  individuals_file = open(ind_name_outfile,"a+")
else:
  individuals_file = open(ind_name_outfile,"a+")
  header = "model\tN\talpha\ttime\tsample\tindividual\tindividual_position\n"
  individuals_file.write(header)

while(len(mrca)!= 1 ):
  #Get individuals we are finding parents for in the time generation t
  individuals = individuals_list[t]
  #Get possible parents
  possible_parents = range(individuals[-1]+1,individuals[-1]+1+N)
  individuals_list.append(possible_parents)
  for i in range(N):
    ind = individuals[i]
    parent_index = np.random.choice(N,p=pij[i])
    parent = possible_parents[parent_index]
    parents[ind] = parent
  allmrca = [parents[m] for m in mrca]
  mrca = list(set(allmrca))
  if(len(mrca)!=len(allmrca)):
    coalescent_time.append(t)
  t+= 1
  data.append((param.run,t,len(mrca)))
#  #individuals_outfile info
#  if(len(mrca)<=3):
#    for mrca_ind in mrca:
#      individuals_outline = param.type+"\t"+str(param.N)+"\t"+str(param.alpha)+"\t"+str(t)+"\t"+str(len(mrca))+"\t"+str(mrca_ind)
#      individuals_outline += "\t"+str(mrca_ind%param.N)+"\n"
#      individuals_file.write(individuals_outline)

#Close individuals outfile
individuals_outline = param.type+"\t"+str(param.N)+"\t"+str(param.alpha)+"\t"+str(t)+"\t"+str(len(mrca))+"\t"+str(mrca[0])
individuals_outline += "\t"+str(mrca[0]%param.N)+"\n"
individuals_file.write(individuals_outline)

individuals_file.close()


##Output tmrca or t
#if(os.path.exists(per_generation_outfile)):
#  fout = open(outfile,"a+")
#  outline = param.type+"\t"+str(param.run)+"\t"+str(N)+"\t"+str(param.alpha)+"\t"+str(t)
#  fout.write(outline+"\n")
#  fout.close()
#else:
#  fout = open(outfile,"w")
#  header = "model\truns\tN\talpha\ttmrca\n"
#  fout.write(header)
#  outline = param.type+"\t"+str(param.run)+"\t"+str(N)+"\t"+str(param.alpha)+"\t"+str(t)
#  fout.write(outline+"\n")
#  fout.close()

##second outfile to store both number of individuals and time in generations
#columns=["run","time","number of individuals"]
#df = pd.DataFrame(data,columns=columns)
#if(os.path.exists(per_generation_outfile)):
#  df.to_csv(per_generation_outfile,index=False,sep="\t",mode="a",header=False)
#else:
#  df.to_csv(per_generation_outfile,index=False,sep="\t")

##To go through every individual one by one
#range(t*N)


#for ind in individuals_list[0]:
#  

##In case I need individuals and corresponding parents
#parent_data = []
#for key in parents:
#  parent_data.append((key,parents[key]))
#pf = pd.DataFrame(parent_data,columns=["individual","parent"])

print("Ran "+param.type+" model for "+str(param.alpha)+" with N = "+str(N))
