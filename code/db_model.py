import sys
import random as rd
import numpy as np
import pandas as pd
import os.path
import argparse
import db_model_functions as db

parser = argparse.ArgumentParser(description="Run distance-based model.")
subparsers=parser.add_subparsers(help="sub-command help")

parser.add_argument('-N',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('-r',action="store",dest="run",help="What run number is this?")
parser.add_argument('--db_model',action="store_true",default=False,help="Run db-model.")
parser.add_argument('--binomial',action="store_true",default=False,help="Run model with binomial distribution.")
parser.add_argument('--poisson',action="store_true",default=False,help="Run model with poisson distribution.")
#per generation outfile
parser.add_argument('--per_generation',action="store_true",default=False,help="Output file with number of individuals per generation.")

#Create parser for binomial distribution
parser_bi = subparsers.add_parser('bi',help='binomial distribution help')
#parser_bi.add_argument('-n',type=int,action="store",dest="binom_n",help="Size/Number of trials (n) for binomial distribution.")
#parser_bi.add_argument('-p',type=float,action="store",dest="binom_p",help="Probability of a success for binomial distribution.")
parser_bi.add_argument('--binomial',action = "store_true", default=False,help="Use binomial distribution to pick the parent's location.")
parser_bi.set_defaults(binomial=True)

#Creat subparser for poisson distribution
parser_poisson = subparsers.add_parser('poisson',help="Poisson distribution help.")
#parser_poisson.add_argument('-mu',type=int,action="store",dest="poisson_mu",help="Mu for poisson distribution.")
parser_poisson.set_defaults(poisson=True)

#Creat subparser for linear db-model
parser_db = subparsers.add_parser('db',help="Linear and exponential db-model help.")
parser_db.add_argument('-alpha',action="store",dest="alpha",type=int,help="Alpha is strength of isolation.")
parser_db.add_argument('-type',action="store",dest="type",choices=['linear-db','exponential-db'],help="Probability of choosing parent: 'linear-db' or 'exponential-db'.")
parser_db.add_argument('--standard',action = "store_true",default=False,help="Run standard Wright Fisher model with the db-model.")
parser_db.set_defaults(db_model=True)

param = parser.parse_args()
N = param.N

#if(param.standard):
#  if(param.type=="linear"):
#    param.alpha = 0
#  elif(param.type=="exponential"):
#    param.alpha = 1
#  else:
#    print("Specify model type 'linear' or 'exponential'")
#    sys.exit()

#Outfile
outfile_path = "../Results/"
if(param.db_model):
  outfile = outfile_path +param.type+"_db-model_"+str(N)+"-N_"+str(param.alpha)+"-alpha.csv"
elif(param.binomial):
  outfile = outfile_path + "binomial_distribution_db-model_"+str(N)+"-N.csv"
elif(param.poisson):
  outfile = outfile_path + "poisson_distribution_db-model_"+str(N)+"-N.csv"

#per generation
if(param.per_generation):
  per_generation_outfile = outfile_path + param.type+"_db-model_"+str(N)+"-N_"
  per_generation_outfile += str(param.alpha)+"-alpha_ind-per-generation.csv"

t = 0
parents = {}
#Keep list of a list of every individual at each time
individuals_list = [range(N)]
#Keep track of mrca for starting sample
mrca = individuals_list[t]
#Keep track of the times a coalescent event occurs
coalescent_time = []

##For linear and exponential db-model
if(param.db_model):
  #Get position of all starting individuals
  pos = range(1,N+1)
  #Get dij matrix: distance between individuals i and individuals j
  dij = db.get_dij_matrix(N=N,pos=pos)
  #Get probability matrix for an individual in each position to get a parent based on position
  pij = db.probability_picking_parent(param.alpha,dij,N,param.type)

data=[]

while(len(mrca)!= 1 ):
  #Get individuals we are finding parents for in the time generation t
  individuals = individuals_list[t]
  #Get possible parents
  possible_parents = range(individuals[-1]+1,individuals[-1]+1+N)
  individuals_list.append(possible_parents)
  for i in range(N):
    ind = individuals[i]
    #For binomial distribution: need sample size (N) and location of the offspring (index i)
    if(param.binomial):
      parent_index = db.choose_parent_location_binomial(N,i)
    #For poisson distribution: need sample size (N) and location of the offspring (index i)
    if(param.poisson):
      parent_index = db.choose_parent_location_poisson(N,i)
    #For db-model
    if(param.db_model):
      parent_index = np.random.choice(N,p=pij[i])
    parent = possible_parents[parent_index]
    parents[ind] = parent
  allmrca = [parents[m] for m in mrca]
  mrca = list(set(allmrca))
  if(len(mrca)!=len(allmrca)):
    coalescent_time.append(t)
  t+= 1
  data.append((param.run,t,len(mrca)))

#Output tmrca or t
#db-model
if(param.db_model):
  header = "model\truns\tN\talpha\ttmrca\n"
  outline = param.type+"\t"+str(param.run)+"\t"+str(N)+"\t"+str(param.alpha)+"\t"+str(t)
elif(param.binomial):
  header = "model\t\runs\t\N\tmrca\n"
  outline = "binomial\t"+str(param.run)+"\t"+str(N)+"\t"+str(t)
elif(param.poisson):
  header  = "model\t\runs\t\N\tmrca\n"
  outline = "poisson\t"+str(param.run)+"\t"+str(N)+"\t"+str(t)

if(os.path.exists(per_generation_outfile)):
  fout = open(outfile,"a+")
  fout.write(outline+"\n")
  fout.close()
else:
  fout = open(outfile,"w")
  fout.write(header)
  fout.write(outline+"\n")
  fout.close()

#second outfile to store both number of individuals and time in generations
if(param.per_generation):
  columns=["run","time","number of individuals"]
  df = pd.DataFrame(data,columns=columns)
  if(os.path.exists(per_generation_outfile)):
    df.to_csv(per_generation_outfile,index=False,sep="\t",mode="a",header=False)
  else:
    df.to_csv(per_generation_outfile,index=False,sep="\t")

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
