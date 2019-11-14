from fractions import Fraction
import pandas as pd
import numpy as np
import os.path
import sys
import argparse
import db_model_functions as fun

parser = argparse.ArgumentParser(description="Run distance-based model.")

parser.add_argument('-alpha',action="store",dest="alpha",type=int,help="Alpha is strength of isolation.")
parser.add_argument('-type',action="store",dest="type",help="Probability of choosing parent: 'linear' or 'exponential'.")
parser.add_argument('-n',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('--standard',action = "store_true",default=False,help="Run standard Wright Fisher model.")
param = parser.parse_args()
N = param.N

outfile_path = "../Results/Analytical/"
coalescence_probability_outfile = outfile_path + param.type+"_probability_of_coalescence.csv"
expectation_coalescence_probability_outfile = outfile_path + param.type+"_expectation_probability_of_coalescence.csv"

if(param.standard):
  if(param.type=="linear"):
    param.alpha = 0
  elif(param.type=="exponential"):
    param.alpha = 1
  else:
    print("Specify model type 'linear' or 'exponential'")
    sys.exit()

#Get position of all starting individuals
pos = range(1,N+1)
#Get dij matrix: distance between individuals i and individuals j
dij = fun.get_dij_matrix(N=N,pos=pos)
#Get probability matrix for an individual in each position to get a parent based on position
pij = fun.probability_picking_parent(param.alpha,dij,N,param.type)


prob=[] #;counter=0
for i in np.arange(N-1):
  for j in np.arange(i+1,N):
#    print(np.dot(pij[n],pij[i]))
    prob.append((param.type,param.alpha,N,i+1,j+1,np.dot(pij[i],pij[j])))
#    counter += 1

#Get dataframe for probabilities
columns=["model_type","alpha","N","individual_i","individual_j","coalescence_probability"]
ds = pd.DataFrame(prob,columns=columns)

#Expectation of probability of a coalescent event
expectation = np.mean(ds["coalescence_probability"])

#Outfile to store the probability of a coalescence event between every individual
if(os.path.exists(coalescence_probability_outfile)):
  ds.to_csv(coalescence_probability_outfile,index=False,sep="\t",mode="a",header=False)
else:
  ds.to_csv(coalescence_probability_outfile,index=False,sep="\t")

#Outfile to store the expectation of the probability of a coalescent event
if(os.path.exists(expectation_coalescence_probability_outfile)):
  fout = open(expectation_coalescence_probability_outfile,"a+")
  outstring = param.type+"\t"+str(param.alpha)+"\t"+str(param.N)+"\t"+str(expectation) +"\n"
  fout.write(outstring)
else:
  fout = open(expectation_coalescence_probability_outfile,"w")
  header = "model_type\talpha\tN\texpectation\n"
  fout.write(header)
  outstring = param.type+"\t"+str(param.alpha)+"\t"+str(param.N)+"\t"+str(expectation) +"\n"
  fout.write(outstring)
fout.close()

outstring = param.type+": "+str(param.alpha)+"-alpha "+str(param.N)+"-N E[Prob(coalescence event)] is "
outstring += str(expectation) +"\n"
print(outstring)
