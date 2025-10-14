import numpy as np
import pandas as pd

import argparse
import os.path
import sys

#Ex: sys_line = 'expectation_tmrca.py -type linear --standard -n 130'
# sys.argv = sys_line.split()

parser = argparse.ArgumentParser(description="Get expectation of tmrca of distance-based model.")

parser.add_argument('-alpha',action="store",dest="alpha",type=int,help="Alpha is strength of isolation.")
parser.add_argument('-type',action="store",dest="type",help="Probability of choosing parent: 'linear' or 'exponential'.")
parser.add_argument('-n',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('--standard',action = "store_true",default=False,help="Run standard Wright Fisher model.")
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

infile_path = "/home/llopez84/qsb/model-collab/Data/"
outfile_path = "/home/llopez84/qsb/model-collab/Results/Tmrca/"
infile = infile_path +param.type+"_db-model_"+str(param.N)+"-N_"+str(param.alpha)+"-alpha.csv"
outfile = outfile_path + "expectations_"+param.type+"_db-model_"+str(param.alpha)+"-alpha.csv"

#Read in file
#file = open(infile,"r")
#lines = file.readlines()
#file.close()
##Get all tmrca
#times = [float(i.split()[0]) for i in lines]
df = pd.read_csv(infile,sep="\t")

#analytical expectation of tmrca
ana_exp = 2.0*param.N*(1.-(1./param.N))

#simulated expectation of tmrca
sim_exp = np.mean(df["tmrca"])
sim_sd = np.std(df["tmrca"])
sim_se = sim_sd / np.sqrt(param.N)

outfile_exists=os.path.exists(outfile)
fout = open(outfile,"a+")
if not outfile_exists:
  header = "population size\truns\tanalytical expectation\tsimulated expectation\t"
  header += "simulated standard deviation\tsimulated standard error\n"
  fout.write(header)
outline = str(param.N)+"\t"+str(np.shape(df["tmrca"])[0])+"\t"+str(ana_exp)+"\t"+str(sim_exp)
outline+= "\t"+str(sim_sd)+"\t"+str(sim_se)+"\n"
fout.write(outline)
fout.close()

print("Analytical expectation is "+str(ana_exp))
print("Simulated expectation is "+str(sim_exp))

