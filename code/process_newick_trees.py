import analysis_limited_dispersal_model_functions as af
import numpy as np
import pandas as pd
import argparse
import os

#Simulate SNPs from Newick trees and get M from simulated data
parser = argparse.ArgumentParser(description="Simulate L SNPs from L number of Newick trees from L runs of limited dispersal model and get M matrix from simulated data.")

parser.add_argument('-N',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('-L',action="store",dest="L",type=int,help="L is the number of loci simulated by one run of the limited dispersal model.")
parser.add_argument('-s',action="store",dest="total_sample_size_of_tij",type=int,help="What is the total sample size of all tij's (ie (N*(N-1))/2 * sample size of one tij?")
parser.add_argument('-alpha',action="store",type=float,dest="alpha",help="Alpha: strength of limited dispersal model. Wright Fisher has alpha of 1.0")
parser.add_argument('--verbose',action="store_true",default=False)
param = parser.parse_args()


tree_path="../data/trees/";results_path="../results/"
results_array_path="../results/arrays/"
newick_file="{}newick_tree_db-model_N-{}_alpha-{}_L-{}_runs-{}.tre".format(tree_path,param.N
             ,param.alpha,param.L,param.total_sample_size_of_tij)
M_matrix_file="{}M_from_snps_db-model_N-{}_alpha-{}_L-{}_runs-{}".format(results_array_path,param.N
             ,param.alpha,param.L,param.total_sample_size_of_tij)
estimated_coalescent_times_file="{}estimated_average_coalescent_times_db-model_N-{}_L-{}_runs-{}.csv".format(
             results_path,param.N,param.L,param.total_sample_size_of_tij)
pc_file="{}principal_components_db-model_N-{}_L-{}_runs-{}.csv".format(
             results_path,param.N,param.L,param.total_sample_size_of_tij)

#Model type based on alpha
model_type_dict=af.create_model_type_label()

if param.verbose:
  print("Starting SNPs simulation.\n")
#Consider making this into functions LESLY. Depends of if need total_branch_lengths
#Run SNP simulation function using newick trees
X_matrix,total_branch_lengths_list=af.get_allelic_matrix(newick_file=newick_file,N=param.N)

if param.verbose:
  print("SNPs simulation done.\n")

#Zero centre the X matrix
Z_matrix=X_matrix-X_matrix.mean(axis=1)[:,None]
#Get M matrix (NxN upper triangular matrix) using Equation LESLY
M_matrix=(1/param.L)*np.matmul(Z_matrix.T,Z_matrix)

#Save M matrix
#np.save(M_matrix_file,arr=M_matrix)

#########################################################################################################
#Compute estimated average coalescent times using M marix
estimates = af.compute_estimated_average_coalescent_times(M_matrix)

if param.verbose:
  print("Coalescent times estimated.\n")
#Make pandas dataframe with estimated coalescent times per distance between individuals i and j
tij_list=[];distance_list=[];distance_label_list=[]
counter=0
for i in range(param.N):
  for j in range(i+1,param.N):
    tij_list.append(estimates[counter])
    distance_list.append(abs(i-j))
    distance_label_list.append(af.create_distance_label(abs(i-j)))
    counter+=1

estimates_df=pd.DataFrame()
estimates_df["estimated_tij_average"]=tij_list
estimates_df["distance"]=distance_list
estimates_df["distance_label"]=distance_label_list
estimates_df["model_type"]=model_type_dict[param.alpha]
estimates_df["alpha"]=param.alpha
estimates_df=estimates_df[["model_type","alpha","estimated_tij_average","distance","distance_label"]]

if param.verbose:
  print("Estimates dataframe done.\n")
#########################################################################################################
#Compute principal components from M_matrix
pc1,pc2 = af.compute_principal_components_from_m_matrix(M_matrix)

if param.verbose:
  print("PCs computed.\n")

#Make pandas dataframe with top PC1 and PC2
pc_df=pd.DataFrame()
pc_df["model_type"]=[model_type_dict[param.alpha] for i in range(param.N)]
pc_df["alpha"]=[param.alpha for i in range(param.N)]
pc_df["PC1"]=pc1
pc_df["PC2"]=pc2
pc_df["Position"]=range(param.N)
pc_df["distance_label"]=pc_df.apply(lambda row: af.create_distance_label(int(row["Position"]))
                             ,axis=1)

if param.verbose:
  print("PC dataframe done.\n")
#########################################################################################################
#Write out estimates dataframe to estimate coalescent times results file
if os.path.exists(estimated_coalescent_times_file):
  estimates_df.to_csv(estimated_coalescent_times_file,index=False,sep="\t",mode="a",header=False)
else:
  estimates_df.to_csv(estimated_coalescent_times_file,index=False,sep="\t")

#Write out pc dataframe to results file
if os.path.exists(pc_file):
  pc_df.to_csv(pc_file,index=False,sep="\t",mode="a",header=False)
else:
  pc_df.to_csv(pc_file,index=False,sep="\t")
