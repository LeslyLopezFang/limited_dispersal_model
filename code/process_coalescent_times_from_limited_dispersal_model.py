import numpy as np
import pandas as pd
import argparse
import analysis_limited_dispersal_model_functions as af
import os

parser = argparse.ArgumentParser(description="Get average tij's from running limited dispersal model X (total sample size of tij -s) number of times.")

parser.add_argument('-N',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('-L',action="store",dest="L",type=int,help="L is the number of loci simulated by one run of the limited dispersal model.")
parser.add_argument('-s',action="store",dest="total_sample_size_of_tij",type=int,help="What is the total sample size of all tij's (ie (N*(N-1))/2 * sample size of one tij?")
parser.add_argument('-alpha',action="store",type=float,dest="alpha",help="Alpha: strength of limited dispersal model. Wright Fisher has alpha of 1.0")
parser.add_argument('--verbose',action="store_true",default=False)
param = parser.parse_args()


matrix_path="../data/arrays/";tmrca_path="../data/tmrca/"
results_path="../results/"

#Data files
G_sum_file="{}coalescent_times_matrix_sum_db-model_N-{}_alpha-{}_L-{}_runs-{}.npy".format(matrix_path,param.N
             ,param.alpha,param.L,param.total_sample_size_of_tij)
tmrca_file="{}tmrca_list_db-model_N-{}_alpha-{}_L-{}_runs-{}.txt".format(tmrca_path,param.N
             ,param.alpha,param.L,param.total_sample_size_of_tij)
#Results files
average_coalescent_times_file="{}average_coalescent_times_db-model_N-{}_L-{}_runs-{}.csv".format(results_path
                  ,param.N,param.L,param.total_sample_size_of_tij)
tmrca_dataframe_file="{}tmrca_dataframe_db-model_N-{}_L-{}_runs-{}.csv".format(tmrca_path
                  ,param.N,param.L,param.total_sample_size_of_tij)

#Get matrix of average coalescent times
G_avg=af.compute_average_coalescence_times_matrix(G_sum_file,param.total_sample_size_of_tij)

if param.verbose:
  print("Coalescent times averaged.\n")

#Model type based on alpha
model_type_dict=af.create_model_type_label()

tij_list=[];distance_list=[];distance_label_list=[]
for i in range(param.N):
  for j in range(i+1,param.N):
    #Works even if not all tij's sampled yet
    if G_avg[i][j] > 0:
      tij_list.append(G_avg[i][j])
      distance_list.append(abs(i-j))
      distance_label_list.append(af.create_distance_label(abs(i-j)))

df=pd.DataFrame()
df["tij_average"]=tij_list
df["distance"]=distance_list
df["distance_label"]=distance_label_list
df["alpha"]=param.alpha
df["model_type"]=model_type_dict[param.alpha]
df=df[["model_type","alpha","tij_average","distance","distance_label"]]

if param.verbose:
  print("Average coalescent times Data frame created.\n")

#Rearrange tmrca data to combine all model types to be ready to graph
tmrca_df=pd.read_csv(tmrca_file,names=["tmrca"])
tmrca_df["model_type"]=model_type_dict[param.alpha]

if param.verbose:
  print("TMRCA Data frame created.\n")

#Write average coalescent times dataframe
if os.path.exists(average_coalescent_times_file):
  df.to_csv(average_coalescent_times_file,index=False,sep="\t",mode="a",header=False)
else:
  df.to_csv(average_coalescent_times_file,index=False,sep="\t")

#Write Tmrca dataframe
if os.path.exists(tmrca_dataframe_file):
  tmrca_df.to_csv(tmrca_dataframe_file,index=False,sep="\t",mode="a",header=False)
else:
  tmrca_df.to_csv(tmrca_dataframe_file,index=False,sep="\t",mode="a",header=False)
