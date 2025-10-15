import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Graph violinplot of Tmrca.")

parser.add_argument('-N',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('-L',action="store",dest="L",type=int,help="L is the number of loci simulated by one run of the limited dispersal model.")
parser.add_argument('-s',action="store",dest="total_sample_size_of_tij",type=int,help="What is the total sample size of all tij's (ie (N*(N-1))/2 * sample size of one tij?")
param = parser.parse_args()

#Files
figures_path="../figures/"
tmrca_path="../data/tmrca/"
tmrca_dataframe_file="{}tmrca_dataframe_db-model_N-{}_L-{}_runs-{}.csv".format(tmrca_path
                  ,param.N,param.L,param.total_sample_size_of_tij)
tmrca_figure_file="{}violinplot_tmrca_db-model_N-{}_L-{}_runs-{}".format(figures_path
             ,param.N,param.L,param.total_sample_size_of_tij)

tmrca_df=pd.read_csv(tmrca_dataframe_file,sep="\t")

#graphing parameters
tmrca_fontsize=14
fontsize=30
dpi_size=550

fig,ax=plt.subplots(figsize=(10,5))

bplot = sns.violinplot(x='model_type', y='tmrca',
                data=tmrca_df, showmeans=True
            ,palette="Pastel1")
##analytical expectation of tmrca
ana_exp = 2.0*param.N*(1.-(1./param.N))
#try color red?
bplot.axhline(ana_exp,color="brown",linestyle="--")

#ax.set_ylim((0,np.max(tmrca.tmrca)+(50-(np.max(tmrca.tmrca)%50))))
#ax.set_yticks(np.arange(0,1800,150))
ax.set_xlabel("")#,fontsize=fontsize)
ax.set_ylabel(r'$T_{MRCA}$',fontsize=fontsize)
plt.tick_params(axis='both', which='major', labelsize=tmrca_fontsize)

plt.tight_layout()
plt.savefig("{}.jpeg".format(tmrca_figure_file),dpi=dpi_size,bbox_inches='tight')
plt.close()

