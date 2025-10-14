import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

#Files
figures_path="../figures/"
tmrca_path="../data/tmrca/"
tmrca_dataframe_file="{}tmrca_dataframe_db-model_N-{}_L-{}_runs-{}.csv".format(tmrca_path
                  ,param.N,param.L,param.total_sample_size_of_tij)
tmrca_figure_file="{}violinplot_tmrca_db-model_N-{}_alpha-{}_L-{}_runs-{}.txt".format(figures_path,param.N
             ,param.alpha,param.L,param.total_sample_size_of_tij)

tmrca_df=pd.read_csv(tmrca_dataframe_file)

#graphing parameters
tmrca_fontsize=14

fig,ax=plt.subplots(figsize=(10,5))

bplot = sns.violinplot(x='model_type', y='tmrca',
                data=tmrca_df, showmeans=True
            ,palette="Pastel1")
#bplot.axhline(100,color="red",linestyle='--')
##analytical expectation of tmrca
ana_exp = 2.0*param.N*(1.-(1./param.N))
bplot_ax.hline(ana_exp,color="red",linestyle="--")

#ax.set_ylim((0,np.max(tmrca.tmrca)+(50-(np.max(tmrca.tmrca)%50))))
ax.set_yticks(np.arange(0,1800,150))
ax.set_xlabel("")#,fontsize=fontsize)
ax.set_ylabel(r'$T_{MRCA}$',fontsize=fontsize)
plt.tick_params(axis='both', which='major', labelsize=tmrca_fontsize)

plt.tight_layout()
plt.savefig("{}.jpeg".format(outfile),dpi=dpi_size,bbox_inches='tight')
plt.close()

