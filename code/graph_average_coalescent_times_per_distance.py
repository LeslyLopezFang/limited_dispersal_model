import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description="Graph violinplot of average tij per distance between individuals i and j.")

parser.add_argument('-N',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('-L',action="store",dest="L",type=int,help="L is the number of loci simulated by one run of the limited dispersal model.")
parser.add_argument('-s',action="store",dest="total_sample_size_of_tij",type=int,help="What is the total sample size of all tij's (ie (N*(N-1))/2 * sample size of one tij?")
param = parser.parse_args()

#Files
results_path = "../results/"
figures_path = "../figures/"

#File to with average tij per distance from limited dispersal model results
average_tij_file="{}average_coalescent_times_db-model_N-{}_L-{}_runs-{}.csv".format(results_path
                  ,param.N,param.L,param.total_sample_size_of_tij)
figure_file="{}violinplot_average_pairwise_coalescent_times_db-model_N-{}_L-{}_runs-{}".format(
                  figures_path,param.N,param.L,param.total_sample_size_of_tij)

df=pd.read_csv(filepath_or_buffer=average_tij_file,sep="\t")

#Graphing parameters
dpi_size=550
fontsize=30
axes_fontsize=40
text_font={'family': 'arial',
        #'color':  'black',
        'weight': 'normal',
        'size': axes_fontsize
        }


#Violing Plot of average tij per distance
plt.figure(figsize=(25,10))

grped_bplot = sns.violinplot(x='distance_label', y='tij_average'
                ,data=df, showmeans=True
                ,hue="model_type", palette="Pastel1")

grped_bplot.axhline(param.N,color="brown",linestyle='--')

grped_bplot.set_xlabel(r'Geographic Distance',font=text_font)
grped_bplot.set_ylabel(r'Pairwise coalescent time',font=text_font)

av_step=50
ytick_space_from_x_axis=10; min_order=10; max_order=100
ytick_start=round(min(df.tij_average)) - (round(min(df.tij_average)) % min_order)
ytick_stop=round(max(df.tij_average)) + (max_order - round(max(df.tij_average)) % max_order)
#LESLY
#df.tij_average.min().round() -> gives like 22.0
#y_ticks=np.arange(np.floor(np.min(df.tij_average))-10,np.max(df.tij_average)+av_step,av_step)
y_ticks=[ytick for ytick in range(ytick_start,ytick_stop+av_step,av_step)]
grped_bplot.set_yticks(y_ticks)

legend_coordinate=0.9
grped_bplot.legend(loc="lower right",frameon=False
                   ,title_fontsize=fontsize,fontsize=fontsize)


plt.tick_params(axis='both', which='both', labelsize=fontsize)

plt.tight_layout()
plt.savefig(fname="{}.svg".format(figure_file),dpi=dpi_size,bbox_inches='tight')
plt.savefig(fname="{}.jpeg".format(figure_file),dpi=dpi_size,bbox_inches='tight')
plt.close()
