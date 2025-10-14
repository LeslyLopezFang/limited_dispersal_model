import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description="Graph a 3 panel PCA of simulated SNPs for the following populations: well mixed, weak limited dispersal, strong limited dispersal.")

parser = argparse.ArgumentParser(description="Graph violinplot of average tij per distance between individuals i and j estimated from simulated SNPs.")
parser.add_argument('-N',action="store",dest="N",type=int,help="Constant population size.")
parser.add_argument('-L',action="store",dest="L",type=int,help="L is the number of loci simulated by one run of the limited dispersal model.")
parser.add_argument('-s',action="store",dest="total_sample_size_of_tij",type=int,help="What is the total sample size of all tij's (ie (N*(N-1))/2 * sample size of one tij?")
param = parser.parse_args()

#Files
results_path = "../results/"
figures_path = "../figures/"

pc_file="{}principal_components_db-model_N-{}_L-{}_runs-{}.csv".format(
             results_path,param.N,param.L,param.total_sample_size_of_tij)
figure_file="{}scatterplot_pca_db-model_N-{}_L-{}_runs-{}".format(
             figures_path,param.N,param.L,param.total_sample_size_of_tij)

pc_df=pd.read_csv(pc_file,sep="\t")

#Graphing parameters
dpi_size=550

dij_palette="Spectral_r"

pca_fontsize=38

fontsize=30
axes_fontsize=40
font = {'family': 'arial',
        'color':  'black',
        'weight': 'bold',
        'size': fontsize
        }

sns.set_style("white")
position_palette="hls"

fig,ax=plt.subplots(ncols=3,figsize=(25,7))
size=100

model_type_label="Well Mixed"
fig=sns.scatterplot(x="PC1",y="PC2",hue="distance_label"
               ,data=pc_df[pc_df.model_type==model_type_label],legend=False
                    ,palette=position_palette,s=size,ax=ax[0])

model_type_label="Weak limited dispersal"
fig=sns.scatterplot(x="PC1",y="PC2",hue="distance_label"
               ,data=pc_df[pc_df.model_type==model_type_label],legend=False
                    ,palette=position_palette,s=size,ax=ax[1])

model_type_label="Strong limited dispersal"
fig=sns.scatterplot(x="PC1",y="PC2",hue="distance_label"
               ,data=pc_df[pc_df.model_type==model_type_label],legend=True
                    ,palette=position_palette,s=size,ax=ax[2])

legend_coordinate=1
fig.legend(loc="center right",bbox_to_anchor=(legend_coordinate, 0.25, 0.5, 0.5)
          ,title=r'$loc(I_{i,t})$',frameon=False
           ,title_fontsize=fontsize,fontsize=fontsize,borderaxespad=0)#,fontsize='large')

letter_coordinate=-0.3#-0.185
ax[0].text(letter_coordinate,1,"A", fontdict=font,transform=ax[0].transAxes)
ax[1].text(letter_coordinate,1,"B", fontdict=font,transform=ax[1].transAxes)
ax[2].text(letter_coordinate,1,"C", fontdict=font,transform=ax[2].transAxes)

for subax in ax:
    subax.tick_params(axis='both', which='both', labelsize=fontsize)
    subax.set_xlabel(subax.get_xlabel(), fontsize=pca_fontsize)
    subax.set_ylabel(subax.get_ylabel(), fontsize=pca_fontsize)

plt.tight_layout()
plt.savefig("{}.jpeg".format(figure_file),dpi=dpi_size,bbox_inches='tight')
plt.savefig("{}.svg".format(figure_file),dpi=dpi_size,bbox_inches='tight')
plt.close()
