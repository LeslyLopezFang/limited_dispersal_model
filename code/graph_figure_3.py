import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import limited_dispersal_model_functions as db
import numpy as np

#Define the outfile name
outfile_dir="../figures/"

strong_alpha=1.08
weak_alpha=1.04
wright_fisher_alpha=1
N=100

#dpi specification for graphing
dpi_size=550

#Get position of all starting individuals
position_range = range(1,N+1)
#Get dij matrix: distance between individuals i and individuals j
dij = db.get_dij_matrix(N=N,pos=position_range)
#Get probability matrix for an individual in each position to get a parent based on position
probability_mass_function_strong = db.probability_picking_parent(strong_alpha,dij,N)
probability_mass_function_weak = db.probability_picking_parent(weak_alpha,dij,N)
probability_mass_function_wright_fisher = db.probability_picking_parent(wright_fisher_alpha,dij,N)


#Specify color palette
weak_model_color="#4a7398ff"
wright_fisher_model_color="#d94a3fff"
strong_model_color="#519043ff"

pmf_colors=(wright_fisher_model_color,weak_model_color
                   ,strong_model_color)


#Plot Figure 3
prob_outfile="{}pij_pmf_center_wf-weak-strong.jpeg".format(outfile_dir)
print(prob_outfile)

legend_labels=["Well Mixed","Weak limited dispersal","Strong limited dispersal"]

left_edge_index=2
right_edge_index=96
center_index=49

yaxis_step=0.01
yaxis_max=np.round(np.max(probability_mass_function_strong)+yaxis_step,2)
yticks=np.arange(0,yaxis_max,yaxis_step)
x_axis_range=np.arange(1,N+1)

probability_mass_function_linewidth=3

sns.set_style("white")
fig,ax=plt.subplots(figsize=(12,5))

#Probability mass function of Strong limited dispersal
sns.lineplot(y=probability_mass_function_strong[center_index]
             ,x=x_axis_range,linewidth=probability_mass_function_linewidth
             ,color=strong_model_color)

#Probability mass function of Weak limited dispersal
sns.lineplot(y=probability_mass_function_weak[center_index]
             ,x=x_axis_range,linewidth=probability_mass_function_linewidth
             ,color=weak_model_color)

#Probability mass function of Wright Fisher - 1/N
fig=sns.lineplot(y=probability_mass_function_wright_fisher[center_index]
             ,x=x_axis_range,linewidth=probability_mass_function_linewidth
             ,color=wright_fisher_model_color)

wf_line=":" #"dashdot"
ax.lines[2].set_linestyle(wf_line)

# '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
xticks=np.arange(0,110,10)
ax.set_xticks(xticks)
ax.set_yticks(yticks)
pmf_fontsize=26;pmf_axes_fontsize=20
pmf_x_label="Position of Individual"
pmf_y_label=r'$P(i=50,p;\alpha)$'
ax.set_ylabel(pmf_y_label,fontsize=pmf_fontsize)
ax.set_xlabel(pmf_x_label,fontsize=pmf_fontsize)
plt.tick_params(axis='both', which='major', labelsize=pmf_axes_fontsize)

pmf_legend_fontsize=18
legend_coordinate=1.03
patches = [ mpatches.Patch(color=pmf_colors[i]
            , label="{:s}".format(legend_labels[i]) ) for i in range(len(legend_labels)) ]
plt.legend(handles=patches,loc="upper right"
                   ,frameon=False,fontsize=pmf_legend_fontsize)

plt.tight_layout()
plt.savefig(prob_outfile,dpi=dpi_size,bbox_inches='tight')
plt.close()

