import numpy as np
import random
import sys

def get_dij_matrix(pos,N):
  """Get matrix for dij
  Input: N: integer, constant population size
  Output: dij: numpy matrix with distance from individual i to individual j
  """
  #Initiate NxN matrix of 0's
  dij = np.zeros((N,N))
  for i in range(N):
    for j in range(i,N):
      dij[i][j] = abs(pos[i]-pos[j])
  return(dij + np.transpose(dij))

def pmf_numerator(dij_number,alpha):
  return(alpha**(-dij_number))

def probability_picking_parent(alpha,dij,N,error=1e-6):
  #Check that the distance matrix takes into account all individuals
  assert(np.shape(dij)[0]==N)
  pij = np.zeros((N,N))
  for i in range(N):
    denominator = 0.0
    #Sum of all distances from individual i to all other individuals (j)
    dsum = sum(dij[i])
    for j in range(N):
      pij[i][j] = pmf_numerator(dij[i][j],alpha)
      denominator += pij[i][j]
    pij[i] /= denominator
    #Probability must add to 1.0 within margin of error
    try:
      assert(abs(1.0 - sum(pij[i]))<error)
    except AssertionError:
      print("Probability not adding to 1.0 within {} error for alpha = {} & N = {}".format(error,alpha,N))
      sys.exit(1)
  return pij

def probability_picking_parent_two_populations(alpha,dij,N,barrier_size,error):
  """Get the probability of picking a parent based on location for an individual in a subpopulation. There are two
  subpopulations of size N with a barrier of length barrier_size. The probability of picking a parent in the location within
  the barrier is 0. Probability is then re-normalized.
  """
  #Probability without a barrier for one population
  pij = probability_picking_parent(alpha=alpha,dij=dij,N=(2*N)+barrier_size,error=1e-6)
  #Want to set pij(location=barrier)=0 and re-normalize using log
  denom=np.log(np.sum(pij[:,:N],axis=1)+np.sum(pij[:,N+barrier_size:],axis=1))
  pij_log=np.log(pij)
  pij_transformed=pij_log-denom
  pij_transformed=np.exp(pij_transformed)
  pij_transformed[:,N:N+barrier_size]=0
  print(np.min(np.sum(pij_transformed,axis=1)),np.max(np.sum(pij_transformed,axis=1)))
#  assert False not in (abs(1-np.sum(pij_transformed,axis=1))<error),"Pij transformed not adding to 1 within 1e-06"
  return pij_transformed



def choose_parent_ind(possible_parents):
  parent = rd.sample(possible_parents,1)[0]
  return(parent)

def branch_length(new_node,internal_node,node_time):
  """Compute the branch length between two nodes.
  """
  return(float(node_time[internal_node]-node_time[new_node]))

def string_node(children_list,inner_node,string_dicc,node_time,scalar_branch_length,triple_coalescence_branch_length=0.0,):
  """Write the children nodes of an inner node in newick file format
     '(child_one:branch length, child_two:branch length, child_three:branch length)'
  """
  #Get the number of times the multiple coalescence branch length needs to added on
  mew=(np.arange(1,(len(children_list))+1) - 2).clip(min=0)
  s="("
  #Go through each child of the inner_node
  for i in range(len(children_list)):
    child = children_list[i]
    #Does the inner node have more than 2 children (like a triple coalescence)
    if(i > 1):
#      branch_length=triple_coalescence_branch_length*mew[i]+((node_time[inner_node]-node_time[child])
#                        *scalar_branch_length)
      branch_length=((node_time[inner_node]-node_time[child])*scalar_branch_length)
#    branch_length_dict[child]=branch_length
      s="("+s[:-1]+"):"+str(triple_coalescence_branch_length)+","
    #If the child node is not going to have multiple coalescence branch length
    else:
      branch_length=float(node_time[inner_node]-node_time[child])*float(scalar_branch_length)
#    branch_length_dict[child]=branch_length
    #Does the child node have children associated with it
    if(string_dicc[child]):
      s+= string_dicc[child]
    #or is it a tip?
    else:
      s+=str(child)
    s+=":"+str(branch_length)+","
  #Exclude last ',' and add closing paranthesis
  return(s[:-1]+")")

