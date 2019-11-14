import numpy as np
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

def linear_prob_numerator(dsum,dij_number,alpha):
  return((dsum-dij_number)**alpha)

def exponential_prob_numerator(dij_number,alpha):
  return(alpha**(-dij_number))

def probability_picking_parent(alpha,dij,N,type):
  #Check that probability add to zero
  error = 1e-10
  assert(np.shape(dij)[0]==N)
  pij = np.zeros((N,N))
  for i in range(N):
    denominator = 0.0
    #Sum of all distances from individual i to all other individuals (j)
    dsum = sum(dij[i])
    for j in range(N):
      if(type=="linear"):
        pij[i][j] = linear_prob_numerator(dsum,dij[i][j],alpha)
      elif(type=="exponential"):
        pij[i][j] = exponential_prob_numerator(dij[i][j],alpha)
      else:
        print("specify probability type 'linear' or 'exponential'")
        sys.exit()
      denominator += pij[i][j]
    pij[i] /= denominator
    assert((1.0 - sum(pij[i]))<error)
  return pij

def choose_parent_ind(possible_parents):
  parent = rd.sample(possible_parents,1)[0]
  return(parent)

