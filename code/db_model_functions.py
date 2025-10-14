import numpy as np
import scipy.stats
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

def binomial_distribution(n,p):
  #Binomial distribution can be used without modification. Get a random number from bi-dist
  return(scipy.stats.binom.rvs(n, p, size=1)[0])

def choose_parent_location_binomial(N,offspring_location):
  n = N - 1
  p = offspring_location/float(n)
  return(scipy.stats.binom.rvs(n,p,size=1)[0])

def choose_parent_location_poisson(N,offspring_location):
  #If offspring's location is smaller than the middle location
  if(offspring_location < (N/2)):
    #use regular poisson distribution
    location = scipy.stats.poisson.rvs(offspring_location, size=1)[0]
    while(location > N - 1):
      location = scipy.stats.poisson.rvs(offspring_location, size=1)[0]
  else:
    #If sample size (N) is even
    if(N%2==0):
      #Use the left side of the poisson distribution
      #To reflect: N - location of parent from poisson distribution
      adjusted_mu = N - offspring_location - 1
      location = N - scipy.stats.poisson.rvs(adjusted_mu,size=1)[0]
    #If sample size (N) is odd
    else:
      #Use the left side of the poisson distribution
      #To reflect: N - location of parent from poisson distribution
      adjusted_mu = N - offspring_location - 1
      location = N - scipy.stats.poisson.rvs(adjusted_mu,size=1)[0] - 1
  return(location)

def probability_picking_parent(alpha,dij,N,model_type):
  #Check that probability add to zero
  error = 1e-10
  assert(np.shape(dij)[0]==N)
  pij = np.zeros((N,N))
  for i in range(N):
    denominator = 0.0
    #Sum of all distances from individual i to all other individuals (j)
    dsum = sum(dij[i])
    for j in range(N):
      if(model_type=="linear-db"):
        pij[i][j] = linear_prob_numerator(dsum,dij[i][j],alpha)
      elif(model_type=="exponential-db"):
        pij[i][j] = exponential_prob_numerator(dij[i][j],alpha)
#      elif(model_type=="binomial"):
#        pij[i][j] = 
#      elif(model_type=="exponential"):
#        pij[i][j] = 
      else:
        print("specify probability type 'linear-db','exponential-db','binomial' or 'exponential'")
        sys.exit()
      denominator += pij[i][j]
    pij[i] /= denominator
    assert((1.0 - sum(pij[i]))<error)
  return pij

def choose_parent_ind(possible_parents):
  parent = rd.sample(possible_parents,1)[0]
  return(parent)

