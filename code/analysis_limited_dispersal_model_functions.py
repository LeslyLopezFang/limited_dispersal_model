import numpy as np
import dendropy
import random

def get_allelic_matrix(newick_file,N):
  """Compute an LxN numpy array of biallelic SNPs where the prescence of the ancestral allele at a loci
     is denoted by 0 and the prescence of a mutation (i.e. the derived allele) is denoted by a 1.
     Input: the newick file that contains a gene tree for each loci of all the N individuals. The file
     has L number of gene trees.
     Output: numpy array of biallelic alleles and a list of the total branch lengths of each tree. 
  """
  tree_list=dendropy.TreeList()
  tree_list.read(path=newick_file,schema="newick")
  #X matrix
  mut_list = []
  #Total tree length
  total_branch_lengths=[]
  for tree in tree_list:
    tl = tree.length()
    total_branch_lengths.append(tl)
    rn = random.uniform(0,1)
    breakpt = tl * rn
    snp_geno = [0]*N
    tree_length_traversed = 0
    for branch in tree.postorder_edge_iter():
      tree_length_traversed += branch.length
      if tree_length_traversed >= breakpt:
        for leaf in branch._head_node.leaf_iter():
          snp_geno[int(leaf.taxon.label)] = 1
        mut_list.append(snp_geno)
        break
  return(np.asarray(mut_list),total_branch_lengths)

def create_distance_label(x,step=10):
  """Create a distance label to group distances by step (default step=10). A distance of 1 gets label
     string "0-10".
     Input: x: distance between two individuals (int or float) and step: bin for distance label (int)
            Default step=10
     Output: string of distance label in format of "int-int" with x within range of both integers
  """
  return("{}-{}".format(int(x-(x%step)),int((x-(x%step))+step)))

def create_model_type_label():
  """Create a label based on the strength of the parameter alpha from the limited dispersal model.
    Output: dictionary with keys of float type alpha and values of string such that
            an alpha of 1 is well mixed(Wright Fisher), 1.08 is strong and 1.04 is weak
  """
  model_type_dict={1.0:"Well Mixed",1.04:"Weak limited dispersal",1.08:"Strong limited dispersal"}
  return(model_type_dict)

def compute_estimated_average_coalescent_times(M_matrix):
  """Estimate average coalescent times tij for all pairs of N individuals using M matrix from simulated data (SNPs).
     Solve equation LESLY linearly
     Input: NxN M matrix from simulated SNPs data
     Output:Numpy array of N*(N-1)/2 solutions (estimated average coalescent times of each pair of individuals) 
  """
  all_equations=[]
  solutions=[]
  N = int(M_matrix.shape[0])
  #Get variables and solutions based off equation LESLY
  for i in range(N):
    for j in range(i+1,N):
      variables=list()
      for k in range(N):
        for l in range(k+1,N):
          if (i in (k,l)) & (j in (k,l)):
            variables.append(((-np.square(N)+(2*N)-2)/np.square(N)))
          elif (i in (k,l)) or (j in (k,l)):
            variables.append(((N-2)/np.square(N)))
          else:
            variables.append((-2/np.square(N)))
      all_equations.append(variables)
      solutions.append(M_matrix[i][j])

  all_equations=np.asarray(all_equations)
  solutions=np.asarray(solutions)
  #Note: solutions.shape[0] == all_equations.shape[0] == all_equations.shape[1] == N*(N-1)/2
  #Solve linear equations
  estim_t=np.linalg.inv(all_equations).dot(solutions)
  return(estim_t)

def compute_principal_components_from_m_matrix(M_matrix):
  """Get PC1 and PC2 from eigenvectors of M marix from SNPs.
     Input: NxN M matrix from simulated SNPs data
     Output: numpy array of PC1, numpy array of PC2
  """
  eigenvalue,eigenvector = np.linalg.eig(M_matrix)
  pc1_index=np.where(eigenvalue==np.max(eigenvalue))[0][0]
  pc2_index=np.where(eigenvalue==np.max(eigenvalue[eigenvalue!=np.max(eigenvalue)]))[0][0]
  pc1=np.sqrt(eigenvalue[pc1_index])*eigenvector[:,pc1_index]
  pc2=np.sqrt(eigenvalue[pc2_index])*eigenvector[:,pc2_index]
  return(pc1,pc2)


def compute_average_coalescence_times_matrix(G_sum_file,total_sample_size_of_tij):
  """Average the matrix of the sum of coalescent times from limited dispersal model
     Input: String of path to G_sum and the total sample size of all the coalescent times
     Output: Numpy matrix containing the average coalescent time of each pair of N individuals from limited dispersal model
  """
  G_sum = np.load(G_sum_file)
  N=G_sum.shape[0]
  sample_size_of_tij=np.floor(total_sample_size_of_tij/((N*(N-1))/2))
  G_average = (1/sample_size_of_tij)*G_sum
  return(G_average)

