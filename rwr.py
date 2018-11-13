#!/usr/bin/env python3
import numpy as np
import pandas as pd

def rwr (restart_prob, tran_prob_matrix, length, tolerance, max_i):
 '''
 Random Walk with Restart function takes in:
  restart probabolity: the probability the walker jump to the start point
   restart matrix will be the identity matrix of restart probability multiplied by the restart probability  
   continue_prob = probability of not restarting and continue on the transition prob = 1-restart_prob
  transition probability matrix: row sums = 1; entry(i,j) = transition prob from j to i
  length = the number of nodes
  tolerance = the threashold at which the function considered to have converged
  max_i = maximum iterations allowed
 output:
  i: number of iterations
  residual: the final residual 
  the computed probability matrix
 '''
 continue_prob = 1 - restart_prob
 restart_matrix = np.identity(length) * restart_prob
 prob_matrix_old = np.full((length,length), 1/length)
 print (1/length)
 residual = 100
 i = 0 # iteration counter
 while (residual >= tolerance) and (i <= max_i):
  prob_matrix = np.matmul(prob_matrix_old, tran_prob_matrix) * continue_prob + restart_matrix
  residual = np.mean(np.absolute(prob_matrix - prob_matrix_old).sum(axis = 0)) #residule will be the mean of relative distance
  i += 1
  prob_matrix_old = prob_matrix.copy()
  print ('residual = ', residual, ' iteration = ', i)
 return i, residual, prob_matrix

def tpm(gg_affinity, gene_names, cc_affinity, cell_names, ppi_weight=0.5):
 '''
 Transition probabilily matrix construction, takes in:
  gene-gene affinity matrix is constructed from PPI with all genes appear in the query expression matrix. Gene-gene relationship not in the PPI should have been filled with 0s
  gene_names is the list containing the index of the gene gene affinity matrix
  cell-cell affinity matrix is constructed from the cell cell distance infered from the original expression matrix query.
  cell_names = the index of the cell cell affinity matrix
Both matrices should have been normalised so the rowsums and colsums equal to 1       
 Outputs: 
  node_names: the ordering of nodes in the form of an array of (cell_i, gene_j) tuples, permutation of each gene of each cells                        tran_prob_matrix: the constructed transition probability matrix  in the order of the node_names
  ppi_weight is defaulted to be 0.5, i.e., assigning gene-gene affinity and cell-cell affinity the same weight, it can also be set to a number between 0 and 1  
 '''
 len_c = len(cell_names)
 len_g = len(gene_names)
 length = len_c * len_g
 #initialise the transition probability to be 0 filled, except for the diagonal being the distance between the same nodes
 tran_prob_matrix = np.identity(length) * gg_affinity[1][1]
 #construct node_name list
 node_names=[]
 for c in range(len_c):
  cell_name = cell_names[c]
  for g in range(len_g):
   node_names.append((cell_name, gene_names[g])) 

 #######fill in gene-gene affinity into the transition probability matrix
 #iterating through the upper triangle of the gene-gene affinity matrix
 for i in range(len_g-1):
  for j in range (i+1, len_g):
   affinity = gg_affinity[i][j] * ppi_weight
   #fill in affinity for all the gene_i and gene_j pairs sharing the same cell, iterating from cell_0
   for x in range(0, length, len_g):
    tran_prob_matrix[x+i][x+j], tran_prob_matrix[x+j][x+i] = affinity, affinity

 ######fill in the cell-cell affinity #############
 for i in range(len_c-1):
  for j in range(i+1, len_c):
   affinity = cc_affinity[i][j]/2
   ii = i*len_g
   jj = j*len_g
   for x in range(len_g):
    tran_prob_matrix[ii+x][jj+x], tran_prob_matrix[jj+x][ii+x] = affinity, affinity
 return node_names, tran_prob_matrix

 
def c_c(exp_matrix)
 '''
 cell_cell affnity matrix construction function, measures cell cell euclidean distance then  normalise using 1/(1+d)
 input:
  expression matrix (cell by gene) in the format of pandas dataframe
 output:
  matrix: cell cell affinity matrix
  c_names: the index of cell names
  g_names: the index of gene names
 '''
 c_names = exp_matrix.index
 g_names = exp_matrix.column
 len_c = len(c_names)
 matrix = np.zeros([len_c,len_c])
 for i in range(len_c-1):
  for j in range(i+1, len_c):
   #calculation of the euclidean distance
   d = np.sum((exp_matrix[i]-exp_matrix[j])**2)
   d = d**(1/2)
   matrix[i][j], matrix[j][i] = d, d
 matrix = 1/(matrix+1)
 ############ matrices inputs are handled in the format of pandas dataframes #########################
 ppi = proper ppi matrix (only contains the genes in the expression matrix)
 expression_matrix is required to be gene by cell (should be handled in io step)
 outputs:
 length = number of cells * number of genes = the number of nodes 
 node_names: the ordering of nodes in the form of an array of (celli, genej) tuples, permutation of each gene of each cells
 tran_prob_matrix: the constructed transition probability matrix  in the ordering of the node_names
 '''
 

def main():

 print ('matrix a\n',a)
 print ('matrix b\n', b)
 result = tpm(a, name_a, b, name_b)
 print (result)
 print ('row sums', np.sum(result[1], axis=1))
if __name__ == '__main__':
 main()
