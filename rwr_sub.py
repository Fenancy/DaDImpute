#!/usr/bin/env python3
import time
import numpy as np
import pandas as pd
from scipy import sparse
import ppi 


def imputation(prob_matrix, exp_matrix):
 imputed = prob_matrix @ exp_matrix
 return imputed

def rwr(adj_matrix, length, restart_prob = 0.5, tolerance=0.0001, max_i=100):
 '''                      
 Random Walk with Restart function takes in:                                  
  restart probabolity: the probability the walker jump back to the start point
   restart matrix will be the identity matrix of restart probability multiplied  by the restart probability
   continue_prob = probability of not restarting and continue on the transition prob = 1-restart_prob
  adj_matrix: row sums = 1; entry(i,j) = transition prob from i to j            
  length = the number of nodes                                                  
  tolerance = the threashold at which the function considered to have converged
  max_i = maximum iterations allowed                                            
 output:                                                                        
  i: number of iterations                                                       
  residual: the final residual of the computed probability matrix               
 '''
 continue_prob = 1 - restart_prob
 restart_matrix = np.identity(length) * restart_prob
 prob_matrix_old = np.full((length,length), 1/length)
 residual = 100
 i = 0 # iteration counter                         
 stime = time.time()
 while (residual >= tolerance) and (i <= max_i):
  prob_matrix = adj_matrix @ prob_matrix_old * continue_prob + restart_matrix
  residual = np.mean(np.absolute(prob_matrix - prob_matrix_old).sum(axis = 0))
  #residule will be the mean of relative distance                            
  i += 1
  prob_matrix_old = np.copy(prob_matrix)
  print ('residual = ', residual, ' iteration = ', i)
 etime = time.time()
 print(etime-stime)
 return prob_matrix

def rowNorm(matrix, length):
 for i in range(length):
  row_sum = np.sum(matrix[i])
  if (row_sum==0): continue
  matrix[i] = matrix[i]/row_sum
 return matrix

def extractExpMatrix(exp, geneBycell=True):
 '''
 expression matrix ( gene by cell in defaut) in the format of pandas dataframe
 output:
  matrix: cell cell affinity matrix
  c_names: the index of cell names
  g_names: the index of gene names
 '''
 if (geneBycell==False):
  exp = exp.transpose()
 g_names = exp.index
 c_names = exp.columns
 exp_matrix = exp.values
 len_c = len(c_names)
 len_g = len(g_names)
 return (exp_matrix, c_names, g_names, len_c, len_g)

def gen_cc_matrix(exp_matrix, len_c):
 '''
 cell_cell affnity matrix construction function, measures cell cell euclidean distance
 ( cell by gene in defaut)
 '''
 matrix = np.zeros([len_c,len_c])
 for i in range(len_c-1):
  for j in range(i+1, len_c):
   #calculation of the euclidean distance
   d = np.sum((exp_matrix[i]-exp_matrix[j])**2)
   d = 1/(1+d**(1/2))
   matrix[i][j], matrix[j][i] = d, d
  #normalisation
  matrix = rowNorm(matrix, len_c)
 return matrix

def gen_g2i(g_names, len_g):
 gene2index = {g_names[i]:i for i in range(len_g)}
 return gene2index 

def main():
############# read expression matrix ############
 exp = pd.read_csv("/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/gliobalstoma_sub_missing20.csv", header=0, index_col=0)
#generate the cell-cell affinity matrix
 (exp_matrix, c_names, g_names, len_c, len_g) = extractExpMatrix(exp)
#read in the porper gene-gene affinity matrix (sparse)
 #ppi_file= "/users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/human.ppi.txt"
 #mapping = "/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/mart_gene_to_protein_mapping.txt"
 #(gg_matrix, len_ppi) = ppi.gg_from_ppi(ppi_file, mapping, g_names, len_g)
 #np.save("/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/g_matrix.npy", gg_matrix)
 gg_matrix = np.load("/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/g_matrix.npy")
 len_ppi = len(gg_matrix)
 gg_rwr = rwr(gg_matrix, len_ppi, 0.5)
 gg_rwr = gg_rwr[:len_g, :len_g]
 gg_rwr = rowNorm(gg_rwr, len_g)
 np.save("/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/ggrwr.npy", gg_rwr)
 imputedgg = imputation(gg_rwr, exp_matrix)
 np.save("/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/imputedgg.npy", imputedgg)
 imputedgg = imputedgg.transpose()
 cc_matrix = gen_cc_matrix(imputedgg, len_c)
 cc_rwr = rwr(cc_matrix, len_c, 0.5)
 cc_rwr = rowNorm(cc_rwr, len_c)
 np.save("/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/ccrwr2.npy", cc_rwr)
 
 imputedcc = imputation(cc_rwr, imputedgg)
 np.save("/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/imputedcc.npy", imputedcc)
 #improve the counts so that the libsize is the same as input
 #for i in range(len_c):
  #col_sum = np.sum(imputed2[:, i])
  #if (col_sum==0): continue
  #imputed2[:,i] = libsize * imputed2[:,i] /col_sum
 
 #np.save("/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/imputed_and_normed.npy", imputed2)
# imputed2 = pd.DataFrame(data=imputed1, index=g_names, columns=c_names)
# imputed2.to_pickle("/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/imputed_and_normed.pkl")
if __name__ == '__main__':
 main()
