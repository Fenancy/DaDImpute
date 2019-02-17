#!/usr/bin/env python3
import time
import numpy as np
import pandas as pd
from scipy import sparse
import ppi 
import argparse


def argPaser():
 parser = argparse.ArgumentParser()
 parser.add_argument('input_matrix', type=str, help='name of the input expression file (csv)', default='./glioblastoma_sub_missing20.csv') #the glioblastoma expression matrix is pre-trimmed to contain only the genes that can be mapped to the proteins present in the ppi network, it's then to be masked (20% in this file)
 parser.add_argument('-ppi', '--ppi_file', type=str, help='name of the ppi file', default='./human.ppi.txt')
 parser.add_argument('-m', '--mapping_file', type=str, help='names of the protein to gene mapping file', default='./mart_gene_to_protein_mapping.txt')
 parser.add_argument('-p', '--prob_restart', type=float, default = 0.5, help='the restart probability to be employed in the random walk with restart process, defaulted to be 0.5')
 parser.add_argument('-t', '--tolerance', type=float, default = 0.0001, help='the tolerance/threashold to resolve the convergence of the matrix in the random walk with restart process, defaulted to be 0.0001')
 parser.add_argument('-i', '--max_iteratiion', type=int, default = 100, help='the maximum iteration to be performed in the random walk with restart process, defaulted to be 100')
 parser.add_argument('--path_gg', type=str, help='file path to store the gene-gene affinity/mimic netSmooth imputed matrix in .npy', default='./ggImpute.npy')
 parser.add_argument('--path_gc', type=str, help='file path to store the then cell-cell affinity/mimic MAGIC imputed matrix in .npy', default='./ccImpute.npy')

 
 args = parser.parse_args()
 return args

def imputation(prob_matrix, exp_matrix):
 #apply the result probability matrix to the input matrix
 imputed = prob_matrix @ exp_matrix
 return imputed

def rwr(adj_matrix, length, restart_prob, tolerance, max_i):
 '''                      
 Random Walk with Restart function, takes in:                                  
  restart probabolity: the probability the walker jump back to the start point
   continue_prob = probability of not restarting and continue on the transition prob = 1-restart_prob
  adj_matrix: row sums = 1; entry(i,j) = transition prob from i to j            
  length = the number of nodes                                                  
  tolerance = the threashold at which the function considered to have converged
  max_i = maximum iterations allowed                                          
 Using equation:
  prob_matrix' = adjacency matrix * prob_matrix *(1-restart probability) + restart probability * identity matrix
               = adjacency matrix * prob_matrix *   continue_prob        +             restart_matrix
 output:                                                                        
  i: number of iterations                                                       
  residual: the final residual of the computed probability matrix               
 '''
 continue_prob = 1 - restart_prob
 restart_matrix = np.identity(length) * restart_prob
 prob_matrix_old = np.full((length,length), 1/length) #to start the diffusion from a probability matrix with equal weight in each entry
 residual = 100
 i = 0 # iteration counter                         
 #stime = time.time() # for referencing computing time
 while (residual >= tolerance) and (i <= max_i):
  prob_matrix = adj_matrix @ prob_matrix_old * continue_prob + restart_matrix
  residual = np.mean(np.absolute(prob_matrix - prob_matrix_old).sum(axis = 0))
  #residule will be the mean of relative distance                            
  i += 1
  prob_matrix_old = np.copy(prob_matrix)
 #print ('residual = ', residual, ' iteration = ', i)
 #etime = time.time()
 #print(etime-stime)
 return prob_matrix

def rowNorm(matrix, length):
 #function that row normalises a matrix, takes in length = number of rows in the matrix 
 for i in range(length):
  row_sum = np.sum(matrix[i])
  if (row_sum==0): continue
  matrix[i] = matrix[i]/row_sum
 return matrix

def extractExpMatrix(exp, geneBycell=True):
 '''
 helper method to extract information from the input DataFrame matrix and convert it to be gene by cell 
 output:
  exp_matrix: the expression matrix (gene by cell) in the format of numpy array
  c_names: the index of cell names
  g_names: the index of gene names
  len_c: number of cells
  len_g: number of genes
 '''
 if (geneBycell==False):
  exp = exp.transpose() # if needed, convert the matrix to be gene by cell first
 g_names = exp.index
 c_names = exp.columns
 exp_matrix = exp.values
 len_c = len(c_names)
 len_g = len(g_names)
 return (exp_matrix, c_names, g_names, len_c, len_g)

def gen_cc_matrix(exp_matrix, len_c):
 '''
 Generating cell-cell matrix: cell_cell affnity matrix construction function, measures the Euclidean distance between each cell cell pair
 '''
 matrix = np.zeros([len_c,len_c]) #initialise a 0 filled matrix
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
 #Function generating the dictionary mapping gene name to its index in the expresssion matrix
 gene2index = {g_names[i]:i for i in range(len_g)}
 return gene2index 

def main():
 args = argPaser()
 exp = pd.read_csv(args.input_matrix, header=0, index_col=0) #read the input expression file to panda DataFrame
 (exp_matrix, c_names, g_names, len_c, len_g) = extractExpMatrix(exp) #convert to numpy array and extract cell/gene names as well as the dimension of the matrix
 #generation of a gene-gene affinity matrix from the ppi file and protein to gene mapping file, using method from python file 'ppi'
 (gg_matrix, len_ppi) = ppi.gg_from_ppi(args.ppi_file, args.mapping_file, g_names, len_g)
 #run random walk with restart method on the gene-gene matrix
 gg_rwr = rwr(gg_matrix, len_ppi, args.prob_restart, args.tolerance, args.max_iteratiion)
 gg_rwr = gg_rwr[:len_g, :len_g] #trim off the genes that are not in the input file (only in the ppi network)
 gg_rwr = rowNorm(gg_rwr, len_g) #normalisation
 #np.save("you can save the gene-gene probability matrix after rwr to avoid repetitive work",gg_rwr)
 #gg_rwr = np.load("the npy file containing the gene-gene matrix after rwr from previous computation")
 imputedgg = imputation(gg_rwr, exp_matrix) #first impute the expression matrix by 'mimic' netSmooth
 np.save(args.path_gg, imputedgg)
 imputedgg = imputedgg.transpose() #for 'mimic' MAGIC the input needs to be cell by gene
 cc_matrix = gen_cc_matrix(imputedgg, len_c) #generating the cell-cell affinity matric
 cc_rwr = rwr(cc_matrix, len_c, args.prob_restart, args.tolerance, args.max_iteratiion) #imputation
 cc_rwr = rowNorm(cc_rwr, len_c) #normalisation
 imputedgc = imputation(cc_rwr, imputedgg) #apply the cell-cell probability matrix after rwr to the previously imputed matrix
 np.save(args.path_gc, imputedgc)

if __name__ == '__main__':
 main()
