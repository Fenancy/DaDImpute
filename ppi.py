#!/usr/bin/env python3                                                                                              
import numpy as np
import pandas as pd

def gg_csv(ppi_file, sep = ' '):
 '''                                                            
 Generates the gene-gene matrix from ppi network
  '''
 #read the ppi file to data frame 
 ppi = pd.read_csv(ppi_file, sep=sep)
 #get the column names
 protein1 = ppi.columns[0]
 protein2 = ppi.columns[1]
 #get the weights
 weight = ppi[ppi.columns[9]] #column 9 is the experiments scores
 #trim off the first 5 characters of the protein names ('9606.')
 p1_names = ppi[protein1].apply(lambda x: x[5:])
 p2_names = ppi[protein2].apply(lambda x: x[5:])
 #get the unique protein names
 set1 = set(p1_names)
 set2 = set(p2_names)
 if (set1==set2):
  protein_names = set1
 #if set 1 and set2 do not completely overlap, take the union of them
 else: protein_names = set1 | set1
 #convert to list
 protein_names = list(protein_names)
 num_protein = len(protein_names)
 #sort the names                                                                
 protein_names = sorted(protein_names)
 #construc the protein-protein affinity matrix
 matrix = np.zeros([num_protein, num_protein])
 for i in range(len(weight)):
  edge_weight = weight[i]
  #if the edge weight is zero, skip
  if edge_weight==0: continue
  #else, find the indices of protein1 and protein2 in this row
  index_row = protein_names.index(p1_names[i])
  index_col = protein_names.index(p2_names[i])
  #fill the matrix
  matrix[index_row][index_col], matrix[index_col][index_row] = edge_weight, edge_weight
 #normalisation
 for i in range(num_protein):
  row_sum = np.sum(matrix[i])
  matrix[i] = matrix[i]/row_sum
 print (matrix)
 print (protein_names) 
 return matrix, num_protein, protein_names

def protein_to_gene(protein_names, mapping, sep='\t'):
 '''
 maps gene names to protein names
 input:
  protein_names = protein names index array
  mapping = protein to gene ensemble name mapping file in csv    
 output: 
  gene_names = mapped gene names index array
 '''
 map = pd.read_csv(mapping, sep=sep)
 #get the column name that contains the protein ids
 col_p = map.columns[1] 
 #the column name that contains the gene ids
 col_g = map.columns[0]
 #removes the rows where the gene does not produce protein
 map = map.dropna()
 #resort the series by protein id
 map = map.sort_values(by=[col_p])
 i = 0
 len_map = len(map[col_p])
 gene_names = []
 for p in protein_names:
  for x in range(i, len_map):
   #iteration through the protein id list (from the last found id)
   if (map[col_p][x] == p):
    i=x
    gene_names.append(map[col_p][x])
 return gene_names


def main():
 ppi_file="/users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/human.ppi.txt"
 (gg_matrix, num_protein, protein_names) = gg_csv(ppi_file)
 mapping = "/Users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/mart_gene_to_protein_mapping.txt"
 gene_names = protein_to_gene(protein_names, mapping)
 #seve the matrix as dataframe with column names
 gg_matrix = pd.DataFrame(data = matrix,  columns=protein_names)
 gg_matrix.to_pickle("/users/Vie/Google Drive/Autumn 2018/COMP 401/scRNA seq imputation/gene_gene_affinity_matrix.pkl")

if __name__ == '__main__':
 main()
