#!/usr/bin/env python3
import numpy as np
import pandas as pd
from scipy import sparse

def gen_g2i(exp_genes, ppi_genes, len_g):
 '''
 generates a dictionary mapping gene names to its index in the later gene-gene matrix:
 genes that are in the input expression matrix
 followed by the extra genes that are only in the ppi network
 which makes it easier to get rid of the unwanted/extra genes after rwr
 '''
 set_exp = set(exp_genes)
 set_ppi = set(ppi_genes)
 extra = set_ppi - set_exp #compare and get the non-overlaped genes in the expression matrix and ppi network
 extra = list(extra)
 len_extra = len(extra)
 g2i_adj = {exp_genes[i]:i for i in range(len_g)} #first map the genes in the input file
 g2i_adj.update({extra[i]:i+len_g for i in range(len_extra)}) #add the extra genes that are only in the ppi network
 return g2i_adj, len(set_ppi)

def gg_from_ppi(ppi_file, mapping_file, g_names, len_g, sep = ' '):
 '''                                                            
 Generates the gene-gene matrix from ppi network
  '''
 #read the ppi file to data frame 
 ppi = pd.read_csv(ppi_file, sep=sep)
 #get the column names
 protein1 = ppi.columns[0]
 protein2 = ppi.columns[1]
 #get the weightss
 weights = ppi[ppi.columns[9]] #column 9 is the experiments scores
 #trim off the first 5 characters of the protein names ('9606.')
 p1_names = ppi[protein1].apply(lambda x: x[5:])
 p2_names = ppi[protein2].apply(lambda x: x[5:])
 #get the unique protein names
 protein_names = set(p1_names)
 #convert to list
 protein_names = sorted(list(protein_names))
 #get the map to genes
 (protein2gene, ppi_genes) = protein_to_gene(protein_names, mapping_file)
 #convert the weight of those proteins that are not in the map to be 0
 #weights[not_in_map] = 0
 g2i_adj, len_m = gen_g2i(g_names, ppi_genes, len_g)
 #################construction the protein-protein affinity matrix##################
 #get the indices (row number) of nonzero edges
 ind_non0edge = weights.nonzero()
 ind_non0edge = ind_non0edge[0]
 #the indices of the p1 and p2 in the sorted unique protein list, to be the indices of the row/col
 ind_row, ind_col = [], []
 #get the protein ids corresponding to each edge
 #np.save("./p2g.npy",protein2gene)
 for ind in ind_non0edge:
  ind_row.append(g2i_adj[protein2gene[p1_names[ind]]])
  ind_col.append(g2i_adj[protein2gene[p2_names[ind]]])

 #construct the sparse matrix
 matrix = sparse.csr_matrix((weights[ind_non0edge].values, (ind_row, ind_col)), shape = (len_m, len_m), dtype=float)
 #normalisation
 for i in range(len_m):
  row_sum = np.sum(matrix[i])
  if (row_sum==0): continue
  matrix[i] = matrix[i]/row_sum
 return matrix, len_m

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
 len_map = len(map[col_p])
 #############
 p2g = {map[col_p].iloc[i]:map[col_g].iloc[i] for i in range(len_map)}
 gene_names = []
 for p in protein_names:
  try:
   gene_names.append(p2g[p])
  except KeyError:
   gene_names.append(p)
   continue
 protein2gene = {protein_names[i]:gene_names[i] for i in range(len(protein_names))}
# gene_sorted = sorted(list(set(gene_names)))
# len_matrix = len(gene_sorted)
 return protein2gene, gene_names

