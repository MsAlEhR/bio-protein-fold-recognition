# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 11:41:08 2018

@author: Saleh, Mir, A.
"""
import pandas as pd

def bigram(pssm_matrix):

    """
    Get features by bigram method

    Input : PSSM matrix

    Output : Feature vector(400)

    """
    bigram_matrix = []
    
    pssm_matrix = pd.read_csv(pssm_matrix)
    amino_acides = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    for x in  amino_acides:
        
        for y in amino_acides:
            
            bigram_sum =0
        
            for i in range(1,pssm_matrix.shape[0]-1):
                
                 bigram_sum = bigram_sum + (pssm_matrix[x][i]*pssm_matrix[y][i+1])
                 
            bigram_matrix.append(bigram_sum)
            
    return bigram_matrix

