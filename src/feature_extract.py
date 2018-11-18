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
    amino_acides = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    pssm_matrix = pd.read_csv(pssm_matrix,sep=',',names =amino_acides)
    
    
    for x in  amino_acides:
        
        for y in amino_acides:
            
            bigram_sum =0
        
            for i in range(1,pssm_matrix.shape[0]-1):
                
                 bigram_sum = bigram_sum + (pssm_matrix[x][i]*pssm_matrix[y][i+1])
                 
            bigram_matrix.append(bigram_sum)
            
    return bigram_matrix


def bigram_dataset(protein_dtframe):
    
    columns = ['F' + str(i) for i in range(1,401)]
    
    data_set = pd.DataFrame(columns=['Fold','Protein name'] + columns)
    
    for num, i in enumerate(protein_dtframe['Protein name']):
        
        pssm_path = r"./dd_PSSM" + "/" + i  + ".csv"
        
        feature_vector = bigram(pssm_path)
        
        data_set.loc[i] =  [protein_dtframe.loc[num]['Fold'], i] + feature_vector
        
        print("%d/%d - Feature vector of protein %s created... " % (num + 1, \
                                        protein_dtframe.shape[0], i))
        
    data_set.to_csv('./dataset/dd_pssm_dataset.csv', index=False)
    
    return data_set

if __name__ == '__main__':
    
    protein_dtfrm = pd.read_csv(r"./dataset/DD_raw.csv")
    
    data = bigram_dataset(protein_dtfrm)    