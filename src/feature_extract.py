# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 11:41:08 2018

@author: Saleh, Mir, A.
"""
import pandas as pd

def bigram(pssm_matrix,k=1):

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


def create_dataset(protein_dtframe,function,k,dataset_name):
    
    columns = ['F' + str(i) for i in range(1,1601)]
    
    data_set = pd.DataFrame(columns=['Fold','Protein name'] + columns)
    
    for num, i in enumerate(protein_dtframe['Protein name']):
        
        pssm_path = r"./PSSM/rdd_PSSM" + "/" + i  + ".csv"
        
        feature_vector = function(pssm_path,k)
        data_set.loc[i] =  [protein_dtframe.loc[num]['Fold'], i] + feature_vector
        
        print("%d/%d - Feature vector of protein %s created... " % (num + 1, \
                                        protein_dtframe.shape[0], i))
        
    data_set.to_csv('./dataset/'+dataset_name+'.csv', index=False)
    
    return data_set

def separated_dimer(pssm_matrix,K):
    
    """
    Get features by Separated dimers method

    Input : PSSM matrix , K (ditance between dimers)

    Output : Feature vector(400)

    """
    spr_matrix = []
    amino_acides = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    pssm_matrix = pd.read_csv(pssm_matrix,sep=',',names =amino_acides)
    
    
    for x in  amino_acides:
        
        for y in amino_acides:
            
            spr_sum =0
        
            for i in range(1,pssm_matrix.shape[0]-K):
                
                 spr_sum = spr_sum + (pssm_matrix[x][i]*pssm_matrix[y][i+K])
                 
            spr_matrix.append(spr_sum)
            
    return spr_matrix    
    
def acc_transformation(pssm_matrix,K):
    
    """
    Get Feature by ACC Transformation 
    
    Input : PSSM matrix , K = LG
    
    Output : Feature vector(20*LG + 380*LG)
    """
    
    amino_acides = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                    'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        
    pssm_matrix = pd.read_csv(pssm_matrix,sep=',',names =amino_acides)
    
    AC_matrix = []        
    CC_matrix = []
    
    for k in range(1,K+1) :
        
        
        # calculate eh mean of any amino acid clomn in every 
    
        average = pssm_matrix[amino_acides].mean()
        
        for x in amino_acides:
                
                ac_sum = 0
                
                for i in range(1,pssm_matrix.shape[0]-k):
                    
                    ac_sum = ac_sum + (pssm_matrix[x][i] - average[x]) * (pssm_matrix[x][i+k]-average[x])/(pssm_matrix.shape[0]-k)
                        
                AC_matrix.append(ac_sum)
        
        
        for x in amino_acides:
          
            for y in amino_acides :
                
                if x == y :
                    continue
                
                cc_sum =0
                
                for i in range(1,pssm_matrix.shape[0]-k):
                    
                    cc_sum = cc_sum + (pssm_matrix[x][i] - average[x]) * (pssm_matrix[y][i+k] - average[y])/(pssm_matrix.shape[0]-k)
                CC_matrix.append(cc_sum)
    
    return AC_matrix + CC_matrix
    
        

if __name__ == '__main__':
    
    protein_dtfrm = pd.read_csv(r"./dataset/RDD/RDD_raw.csv")
    
    data = create_dataset(protein_dtfrm,acc_transformation,4,'rdd_acc_dataset')  
    
    
    