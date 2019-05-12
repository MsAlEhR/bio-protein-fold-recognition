# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 22:55:13 2019

@author: SaLeH
"""
import pandas as pd 

# write function that comppute number of each dataset 
# write function to under stand what ffeature more hetrohgenous 
# write some thing that define why in data sert 7 , 4, distance between deatrue of pssm can help
# to imporve accuracy 


def compute_distribution(path,rg):
    
    dataset = pd.read_csv('./dataset/RDD/'+ path,'	',header=None)
    
    dataset=dataset.loc[dataset[1]>rg]
    
    dataset[0]=dataset[0].apply(lambda x : int(x[1:]))
    
    numberACC=dataset[dataset[0]>1600][0].count()
    
    numberSEP=dataset[dataset[0]<1600][0].count()
    
    return dataset,numberACC,numberSEP

dataset,acc,sep=compute_distribution('IG_RDD_FEATURE.txt',0.5)

    