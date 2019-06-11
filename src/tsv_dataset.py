# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 15:07:56 2018

@author: Saleh , Mir, A. 

"""

from os.path import splitext
import pandas as pd 

def create_tsv(data_path, dataset_name):
   
    """
      Create .tsv files for label and features
      
      Input : Dataset name , and number of dataset
      
      Output : Label , feature in .tsv files
    
    """

    
    dataset = pd.read_csv(data_path + dataset_name)
    
    no_samples, no_features = dataset.shape
    no_features = no_features - 2
    
     
    y_true, labels = pd.factorize(dataset.Fold)
    dataset.insert(1, 'class labels', y_true)
    # In order to remove white space in protein name 
    dataset['Protein name'] = dataset['Protein name'].map(lambda x : x.replace(" " ,""))
    
    dataset[['Protein name','class labels']].to_csv(data_path + splitext(dataset_name)[0] + \
           '_label' + '.txt', index=False, header=False, sep="\t")
    
    dataset[['Protein name']+['F%d' % i for i in range(1, no_features + 1)]].to_csv(data_path + \
            splitext(dataset_name)[0] + '_feature' + '.tsv', index=False, sep="\t")
    
    return dataset,dataset[['Protein name','class labels']],dataset[['Protein name']+ \
                           ['F%d' % i for i in range(1, no_features + 1)]]

if __name__ == "__main__":
         
    dataset, label , feature = create_tsv('./dataset/EDD/', 'edd_fusion_accK4_dimerK8_FULL.csv')