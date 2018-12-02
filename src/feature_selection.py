# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 11:11:38 2018

@author: Saleh, Mir, A.
"""

import csv
import re
import pandas as pd


def gen_slc_feature(data_name,slc_feature,output_name,no_feature):
    
    """
     Generate selected feature from dataset
     
     Input : Dataset name , Selected feature , Output name, Number of most possible feature
     
     Output :  Dataset with selected feature in .csv
    
    """
    dataset = pd.read_csv('./dataset/'+data_name,sep=",")
    
    slct_feature = pd.read_csv('./dataset/DD/feature_selection/'+slc_feature,sep='\t',header=None)
    
    dataset[['Fold','Protein name']+list(slct_feature[0][:no_feature])].to_csv('./dataset/'+output_name,sep=',',index=False)
    
    return  dataset[['Fold','Protein name']+list(slct_feature[0][:no_feature])]


def fusion_feature(data_name1, data_name2, merged_data_name):
    
    
    dataset1 = pd.read_csv('./dataset/' + data_name1, sep=",", skiprows=[0],
                           header=None)
    
    dataset2 = pd.read_csv('./dataset/' + data_name2, sep=",")
    
    # Find number of features in dataset2
    d2_features = [col for col in list(dataset2.columns.values) if bool(re.match('F\d+',
                   str(col)))]
    
    # Number of features in merged dataframe    
    new_dim = (dataset1.shape[1] - 2) + len(d2_features)
    
        
    result = pd.concat([dataset1, dataset2[d2_features]], axis=1, sort=False)
    
    result.to_csv('./dataset/'+ merged_data_name, index=False, \
                  header=['Fold', 'Protein name'] + ['F%d' % i for i in range(1, new_dim + 1)])
    
    return result 

if __name__ == "__main__":
    
#    df = fusion_feature("ACC_dataset.csv","dd_separated_dimer_dataset.csv", 'dd_fusion_acc_dimer_FULL.csv')
    
    dtst = gen_slc_feature('dd_fusion_acc_dimer_FULL.csv','IG_feature_acc_dimer.txt','ACC_dimer_dataset_IG.csv',500)