# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 11:11:38 2018

@author: Saleh, Mir, A.
"""

import pandas as pd


def gen_slc_feature(data_name,slc_feature,output_name,no_feature):
    
    """
     Generate selected feature from dataset
     
     Input : Dataset name , Selected feature , Output name, Number of most possible feature
     
     Output :  Dataset with selected feature in .csv
    
    """
    dataset = pd.read_csv('./dataset/'+data_name,sep=",")
    
    slct_feature = pd.read_csv('./dataset/DD/feature_selection/'+slc_feature,sep='\t',header=None)
    
    dataset[['Fold','Protein name']+list(slct_feature[0][:no_feature])].to_csv('./dataset/'+output_name,sep=',')
    
    return  dataset[['Fold','Protein name']+list(slct_feature[0][:no_feature])]


if __name__ == "__main__":
    
    dtst = gen_slc_feature('ACC_dataset.csv','ACC_IG_feature.txt','ACC_dataset_IG.csv',500)