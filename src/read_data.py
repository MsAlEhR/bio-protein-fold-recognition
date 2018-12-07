# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 14:30:43 2018

@author:  Saleh, Mir, A.
"""


import pandas as pd
import numpy as np
import bio_pre as bp
import time
import re # Regular expression for text processing in dataset

start_t = time.time()

# **************** Convert text data file to dictionary ************

with open("../data/TG.dataset") as f:
    content = f.readlines()

content = [x.strip() for x in content]

type_protein = []
protein_num = []
protein_name = []

protein = {}

# protein name status
state_pn = False

for num, line in enumerate(content):

    if bool(re.match('TYPE|FOLD\s*\(\d+\)', line)):

        fold_type = line.replace(" ", "")
        
        protein[fold_type] = {}

        state_pn = False

    if ">" in line:

        state_pn = True
        name_pr = line.replace(" ", "")

        protein[fold_type][name_pr] = ''

    elif state_pn is True:

        protein[fold_type][name_pr] = protein[fold_type][name_pr] + line
        
    #time.sleep(1)
        
# ********************************************************************

# ************ To convert dictionary data to Pandas DataFrame ********
# Added by Mir, A. 
       
# Convert dictionary data to list
# Name of each fold is changed so that Type (?) is removed from the string
        
list_data = [[re.sub('TYPE|FOLD\s*\(\d+\)', '', fold).strip(), pr[1:], \
              protein[fold][pr]] for fold in list(protein.keys()) \
              for pr in list(protein[fold].keys())]

data_frame = pd.DataFrame(list_data, columns=['Fold', 'Protein name', \
                                              'Protein sequence'])
    
# Save this dataframe to CSV file
data_frame.to_csv('./dataset/TG_raw.csv', index=False, header=True)

# ******************************************************************
    
# ********** Feature extraction section *************************
    
#data_frame['FV'] = data_frame['Protein sequence'].map(lambda x: bp.generate_FV(x))
#
#data_frame['count']= data_frame['FV'].map(lambda x: len(x))
#
#for i in range(21):
#
#    column_name = "Feature" + str(i)
#    data_frame[column_name]= data_frame['Protein sequence'].map(lambda x: bp.generate_FV(x)[i])
    

# Saving Pandas dataframe to CSV
#data_frame[['Fold', 'Protein name' ] + ['Feature%d' % i for i in range(21)]].to_csv('./dataset/DD_dataset.csv',
#           index=False , header=True)

# ****************************************************************


print("Finished: %.2f sec" % (time.time() - start_t))
