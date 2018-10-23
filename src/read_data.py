# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 14:30:43 2018

@author:  Saleh, Mir, A.
"""

import pandas as pd

with open("DD-train.dataset.txt") as f:
    content = f.readlines()

content = [x.strip() for x in content]

type_protein = []
protein_num = []
protein_name = []

protein = {}

# protein name status
state_pn = False

for num, line in enumerate(content[7:]):

    if "TYPE" in line:

        fold_type = line

        protein[line] = {}

        state_pn = False

    if ">" in line:

        state_pn = True
        name_pr = line

        protein[fold_type][name_pr] = ''

    elif state_pn is True:

        protein[fold_type][name_pr] = protein[fold_type][name_pr] + line

# Added by Mir, A.
# To convert dictionary data to Pandas DataFrame
        
# Convert dictionary data to list
list_data = [[fold, pr, protein[fold][pr]] for fold in list(protein.keys()) \
             for pr in list(protein[fold].keys())]

data_frame = pd.DataFrame(list_data, columns=['Fold', 'Protein name', \
                                              'Protein sequence'])
################


#for num,line in enumerate(content):
#
#    if "TYPE" in line:
#        
#        type_protein.append(line)
#        continue
#    
#    if ">" in line:
#        
#        protein_num.append(num)
#        protein_name.append(line)
#        continue
#    
#for i, e in enumerate(protein_num):
#    
#    if i + 1 != len(protein_num):
#        
#        protein[content[e]] = ''.join(content[e + 1:protein_num[i + 1]])
