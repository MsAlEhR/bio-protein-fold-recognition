# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 14:30:43 2018

@author:  Saleh, Mir, A.
"""


import pandas as pd
import bio_pre as bp

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


# Added by Saleh, R. 

# Convert  Protein Sequence   
data_frame['Protein list']=data_frame['Protein sequence'].map(lambda x : bp.Sequence_List(x) )   

# Calculate Compostion 
data_frame['Composition_polar'] = data_frame['Protein list'].map(lambda x: bp.Composition_Cal(x)[1][0])
data_frame['Composition_neutral'] = data_frame['Protein list'].map(lambda x: bp.Composition_Cal(x)[1][1])
data_frame['Composition_hydrophoblic'] = data_frame['Protein list'].map(lambda x: bp.Composition_Cal(x)[1][2])

data_frame['Composition_polar_count'] = data_frame['Protein list'].map(lambda x: bp.Composition_Cal(x)[0][0])
data_frame['Composition_neutral_count'] = data_frame['Protein list'].map(lambda x: bp.Composition_Cal(x)[0][1])
data_frame['Composition_hydrophoblic_count'] = data_frame['Protein list'].map(lambda x: bp.Composition_Cal(x)[0][2])


# Calculation Transition 
data_frame['Transition_polar'] = data_frame['Protein list'].map(lambda x: bp.Transition_Cal(x)[0])
data_frame['Transition_neutral'] = data_frame['Protein list'].map(lambda x: bp.Transition_Cal(x)[1])
data_frame['Transition_hydrophoblic'] = data_frame['Protein list'].map(lambda x: bp.Transition_Cal(x)[2])

# Calculation Hydrophoblic 
data_frame['Distribution_polar'] = data_frame[['Protein list','Composition_polar_count']].map(lambda x,y: bp.Distributon_Cal(x,y)[0])


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
