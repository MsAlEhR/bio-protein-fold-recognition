# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 13:37:12 2018

@authors: Saleh, Mir, A.
"""


pro_seq = "ASFSEAPPGNPKAGEKIFKTKCAQCHTVDKGAGHKQGPNLNGLFGRQSGTTPGYSYSTADKNMAVI" \
          "WEENTLYDYLLNPKKYIPGTKMVFPGLKKPQERADLISYLKEATS"
          
def Sequence_List(pro_seq):
    
        """
        This function genrates List of amino acid in the protein sequence
        
        +++ Only for Hydrophobicity property +++
    
        Input:Protein_sequence
    
        Return:Protein_list
        """

        pro_list = []        
        polar = ['R', 'K', 'E', 'D', 'Q', 'N']
        neutral = ['G', 'A', 'S', 'T', 'P', 'H', 'Y']
        hydrophoblic = ['C', 'V', 'L', 'I', 'M', 'F', 'W']
        
        for i in pro_seq:
       
            if i in polar:
  
                pro_list.append(1)

            if i in neutral:

                pro_list.append(2)

            if i in hydrophoblic:

                pro_list.append(3)

        return pro_list


def Composition_Cal(pro_list):

        """
        Composition Calculatuin

        Input : Protein_list

        Retrun : Count and Percnet of each Polar, Neutral, Hydrophoblic
        """

        polar_count = 0
        neutral_count = 0
        hydrophoblic_count = 0

        for i in pro_list:

            if i == 1:

                polar_count = polar_count + 1

            if i == 2: 

                neutral_count = neutral_count + 1

            if i == 3:
      
                hydrophoblic_count = hydrophoblic_count + 1

        percent = [(polar_count / len(pro_list)) * 100,
                   (neutral_count / len(pro_list)) * 100,
                   (hydrophoblic_count / len(pro_list)) * 100]

        count = [polar_count, neutral_count, hydrophoblic_count]

        return count, percent


def Transition_Cal(pro_list):

        """
        Transition Calculatuin

        Input : Protein_List

        Retrun : Transition Parameter
        """

        transition=[]

        for i, e in enumerate(pro_list):
            
            if i+1 != len(pro_list):
                if pro_list[i] == 1 and pro_list[i+1] == 2:
                    transition.append("pn")
                elif pro_list[i] == 1 and pro_list[i+1] == 3:
                    transition.append("ph")
                elif pro_list[i] == 2 and pro_list[i+1] == 3:
                    transition.append("nh")

        trans_list=[(transition.count("pn") / (len(pro_list)-1)) * 100,
                    (transition.count("ph") / (len(pro_list)-1)) * 100,
                    (transition.count("nh") / (len(pro_list)-1)) * 100]


        return trans_list

def Distributon_Cal(pro_list, group_count):

        """
        Distribution Calculatuin

        Input : Protein_list , Group_count(Calculate from Compostion_cal)

        Return : Distribution Parameter for each group
        """

        group_distribution = [i for i, e in enumerate(pro_list) if e == 1]
        
        group_distribution_percent = [((e+1)/len(pro_list))*100 for i, e in enumerate(group_distribution)
                                      if (i == 0 or i == int(group_count/4)-1 or
                                          i == int((group_count/4)*2)-1 or
                                          i == int((group_count/4)*3)-1 or
                                          i == group_count-1)]
        return group_distribution_percent


def generate_FV(pro_seq):
    
    """
    Generate a feature vector from a protein sequence
    Input : Protein sequence (String)
    
    Retrun: Feature vector (List)
    
    """
    

if __name__ == '__main__':
    
    pro_list=Sequence_List(pro_seq)
    x,y=Composition_Cal(pro_list)
    d=Transition_Cal(pro_list)
    k=Distributon_Cal(pro_list,x[0])
    
