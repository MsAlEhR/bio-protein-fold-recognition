# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 13:37:12 2018

@author: User
"""


pro_seq="ASFSEAPPGNPKAGEKIFKTKCAQCHTVDKGAGHKQGPNLNGLFGRQSGTTPGYSYSTADKNMAVIWEENTLYDYLLNPKKYIPGTKMVFPGLKKPQERADLISYLKEATS"

pro_list=[]

polar = ['R','K','E','D','Q','N']

neutral = ['G','A','S','T','P','H','Y']

hydrophoblic = ['C','V','L','I','M','F','W']


for i in pro_seq :
    
    if i in polar :
       pro_list.append(1)
    if i in neutral:
        pro_list.append(2)
    if i in hydrophoblic :
        pro_list.append(3)
        
polar_count=0
neutral_count=0
hydrophoblic_count=0

for i in pro_list :
     if i==1 :
        polar_count=polar_count+1
     if i==2 :   
         neutral_count=neutral_count+1
     if i==3:
         hydrophoblic_count=hydrophoblic_count+1


print("poral :",(polar_count/len(pro_list))*100,"\n neutral :",(neutral_count/len(pro_list))*100,
      "\n hydrophoblic :",(hydrophoblic_count/len(pro_list))*100)

         