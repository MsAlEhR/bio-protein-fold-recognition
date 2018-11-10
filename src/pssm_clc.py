# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 10:24:57 2018

@author:  Saleh, Mir, A.
"""


from Bio import SearchIO
from Bio.Blast import NCBIXML

pro_seq = """PIVDTGSVAPLSAAEKTKIRSAWAPVYSTYETSGVDILVKFFTSTPAAQEFFPKFKGLTT
ADELKKSADVRWHAERIINAVDDAVASMDDTEKMSMKLRNLSGKHAKSFQVDPEYFKVLA
AVIADTVAAGDAGFEKLMSMICILLRSAY"""

blast_qresult = open("../data/2LHB-BLAST.xml")
blast_record = NCBIXML.read(blast_qresult)
blast_records = NCBIXML.parse(blast_qresult)

E_VALUE_THRESH = 0.04

alignments = []

for alignment in blast_record.alignments:

      for hsp in alignment.hsps:

           if hsp.expect < E_VALUE_THRESH:
                
                alignments.append(str(hsp.query[0:]))
                

                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print("e value:", hsp.expect)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")