# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 10:24:57 2018

@author:  Saleh, Mir, A.
"""


from Bio import SearchIO
from Bio.Blast import NCBIXML

pro_seq = """PIVDTGSVAPLSAAEKTKIRSAWAPVYSTYETSGVDILVKFFTSTPAAQEFFPKFKGLTTADELKKSADVRWHAERIINAVDDAVASMDDTEKMSMKLRNLSGKHAKSFQVDPEYFKVLAAVIADTVAAGDAGFEKLMSMICILLRSAY"""


def find_aligments(blast_xml, e_thresh=0.04):
    
    """
    It finds Aligments for a givem BLAST XML file
    
    Input:
        blast_xml: BLAST XML file
        
    return:
        Alighments Query and Match (List, List)
    """
    
    blast_qresult = open(blast_xml)
    
    blast_record = NCBIXML.read(blast_qresult)
    blast_records = NCBIXML.parse(blast_qresult)
    
    E_VALUE_THRESH = e_thresh
    
    alignments = []
    alignment_match = []
    alignment_query = []
    
    # Find alignments of protein sequence and extract 
    for alignment in blast_record.alignments:
    
          for hsp in alignment.hsps:
    
               if hsp.expect < E_VALUE_THRESH:
                   
               # only consider protein sequence query   
                   if hsp.query == pro_seq:
    
                    alignments.append(str(hsp.sbjct[0:]))
                    alignment_match.append(hsp.match)
                    alignment_query.append(hsp.query)
                    
    blast_qresult.close()
                    
    return alignment_match, alignment_query


if __name__ == '__main__':
    
    query, match = find_aligments("../data/2LHB-BLAST.xml")
