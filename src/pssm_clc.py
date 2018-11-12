# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 10:24:57 2018

@author:  Saleh, Mir, A.
"""


from Bio import SearchIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
import pandas as pd 
import numpy as np


def fasta_string(protein_name,protein_seq):
    
    """
    Generate FASTA string format 
    
    """
    fasta = ">"+protein_name+"\n"+protein_seq
    
    return fasta
    


def download_bxml(fasta_string,protein_name):
    
    """
    Download XML of protein  alignments 
    
    Input : Protein sequence in FASTA format , Protein name
    
    output : Protein alignments in XML format 
    
    """
    result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string)    

    with open(protein_name + ".xml", "w") as out_handle:
        out_handle.write(result_handle.read())


def find_aligments(blast_xml, pro_seq, e_thresh=0.04):
    
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
    alignment_sbjct = []
    
    # Find alignments of protein sequence and extract 
    for alignment in blast_record.alignments:
    
          for hsp in alignment.hsps:
    
               if hsp.expect < E_VALUE_THRESH:
                   
               # only consider protein sequence query   
                   if hsp.query == pro_seq:
    
                    alignments.append(str(hsp.sbjct[0:]))
                    alignment_match.append(hsp.match)
                    alignment_sbjct.append(hsp.sbjct)
                    
    blast_qresult.close()
                    
    return alignment_match, alignment_sbjct

def pfm(alignment_sbjct,pro_seq):
    
    """
     Create Postion frequency matrix 
     
     Input : Alignments match , Protein sequence 
     
     Output : PFM matrix 
    
    """
    
    protein_column = list(set(pro_seq))
    
    pfm_matrix = pd.DataFrame(np.zeros((len(pro_seq),len(protein_column))),columns=protein_column )
    
    seq_len = len(pro_seq)

    for amino in range(0, seq_len):    

        for alignm in alignment_sbjct:
            
            if alignm[amino] in protein_column:
            
                pfm_matrix[alignm[amino]][amino] = pfm_matrix[alignm[amino]][amino] + 1
        
    
    return pfm_matrix      

def pssm(pfm_matrix,alignment_sbjct):
    
    """
    Generate PSSM matrix
    
    Input : PFM , Alignment subject 
    
    Output : PSSM
    
    
    """
    # Amino acid frequency at every position     
    pfm_matrix = pfm_matrix *(1/len(alignment_sbjct))
    
    # Score(x)=(frequency+pseudocount)/N+20*(pseudocount)    
    
    # PSSM = log(Score)
    
    pssm_matrix = np.log10((pfm_matrix + 1)*(1/(len(alignment_sbjct)+20)))
    
    return pssm_matrix


if __name__ == '__main__':
    
    protein_dtfrm = pd.read_csv(r"./dataset/DD_raw.csv")
    
    fasta_str = fasta_string(protein_dtfrm['Protein name'][1],protein_dtfrm['Protein sequence'][1])
    
    download_bxml(fasta_str,protein_dtfrm['Protein name'][1])
    
    match, subject = find_aligments(protein_dtfrm['Protein name'][1]+".xml",protein_dtfrm['Protein sequence'][1])
    
    pfm_matrix = pfm(subject,protein_dtfrm['Protein sequence'][1])
    
    pssm_matrix = pssm(pfm_matrix,subject)
    
    


