# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 10:24:57 2018

@author:  Saleh, Mir, A.
"""


from Bio import SearchIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW

pro_seq = """PIVDTGSVAPLSAAEKTKIRSAWAPVYSTYETSGVDILVKFFTSTPAAQEFFPKFKGLTTADELKKSADVRWHAERIINAVDDAVASMDDTEKMSMKLRNLSGKHAKSFQVDPEYFKVLAAVIADTVAAGDAGFEKLMSMICILLRSAY"""

blast_qresult = open("../data/2LHB-BLAST.xml")
blast_record = NCBIXML.read(blast_qresult)
blast_records = NCBIXML.parse(blast_qresult)

E_VALUE_THRESH = 0.04

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


def download_bxml(fasta_string,protein_name):
    
    """
    Download XML of protein  alignments 
    
    Input : Protein sequence in FASTA format , Protein name
    
    output : Protein alignments in XML format 
    
    """
    result_handle = NCBIWWW.qblast("blastp", "pdb", fasta_string)    

    with open(protein_name + ".xml", "w") as out_handle:
        out_handle.write(result_handle.read())
