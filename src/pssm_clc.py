# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 10:24:57 2018

@author:  Saleh, Mir, A.
"""

from os.path import join, isfile
from os import listdir
from threading import Thread
from queue import Queue
from Bio import SearchIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
import threading
import pandas as pd 
import numpy as np


def fasta_string(protein_name, protein_seq):
    
    """
    Generate FASTA string format 
    
    """
    
    fasta = ">" + protein_name + "\n" + protein_seq
    
    return fasta
    

def download_bxml(fasta_string, protein_name, save_path):
    
    """
    Download XML of protein  alignments 
    
    Input : 
        Protein sequence in FASTA format
        Protein name
        save_path: A path for saveing an XML file
    
    output : Protein alignments in XML format 
    
    """
    
    result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string)    

    with open(join(save_path, protein_name + ".xml"), "w") as out_handle:
        
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


def pfm(alignment_sbjct, pro_seq):
    
    """
     Create Postion frequency matrix 
     
     Input : Alignments match , Protein sequence 
     
     Output : PFM matrix 
    
    """
    
    protein_column = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
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
   # pfm_matrix = pfm_matrix *(1/len(alignment_sbjct))
    
    # Score(x)=(frequency+pseudocount)/N+20*(pseudocount)    
    
    # PSSM = log(Score)
    
#    pssm_matrix = np.log10((pfm_matrix + 1)*(1/(len(alignment_sbjct)+20)))

    # PSSM matrix can be computed using log2 or log10
    
    pssm_matrix = np.log2((pfm_matrix * 20 + 1) / (len(alignment_sbjct) + 1))
    
    return pssm_matrix


def dl_bxml_dataset(dataset, save_path):
    
    """
    Downloads all the XML files for a dataset
    
    Input:
        dataset: A pandas dataframe that contains protein name and its sequence
        save_path: A path for saving XML files
    """
    
    # Number of proteins in dataset
    num_protien = dataset.shape[0]
    
    # Find XML files if any
    xml_files = [f for f in listdir(save_path) if isfile(join(save_path, f))]
    
    for i in range(0, num_protien):
        
        protein_name = dataset['Protein name'][i]
        
        # Check if the XML file of protein is already downloaded
        if protein_name + '.xml' in xml_files:
            
            print("%d/%d - XML file of Protein %s ALREADY downloaded... " % (i + 1,\
                                num_protien, protein_name))
            
        else:
        
            # Step 1: Converts to FASTA fomrat in order to download from BLAST
            fasta_str = fasta_string(protein_name, dataset['Protein sequence'][i])
            
            while(True):
            
                try:
                    
                    # Step 2: Download XML file for protien
                    download_bxml(fasta_str, protein_name, save_path)
                    
                    break
                    
                # It may fail to donwload. so we try again!   
                except:
                
                    print("Failed to download protein %s! Try agian." % (protein_name))
            
            print("%d/%d - XML file of Protein %s downloaded... " % (i + 1,\
                                num_protien, protein_name))
            
            
def m_thread_dl_bxml(dataset, input_file, save_path, no_workers=1):
    
    """
    It download BLAST XML files using multi-threading with Queue
    Input:
        dataset: pandas dataframe that contains proteins
        input_file: A text file that contains proteins names
        save_path: A path for saving XML files
        no_workers: Number of threads for downloading
    """
    
    def worker(queue):
        
        """
        One thread for downloading BLAST XML file
        queue: Queue that contains protein names
        """
        
        while True:
            
            thread_id = threading.get_ident()
            
            protein_name = queue.get()
            protein_seq = dataset[dataset['Protein name'] == protein_name]['Protein sequence'].iloc[0]
            
            # Step 1: Converts to FASTA fomrat in order to download from BLAST
            fasta_str = fasta_string(protein_name, protein_seq)
            
            while(True):
            
                try:
                    
                    # Step 2: Download XML file for protien
                    download_bxml(fasta_str, protein_name, save_path)
                    
                    break
                    
                # It may fail to donwload. so we try again!   
                except:
                
                    print("Thread: %d Failed to download protein %s! Try agian." \
                          % (thread_id, protein_name))
            
            print("Thread: %d - XML file of Protein %s downloaded... " % (thread_id,\
                                 protein_name))
            
            queue.task_done()
         
    xml_list = open(input_file, 'r')
    
    # First find the number of XML files that need to be downloaded
    xml_files = [p_name.rstrip('\n') for p_name in xml_list.readlines()\
                 if not isfile(join(save_path, p_name.rstrip('\n') + '.xml'))]
    
    print("%d proteins need to be downloaded from the file." % (len(xml_files)))
    
    # A FIFO queue
    q = Queue()
    
    # Put all the XML file into Queue
    for f in xml_files:
        
        q.put(f)
        
    for i in range(no_workers):
        
        thread_i = Thread(target=worker, args=(q,))
        thread_i.setDaemon(True)
        thread_i.start()
        
    # Wait until the queue is empty
    q.join()
    xml_list.close()
     
        
def dataset_PSSM(dataset, XML_path, out_path):
    
    """
    It calculates the PSSM matrix of every protein in dataset and saves it in a
    CSV file in the specified path
    
    Input:
        dataset: A pandas dataframe that contains protein name and its sequence
        XML_path: Path to the XML files of the dataset
        out_path: A path for saving CSV files of PSSM matrix
        
    """
    
    # Number of proteins in dataset
    num_protien = dataset.shape[0]
    
    # Find XML files if any
    CSV_files = [f for f in listdir(out_path) if isfile(join(out_path, f))]

    for i in range(0, num_protien):
        
        protein_name = dataset['Protein name'][i]
        
        # Check if the CSV file of protein is already created
        if protein_name + '.csv' in CSV_files:
            
            print("%d/%d - CSV file of protein %s ALREADY created... " % (i + 1,\
                                num_protien, protein_name))
        
        else:
            
            # Step 1: Find aligments of the protein
            match, subject = find_aligments(join(XML_path, protein_name + ".xml"), \
                                            dataset['Protein sequence'][i])
            
            # Step 2: Create Postion frequency matrix of the protein
            pfm_matrix = pfm(subject, dataset['Protein sequence'][i])
            
            # Step 3: Generate the PSSM matrix of the protein
            pssm_matrix = pssm(pfm_matrix, subject)
            
            # Step 4: Save the CSV file of the PSSM matrix
            pssm_matrix.to_csv(join(out_path, protein_name + ".csv" ), sep=',', \
                               float_format='%f', index=False, header=False)
            
            print("%d/%d - CSV file of PSSM matrix for the Protein %s created... " % \
                  (i + 1, num_protien, protein_name))
        

if __name__ == '__main__':
    
    protein_dtfrm = pd.read_csv(r"./dataset/RDD/RDD_raw.csv")
    
    
    
#    dataset_PSSM(protein_dtfrm, './BXML/tg_bxml/', './PSSM/tg_PSSM/')
    
    #dl_bxml_dataset(protein_dtfrm, './BXML/tg_bxml/')
    
    m_thread_dl_bxml(protein_dtfrm, 'RDD-bxml.txt', './BXML/rdd_bxml/', 5)
    
    #fasta_str = fasta_string(protein_dtfrm['Protein name'][1],protein_dtfrm['Protein sequence'][1])
    
    #download_bxml(fasta_str,protein_dtfrm['Protein name'][1])
    
#    match, subject = find_aligments(join('./dd_bxml/', \
#                                    protein_dtfrm['Protein name'][1] + ".xml"),
#                                    protein_dtfrm['Protein sequence'][1])
#    
#    pfm_matrix = pfm(subject,protein_dtfrm['Protein sequence'][1])
#    
#    pssm_matrix = pssm(pfm_matrix,subject)
#    pssm_matrix.to_csv(protein_dtfrm['Protein name'][1]+".csv", sep=',',
#                       float_format='%f', index=False, header=False)

    
    


