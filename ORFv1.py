#!/usr/bin/env python3

from Bio import Entrez,SeqIO,SeqUtils, SeqRecord
from Bio.Seq import Seq
import pandas as pd
import os
import numpy as np

current_dir = os.getcwd()
print(current_dir)
try:
    with open('/home/skeja/Luka_work/code_area/bioinformatics_repo/bioinformatics-projects- Python/BIopython/ORFs & CDS/experimental_sequence.txt','r') as f:
        result = SeqIO.read(f,'fasta')
        sequence = result.seq
        print(len(sequence))                
except FileNotFoundError:
    print('Error! Fasta file not found!\nPlease check the file path.')
    exit()
except ValueError:
    print('Error! This file does not appear to be a valid fasta!')
def clean_sequence(): # This needs to be fixed
        valid_nucleotides = set('TAGC')
        invalid_characters = (set(sequence.upper())-valid_nucleotides)
        if invalid_characters:
            print(f"Sequence has invalid characters {invalid_characters}")
            instead_of_N = np.random.choice(valid_nucleotides)
            for nt in sequence:
                if nt in valid_nucleotides:
                    print('no invalid characters')
                else:
                    sequence = "".join(nt for np.random.choice(valid_nucleotides) in sequence )
                          
#Some Error handling
def reverse_complement(sequence):
    sequence= Seq(sequence).upper()
    sequence_reverse = sequence.reverse_complement()
    print(type(sequence_reverse))
    return sequence_reverse
reverse_complement(sequence)

def orf_finder(sequence,strand = ['+','-'],info=False):
    no_stop_orf = []
    start_codon = ['GTG','TTG','ATG']
    stop_codon = ['TGA','TAG','TAA']
    if strand == '-':
        print(f'Strand: "{strand}"\nreverse_complement() method succesfull\n{"-"*60}')
        sequence = reverse_complement(sequence)
    all_orfs= []
    for frame in [0,1,2]:
        active_orf = []
        for j in range(frame,len(str(sequence)),3):
            codon = sequence[j:j+3]
            if codon in start_codon:
                active_orf.append({
                    'start_pos':j,
                    'start_codon':codon,
                    'sequence':[]
                })
            for orf in active_orf: #za svaki ORF otvoreni, apenduju se kodoni
                orf['sequence'].append(codon) 
            if codon in stop_codon and active_orf:
                for orf in active_orf:
                    coding_sequence = orf['sequence'][:-1] # remove the last triplet -> it's a stop codon.
                    len_coding_sequence = len(coding_sequence)
                    orf_entry = {
                        "strand":strand,
                        "frame":frame,
                        "start_pos":orf['start_pos'],
                        "stop_pos":j+3,
                        "start_codon":orf['start_codon'],
                        "stop_codon":codon,
                        "len_codons":len_coding_sequence,
                        "len_codons_nt":len_coding_sequence*3,
                        "sequence_codons":coding_sequence,
                        "complete":True
                    }
                    all_orfs.append(orf_entry)
                active_orf = []
        for orf in active_orf: #this is triggered after the loop finishes iterating over the sequence (no stop codon)
                coding_sequence_nostop = orf['sequence'][:] #complete sequnce from the START codon to the end of the sequence
                input = {
                    "strand": strand,
                    "frame": frame,
                    "start_pos": orf['start_pos'],
                    "stop_pos": None,
                    "start_codon": orf['start_codon'],
                    "stop_codon": None,
                    "len_codons": len(coding_sequence_nostop),
                    "len_codons_nt": len(coding_sequence_nostop) * 3,
                    "sequence_codons": coding_sequence_nostop,
                    "complete":False
                }
                no_stop_orf.append(input)
            
    verysmallORFs = [orf for orf in all_orfs if orf['len_codons']<=15]
    smORFs= [orf for orf in all_orfs if orf['len_codons']>15 and orf['len_codons']<=100]
    longORFs = [orf for orf in all_orfs if orf['len_codons']>100]   
    no_stop_orf = [orf for orf in no_stop_orf]
    print(f'Very small ORFs: {verysmallORFs}\nTotal number of frames: [{len(verysmallORFs)}]\n{'-'*28}')
    print(f'Small ORFs:\n {smORFs}\nTotal number of frames: [{len(smORFs)}]\n{'-'*28}')
    print(f'Long ORFs:\n {longORFs}\nTotal number of frames: [{len(longORFs)}]\n{'-'*28}')
    print(f'no STOP ORFs:\n {no_stop_orf}\nTotal number of frames: [{len(no_stop_orf)}]\n{'-'*28}')

    if info:
        print(f'total number of ORFs: {len(all_orfs)+len(no_stop_orf)} of which:\nVery small ORFs: {len(verysmallORFs)}\nSmall ORFs: {len(smORFs)}\nLong ORFs {len(longORFs)}\nno STOP ORFs {len(no_stop_orf)}')
    return all_orfs, no_stop_orf

    
if __name__ == "__main__":
    all_orfs,no_stop_orf = orf_finder(sequence=sequence,strand="+",info=True)


