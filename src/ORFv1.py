#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import os
import numpy as np  
from pathlib import Path
import argparse
from utils.get_output_path import parser,load_fasta,get_output_path

args = parser() #ok
sequence_storage = load_fasta(args)
csv_path = Path("..")/ "results" / "fasta.csv"
def clean_sequence(sequence_storage): # IUPAC nucleotide code [https://www.bioinformatics.org/sms/iupac.html]
        for i in sequence_storage:
            sequence = i['Seq']
            valid_nucleotides = set('TAGC')
            cleaned = []
            for nt in str(sequence).upper():
                if nt in valid_nucleotides:
                    cleaned.append(nt)
                elif nt == 'N':
            # unknown, pick any nucleotide
                    cleaned.append(np.random.choice(list(valid_nucleotides)))
                elif nt == 'R':
                    cleaned.append(np.random.choice(['A', 'G']))
                elif nt == 'Y':
                    cleaned.append(np.random.choice(['C', 'T']))
                elif nt == 'S':
                    cleaned.append(np.random.choice(['G', 'C']))
                elif nt == 'W':
                    cleaned.append(np.random.choice(['A', 'T']))
                elif nt == 'K':
                    cleaned.append(np.random.choice(['G', 'T']))
                elif nt == 'M':
                    cleaned.append(np.random.choice(['A', 'C']))
                else:
            # any truly invalid character
                    cleaned.append(np.random.choice(list(valid_nucleotides)))

            yield "".join(cleaned)
cleaned = clean_sequence(sequence_storage)


#Some Error handling

def reverse_complement(sequence):
    for i in sequence:
        sequence = i['Seq']
        sequence= Seq(sequence).upper()
        sequence_reverse = sequence.reverse_complement()
        print(type(sequence_reverse))
        yield sequence_reverse
reverse_complement(cleaned)

def orf_finder(cleaned,strand = '+',info=True):
    for sequence in cleaned:
        print(sequence)
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
                for orf in active_orf: # appending codons for every open ORF
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
    print(f'Very small ORFs: {verysmallORFs}\nTotal number of frames: [{len(verysmallORFs)}]\n' + '-'*28)
    print(f'Small ORFs:\n {smORFs}\nTotal number of frames: [{len(smORFs)}]\n' + '-'*28)
    print(f'Long ORFs:\n {longORFs}\nTotal number of frames: [{len(longORFs)}]\n' + '-'*28)
    print(f'no STOP ORFs:\n {no_stop_orf}\nTotal number of frames: [{len(no_stop_orf)}]\n' + '-'*28)


    if info:
            print(f'total number of ORFs: {len(all_orfs)+len(no_stop_orf)} of which:\nVery small ORFs: {len(verysmallORFs)}\nSmall ORFs: {len(smORFs)}\nLong ORFs {len(longORFs)}\nno STOP ORFs {len(no_stop_orf)}')
            if strand == ['-']:
                print("Reverse complement generated successfully.")

    return all_orfs, no_stop_orf


  
if __name__ == "__main__":
    fasta_output_path,summary_output_path = get_output_path(args)
    standard_path_fasta = Path("..")/ "results" / "ORF.fasta"
    standard_path_summary =  Path("..")/ "results" / "ORF_summary.txt"
    all_orfs,no_stop_orf = orf_finder(cleaned,strand="+",info=False)
    complete_sequence = []
    sequence_metadata = []
    
    for i,orf in enumerate(all_orfs):
        joined_seq = "".join(orf['sequence_codons'])
        record_info = ({
            'seq':orf['sequence_codons'],
            'label':f'ORF_{i+1}',
            'strand':orf['strand'],
            'frame':orf['frame'],
            'start':orf['start_pos'],
            'stop':orf['stop_pos'],
            'len_codons':orf['len_codons'],
            'complete':orf['complete']
        }) 
        sequence_metadata.append(record_info)

        seq_record = SeqRecord(Seq(
            joined_seq),
            id = record_info['label'],
            description=f"strand={record_info['strand']} frame={record_info['frame']} start={record_info['start']} stop={record_info['stop']} len_codons={record_info['len_codons']} complete={record_info['complete']}"
        )
        complete_sequence.append(seq_record)

with open(fasta_output_path,'w') as f:
    SeqIO.write(seq_record,f,'fasta')
with open(summary_output_path,'w') as s:
    s.write(f'Number of detected ORFs:{len(all_orfs)+len(no_stop_orf)} of which ORFs without stop codon = {len(no_stop_orf)}')       
   


