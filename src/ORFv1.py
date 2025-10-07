#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import os
import numpy as np  
from pathlib import Path
import argparse


def parser():
    parser = argparse.ArgumentParser(
        prog='CLI-find-FASTA',
        description='Easy to use',
        usage='CLI friendly argument parser'  
    )
    parser.add_argument('--fasta',required=True,help='Required argument for Path to input FASTA') #CLI user input
    parser.add_argument('--out',required=False,help='Optional arg. to save the output')
    args=parser.parse_args()
    return args
args = parser()
def get_output_path(args): # User defined. njih samo prosledjujem dole da ih cekiram u sistemu.
    standard_path_fasta = Path("..")/ "results" / "ORF.fasta"
    standard_path_summary =  Path("..")/ "results" / "ORF_summary.txt"
    if args.out:
        base = Path(args.out)
        fasta_out = base.with_suffix('.fasta')
        summary_out = base.with_suffix('.txt')
        if base.suffix.lower() not in ['.fasta','.fa','.fas']:
            print('Warning: file may not be a FASTA type.')     
            if base.exists() and base.is_dir(): #if user provides a dir path, default filename will be used
                print(f"Defined path {base} is a dir.\nRegular filename {standard_path_fasta} and {standard_path_summary} will be used.\n")
                return standard_path_summary,standard_path_fasta
            if base.exists() and base.is_file(): #OVDE SAM STIGAO!
                print(f"\nPath exists.Do you want to overwrite it? Y/N?") 
                choice = input().strip().upper()
                if choice == "Y":
                    print(f"\File will be saved at \t{fasta_out}|\t{summary_out}.")
                    return fasta_out,summary_out
                if choice == "N":
                    print(f"\nFile will be saved at \t{standard_path_fasta}|\t{standard_path_summary}")
                    return standard_path_summary,standard_path_fasta
            if not base.parent.exists(): #--out ..\wtf\text1.txt check if 'wtf' exists (parent)
                print(f'\nProvided file path: {base.parent} does not exist.\nDo you want me to create directory and file?\tY or N')
                input2 = str(input()).strip().upper()
                while input2 not in ("Y","N"):
                    print("Please choose Y or N only.")
                    input2 = str(input()).strip().upper()
                if input2 == "Y": 
                    base.parent.mkdir(parents=True,exist_ok=True)
                    print(f" {base} successfully created!")
                    return base
                if input2 =="N":
                    print(f"File will be saved at regular path:\n{standard_path_fasta.parent}")
                    return standard_path_fasta,standard_path_summary
            else:
                print(f"Parent directory {base.parent} exists, file will be created as: {standard_path_fasta.parent}")
                return standard_path_fasta.parent
    else:   
        print(f'\nNo output path provided.\nUsing default: {standard_path_fasta.parent}')
        return standard_path_fasta.parent
output_path=get_output_path(args)

def load_fasta(args):
    path = Path(args.fasta)
    if not path.exists():
        print(f"\nFile {path} not Found! Please check path to FASTA again")
    else:
        print("File exists, method continues!")
        sequence = SeqIO.read(path,'fasta')
        return sequence
sequence = load_fasta(args)

def clean_sequence(seq): # IUPAC nucleotide code [https://www.bioinformatics.org/sms/iupac.html]
        valid_nucleotides = set('TAGC')
        invalid_nucleotides = (set(str(seq.upper()))-valid_nucleotides)
        cleaned = []
        for nt in seq.upper():
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

        return "".join(cleaned)
sequence = clean_sequence(sequence)
print(sequence)
                  
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
    standard_path_fasta = Path("..")/ "results" / "ORF.fasta"
    standard_path_summary =  Path("..")/ "results" / "ORF_summary.txt"
    all_orfs,no_stop_orf = orf_finder(sequence=sequence,strand="+",info=False)
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

    SeqIO.write(complete_sequence,fasta_output_path,'fasta')
    with open(summary_output_path,'w') as f1:
        f1.write(f'Total ORFs: {len(all_orfs)}')


