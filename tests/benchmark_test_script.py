#!/usr/bin/env python3
import os
from src.ORFv1 import load_fasta, clean_sequence, reverse_complement, orf_finder

#paths

fasta = 'synthetic_sequence.txt'
fasta_path = os.path.join(os.getcwd(),'tests_sequences',fasta)
output_path = os.path.join(os.getcwd(),'results')
sequence = load_fasta(fasta_path)
sequence = clean_sequence(sequence)
all_orfs, no_stop_orfs = orf_finder(sequence=sequence,strand = '+',info=True)
