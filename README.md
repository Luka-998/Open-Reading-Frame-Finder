### ORF_Finder — Prokaryotic DNA 

#### This repository serves as a clean and extensible bioinformatics workflow, suitable for further data science and machine learning applications in genomic and protein function analysis.
#### The idea is the recreation of already existing open reading frame finders (such as [https://www.ncbi.nlm.nih.gov/orffinder/]) from scratch as my personal project.
#### Knowledge learned along the way is invaluable for my career development in the field of bioinformatics and data analysis

- Minimal usage of AI for code writing. Complete code is written by myself
- AI is used only for synthetic sequence and debugging.
- Logic of the ORF_finder is related to my knowledge and deep understaning of molecular process in molecular biology.


### Project overview

__ORF_Finder__ is a compact, well-documented Python tool to find all possible open reading frames (ORFs) in a multiple prokaryotic FASTA sequence.

### `Project Goal`: Detect every start and stop combination (including ORFs that start but have no in-frame stop), capture small ORFs, translate predicted ORFs to protein sequences, and later compare translated products to known proteins (e.g., BLAST+).

__Current Phase__:
1. Improving accuracy of codon identification and storage
2. Testing workflow with long synthetic prokaryotic DNA sequences containing known open reading frames and stop codons

`Next Phase`:
3. Batch-translation of identified open reading frames
4. Comparing translated proteins to curated sequences from the Swiss-Prot database
5. Identifying potential novel genes and proteins based on sequence alignment and similarity metrics


Data Analysis & ML part:
1. Compare ORF_finder proteins to currated proteins
2. Feature engineering and data cleaning
3. Building machine learning models to predict protein functions based on extracted features


