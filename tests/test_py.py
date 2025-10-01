from src.ORFv1 import load_fasta, clean_sequence, reverse_complement, orf_finder
path = 'C:\\Users\\Admin\\myorf\\test_sequences\\Escherichia coli str. K-12 substr. MG1655, CDS'
sequence = load_fasta(path)
sequence = clean_sequence(sequence)
all_orf, no_stop_orfs = orf_finder(sequence=sequence,strand='+',info=True)

