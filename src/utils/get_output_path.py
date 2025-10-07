import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# File convention
# {base}.{extension}

def parser():
    parser = argparse.ArgumentParser(
        prog='CLI-find-FASTA',
        description='--out: Path to save summary.txt of the analysed FASTA',
    )
    parser.add_argument('--fasta',required=True,help='Required argument for Path-to-FASTA') #CLI user input
    parser.add_argument('--out',required=False,help='Optional output path')
    args=parser.parse_args()
    return args
args=parser()
def load_fasta(args):
    path = Path(args.fasta)
    print(path)
    while not path.exists():
        print(f"\nFile {path} not Found! Please check path to FASTA again")
        print("please enter valid path!")
        path = Path(input())
    if path.exists():
        print(f"\nReading input file from: {path}")
        print(f"\nFile suffix: {path.suffix}")
        print(f'This is the path{path}')
        try:
            for record in SeqIO.parse(path,'fasta'):
                sequence_storage = {
                     "ID":record.id,
                     "Seq":SeqRecord(Seq(record.seq)),
                     "Desc.":record.description,
                }
                yield sequence_storage
            
        except ValueError:
            print('Provided file is not a valid FASTA.')
            return None
        except Exception as e:
             print(f"Error reading {path.suffix} file: {e}")
             return None
sequence_storage = load_fasta(args)
print(type(sequence_storage))


def get_output_path(args): # User defined. njih samo prosledjujem dole da ih cekiram u sistemu.
    standard_path_fasta = Path("..")/ "results" / "ORF.fasta"
    standard_path_summary =  Path("..")/ "results" / "ORF_summary.txt"
    if args.out:
        base = Path(args.out)
        fasta_out = base.with_suffix('.fasta')
        summary_out = base.with_suffix('.txt')
        if base.suffix not in('.fasta','.fa','.fas','.txt') and not base.is_dir(): # optional case
                print(f'Warning: Provided file path suffix ({base}) may not be appropriate for further computation.\nProgram exit()')
                exit()
        if base.parent.exists() and base.is_dir(): #if user provides a dir path, default filename will be used
                print(f"Output directory detected.\nFiles will be saved as:\n{fasta_out}\n{summary_out}\n")
                fasta_out = base / "ORF.fasta"
                summary_out = base/ "ORF_summary.txt"
        if base.exists() and base.is_file():  # ovde
                print(f"\nPath exists. Do you want to overwrite it? Y/N?") 
                choice = input().strip().upper()
                while choice not in ['Y','N']:
                     print("Please choose only Y or N!")
                     choice = input().strip().upper()
                if choice == "Y":
                    print(f"Files will be saved at: {base.parent}")
                    return fasta_out,summary_out
                if choice == "N":
                    print(f"\nFile will be saved at default result directory:{standard_path_fasta.parent}\\")
                    return standard_path_fasta,standard_path_summary
        if not base.parent.exists(): #--out ..\wtf\text1.txt check if 'wtf' exists (parent)
            print(f'\nProvided file path: {base.parent} does not exist.\nDo you want me to create directory and file?\tY or N')
            input2 = str(input()).strip().upper()
            while input2 not in ("Y","N"):
                print("Please choose Y or N only.")
                input2 = str(input()).strip().upper()
            if input2 == "Y": 
                base.parent.mkdir(parents=True,exist_ok=True)
                print(f"Directory {base.parent} successfully created!")
                print(f"Files are saved inside {base.parent} directory.")
                return fasta_out,summary_out
            if input2 =="N":
                    print(f"Files will be saved at default result directory:\n{standard_path_fasta.parent}\\")
                    return standard_path_fasta,standard_path_summary
        else:
            print(f"Parent directory '{base.stem}' exists, file will be saved as: {fasta_out} & {summary_out}")
            return fasta_out,summary_out
    else:   
        print(f'\nNo output path provided.\nUsing standard result directory at {standard_path_fasta.parent}\\')
        return standard_path_fasta,standard_path_summary

#tested on CLI for majority of cases.
            