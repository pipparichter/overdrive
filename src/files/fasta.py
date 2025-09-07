import os 
import re 
from typing import List, Dict, Tuple, NoReturn
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd 

# For the HMMer hits, indexing is always relative to the forward strand, so start and stop can be determined by taking the max of the coordinates. 

def get_seq_from_fasta(df:pd.DataFrame, id_:str=None, coords:tuple=None, strand:int=None):
    # When using this for extracting sequences from HMM hits, the coordinates are always relative to the forward strand, so 
    # only take the reverse complement after extracting the subsequence. 
    seq = df.loc[id_].seq 
    if coords is not None:
        # start = start - 1 # I think just do the same start adjustment we needed to do before, seems like HMMer is also 1-indexed. 
        seq = seq[coords[0]:coords[1]] 

    if (strand == '-'):
        seq = str(Seq(seq).reverse_complement())

    return seq 

def parser_default(description:str):
    return {'description':description}


class FASTAFile():
    def __init__(self, path:str=None, df:pd.DataFrame=None):
        '''Initialize a FASTAFile object.'''

        if (path is not None):
            f = open(path, 'r')
            self.seqs, self.ids, self.descriptions = [], [], []
            for record in SeqIO.parse(path, 'fasta'):
                self.ids.append(record.id)
                self.descriptions.append(record.description.replace(record.id, '').strip())
                self.seqs.append(str(record.seq))
            f.close()
            
        if (df is not None):
            self.seqs = df.seq.values
            self.ids = df.index.values 
            self.descriptions = df.description.values if ('description' in df.columns) else [''] * len(self.ids)
        
        self.seqs = [seq.replace(r'*', '') for seq in self.seqs] # Remove the terminal * character if present.

    def __len__(self):
        return len(self.seqs)
            
    def to_df(self) -> pd.DataFrame:

        parser = parser_default

        df = []
        for id_, seq, description in zip(self.ids, self.seqs, self.descriptions):
            row = parser(description)
            row['id'] = id_
            row['seq'] = seq
            df.append(row)
        df = pd.DataFrame(df).set_index('id')
        return df

    def write(self, path:str, mode:str='w') -> NoReturn:
        f = open(path, mode=mode)
        records = []
        for id_, seq, description in zip(self.ids, self.seqs, self.descriptions):
            record = SeqRecord(Seq(seq), id=id_, description=description)
            records.append(record)
        SeqIO.write(records, f, 'fasta')
        f.close()