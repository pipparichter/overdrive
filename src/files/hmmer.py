import pandas as pd 
import numpy as np 
import os 
import re
import io
import glob

# http://eddylab.org/software/hmmer/Userguide.pdf

# NZ_AFDM01000009.1    -          T-DNA_right_border   -                1      15  762824  762810  762824  762808 1626953    -         5.5   10.0   1.4  Acinetobacter baumannii OIFC189 Contig1221746139, whole genome shotgun sequence

class HMMerFile():

    # The "from" and "to" fields encode strand information, i.e. if from > to, then the hit is on the opposite strand. 
    # I think env corresponds to "envelope," but not completely sure what these are. 
    # See page 68 in user guide for full descriptions of these fields. "query_name" is populated with the name of the HMM. 
    fields = ['target_name', 'target_accession',  'query_name', 'query_accession', 'query_from', 'query_to', 'target_from', 'target_to', 'env_from', 'env_to', 'length', 'strand', 'e_value' , 'score', 'bias', 'target_description']

    @staticmethod
    def parse_line(line:str): # Need to manually parse because of inconsistent output format.
        n = len(HMMerFile.fields)
        # pattern = '^' + ''.join([r'(^\s)\s+'] * (n - 1)) + '(.+)$'
        pattern = ''.join([r'(\S+)\s+'] * (n - 1)) + '(.+)'
        line = re.search(pattern, line, flags=re.DOTALL)
        line = [line.group(i + 1) for i in range(n)]
        return line

    def __init__(self, path:str):
        self.id_ = os.path.basename(path).replace('.tsv', '')
        
        with open(path, 'r') as f:
            self.lines = [line for line in f.readlines() if (line[0] != '#')]
        # self.df = pd.read_csv(path, comment='#', sep=r'\s{2,}', header=None, names=self.fields, engine='python')

    def to_df(self, max_e_value:float=0.05):
        lines = [HMMerFile.parse_line(line) for line in self.lines]
        content = '\n'.join(['\t'.join(line) for line in lines])
        df = pd.read_csv(io.StringIO(content), sep='\t', header=None, names=HMMerFile.fields)
        df['id'] = self.id_
        # df = df[df.e_value < max_e_value].copy() # Remove weak hits. 
        return df
    

def load_hmmer(hmmer_dir:str='../data/data-1/hmmer', max_e_value:float=0.05):
    hmmer_dir = os.path.join(hmmer_dir, '*')
    hmmer_df = [HMMerFile(path).to_df() for path in glob.glob(hmmer_dir)]
    hmmer_df = pd.concat([df for df in hmmer_df if (len(df) > 0)])
    hmmer_df = hmmer_df[hmmer_df.e_value < max_e_value].copy() # Get only significant hits. 
    # hmmer_df = hmmer_df[hmmer_df.target_description.str.contains('agro|Agro|Rhizo|rhizo', regex=True)].copy()

    for query_name, df in hmmer_df.groupby('query_name'):
        print(f'load_hmmer: Num. hits for query {query_name}:', len(df))

    return hmmer_df
