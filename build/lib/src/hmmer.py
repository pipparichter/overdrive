import pandas as pd 
import numpy as np 
import os 

# http://eddylab.org/software/hmmer/Userguide.pdf

class HMMerFile():

    # The "from" and "to" fields encode strand information, i.e. if from > to, then the hit is on the opposite strand. 
    # I think env corresponds to "envelope," but not completely sure what these are. 
    # See page 68 in user guide for full descriptions of these fields. "query_name" is populated with the name of the HMM. 
    fields = ['target_name', 'target_accession',  'query_name', 'query_accession', 'query_from', 'query_to', 'target_from', 'target_to', 'env_from', 'env_to', 'length', 'strand', 'e_value' , 'score', 'bias', 'target_description']

    def __init__(self, path:str):
        self.plasmid_id = os.path.basename(path).replace('.tsv', '')
        self.df = pd.read_csv(path, comment='#', sep=r'\s+', header=None, names=self.fields)

    def to_df(self, max_e_value:float=0.05):
        df = self.df.copy()
        df['plasmid_id'] = self.plasmid_id
        df = df[df.e_value < max_e_value].copy() # Remove weak hits. 
        return df
