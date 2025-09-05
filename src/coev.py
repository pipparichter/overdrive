import pandas as pd 
import numpy as np 

def coev_mi(paired_msa_df:pd.DataFrame):

    length_a, length_b = len(paired_msa_df.seq_a.iloc[0]), len(paired_msa_df.seq_a.iloc[0])


    # f_ab = np.empty(shape=(length_a, length_b))
    
    # First compute the frequency of each amino acid in the alphabet.
    
    m_a = np.array([list(seq) for seq in paired_msa_df.seq_a])
    m_b = np.array([list(seq) for seq in paired_msa_df.seq_b])

    f_a = [{token.item():(col == token).mean() for token in np.unique(col)} for col in m_a.T]
    f_b = [{token.item():(col == token).mean() for token in np.unique(col)} for col in m_b.T]
    # Then compute the co-ocurrence for each pair of tokens at each position. 

    f_ab = [[None] * length_b] * length_a
    for i in range(length_a):
        for j in range(length_b):
            pairs = np.array([''.join(pair) for pair in (zip(m_a.T[i], m_b.T[j]))])
            f_ab[i][j] = {pair:(pairs == pair).mean() for pair in np.unique(pairs)}

    def h_i(i, f:np.ndarray, k:int=None):
        # k = len(f[i]) # Get the alphabet size. 
        return - np.sum([p * np.log2(p) for p in f[i].values()])

    def h_ij(i, j, f_ab:np.ndarray=f_ab, f_a:np.ndarray=f_a, f_b:np.ndarray=f_b):
        # k = len(f_ab[i][j]) # Get the alphabet size. 
        return - np.sum([p * np.log2(p) for p in f_ab[i][j].values()])
    
    h = np.empty(shape=(length_a, length_b))
    for i in range(length_a):
        for j in range(length_b):
            h[i, j] = h_i(i, f_a) + h_i(j, f_b) - h_ij(i, j)
    
    return h
