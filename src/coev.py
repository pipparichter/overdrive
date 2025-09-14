import pandas as pd 
import numpy as np 

ry2 = {
    'A':'R','G':'R',   # Purines
    'C':'Y','T':'Y',   # Pyrimidines
    'U':'Y'            # RNA pyrimidine
}

# 2. Strong (S) vs Weak (W) base pairs
# Strong = G/C (3 H-bonds), Weak = A/T (2 H-bonds)
sw2 = {
    'G':'S','C':'S',
    'A':'W','T':'W',
    'U':'W'
}

# 3. Keto (K) vs Amino (M)
# Keto: G/T(U); Amino: A/C
km2 = {
    'G':'K','T':'K','U':'K',
    'A':'M','C':'M'
}

# 4. 3-letter reduced alphabet (IUPAC-style groupings)
# (a) Purines, Pyrimidines, and G/T distinction
ry3 = {
    'A':'R','G':'R',     # Purines
    'C':'Y',             # Pyrimidine C
    'T':'T','U':'T'      # Pyrimidine T/U
}

# (b) Transition vs Transversion scheme
# Group 1: A/G (transition), Group 2: C/T (transition), Group 3: everything else (N)
transition3 = {
    'A':'R','G':'R',
    'C':'Y','T':'Y','U':'Y',
    'N':'N'
}

# 5. Binary encoding by hydrogen bonds
# 1 = strong (G/C), 0 = weak (A/T/U)
hb2 = {
    'G':'1','C':'1',
    'A':'0','T':'0','U':'0'
}

dayhoff6 = {
    'A':'A','G':'A','P':'A','S':'A','T':'A',                # Small
    'D':'B','E':'B','N':'B','Q':'B',                        # Acidic / amide
    'R':'C','H':'C','K':'C',                                # Basic
    'M':'D','I':'D','L':'D','V':'D',                        # Hydrophobic
    'F':'E','W':'E','Y':'E',                                # Aromatic
    'C':'F'                                                 # Cysteine
}

murphy10 = {
    'A':'A','G':'A','P':'A','S':'A','T':'A',                # Tiny/polar
    'D':'B','E':'B',                                        # Acidic
    'N':'C','Q':'C',                                        # Amide
    'R':'D','H':'D','K':'D',                                # Basic
    'M':'E','I':'E','L':'E','V':'E',                        # Aliphatic
    'F':'F','W':'F','Y':'F',                                # Aromatic
    'C':'G'                                                 # Cysteine
}

li5 = {
    'A':'A','G':'A','V':'A','L':'A','I':'A','P':'A',        # Hydrophobic
    'F':'B','Y':'B','W':'B',                                # Aromatic
    'D':'C','E':'C',                                        # Acidic
    'K':'D','R':'D','H':'D',                                # Basic
    'S':'E','T':'E','N':'E','Q':'E','M':'E','C':'E'         # Polar
}

hp2 = {
    'A':'H','V':'H','L':'H','I':'H','P':'H','F':'H',
    'W':'H','M':'H','Y':'H',                                # Hydrophobic
    'G':'P','S':'P','T':'P','C':'P','N':'P','Q':'P',
    'D':'P','E':'P','K':'P','R':'P','H':'P'                 # Polar
}

zhou7 = {
    'A':'1','G':'1','V':'1',                                # Small hydrophobic
    'I':'2','L':'2','M':'2',                                # Large hydrophobic
    'F':'3','W':'3','Y':'3',                                # Aromatic
    'C':'4',                                                # Cysteine
    'P':'5',                                                # Proline
    'T':'6','S':'6','N':'6','Q':'6',                        # Polar uncharged
    'D':'7','E':'7','K':'7','R':'7','H':'7'                 # Charged
}


def reduce(seq:str, alphabet:dict):
    return ''.join(alphabet.get(residue, 'X') for residue in seq)


def check_f(f:list):

    for i, f_ in enumerate(f):
        total_probability = sum(list(f_.values()))
        assert np.isclose(total_probability, 1), f'check_f: Expected total probability to be close to 1, but got {total_probability} for the position {i} distribution.'

def mutual_information(paired_msa_df:pd.DataFrame):

    length_a, length_b = len(paired_msa_df.seq_a.iloc[0]), len(paired_msa_df.seq_b.iloc[0])


    # f_ab = np.empty(shape=(length_a, length_b))
    
    # First compute the frequency of each amino acid in the alphabet.
    
    m_a = np.array([list(seq) for seq in paired_msa_df.seq_a])
    m_b = np.array([list(seq) for seq in paired_msa_df.seq_b])

    f_a = [{token.item():(col == token).mean() for token in np.unique(col)} for col in m_a.T]
    f_b = [{token.item():(col == token).mean() for token in np.unique(col)} for col in m_b.T]
    check_f(f_a)
    check_f(f_b)

    # Then compute the co-ocurrence for each pair of tokens at each position. 

    # f_ab = [[None] * length_b] * length_a
    f_ab = [[None for _ in range(length_b)] for _ in range(length_a)]

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
            h[i, j] = (h_i(i, f_a) + h_i(j, f_b) - h_ij(i, j)) / h_ij(i, j)
    
    return h
