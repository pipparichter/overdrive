import os 
import io 
import pandas as pd
import numpy as np 
import re 
import xml

# def read_csv_chunk(path:str, start_line:int=309, end_line:int=513, sep:str=r'\s+', names:list=None):

#     with open(path, 'r') as f:
#         lines = f.readlines()
#         lines = lines[start_line:end_line]
#         content = ''.join(lines)
#     return pd.read_csv(io.StringIO(content), sep=sep, names=names, header=None)


def parse_motif_diagram(motif_diagram:str, motif_length:int=10):

    matches = re.findall(r'\[[+-][\d]+\]', motif_diagram)
    matches = [(match_[1], match_[2:len(match_) - 1]) for match_ in matches]
    spacers = re.split(r'[-]*\[[+-][\d]+\][-]*', motif_diagram)
    spacers = [0 if (len(spacer) == 0) else spacer for spacer in spacers]
    if motif_diagram[0] == '[':
        spacers = [0] + spacers
    assert len(matches) <= len(spacers), 'parse_motif_diagram: There should be at least one spacer per match.'
    
    parsed_motif_diagram = list()
    loc = 0
    for i, (strand, motif_id) in enumerate(matches):
        row = dict()
        row['start'] = loc + int(spacers[i])
        row['stop'] = loc + int(spacers[i]) + motif_length
        row['strand'] = strand 
        row['motif_id'] = motif_id
        loc += int(spacers[i]) + motif_length
        parsed_motif_diagram.append(row)
    return parsed_motif_diagram


 
class MEMEFile():

    def __init__(self, path:str):

        pass 


class MASTFile():

    def __init__(self, path:str, motif_length:int=None):

        with open(path, 'r') as f:
            lines = f.readlines()

        section_1_start = np.where(['SECTION I' in line for line in lines])[0][0]
        section_2_start = np.where(['SECTION II' in line for line in lines])[0][0]
        section_3_start = np.where(['SECTION III' in line for line in lines])[0][0]

        section_2 = '\n'.join(lines[section_2_start:section_3_start])
        section_2 = re.split(r'[*]{2,}', section_2, flags=re.MULTILINE)[2] # The third thing in this list is the actual data. 
        section_2 = '\n'.join(section_2.split('\n')[7:]) # Remove the first several lines, which are the headers.
        section_2_df = pd.read_csv(io.StringIO(section_2), names=['id', 'e_value', 'motif_diagram'], sep=r'\s+', header=None)
        section_2_df['id'] = section_2_df['id'].str.strip()
        self.section_2_df = section_2_df
        self.motif_length = motif_length

    
    def to_df(self, max_e_value:float=10):
        # max_e_value = self.section_2_df.e_value.mean() 
        df = list()

        for row in self.section_2_df.itertuples():
            motif_diagram = row.motif_diagram
            row = row._asdict()
            for match_ in parse_motif_diagram(motif_diagram, motif_length=self.motif_length):
                # print(match_)
                row.update(match_)
                df.append(row.copy())
        df = pd.DataFrame(df)
        df['length'] = int(self.motif_length)
        df = df[df.e_value < max_e_value].copy()
        return df
    
class MASTXMLFile():

    def __init__(self, path:str):

        with open(path, 'r') as f:
            content = f.read()

        tree = xml.etree.ElementTree.fromstring(content)

        motifs = tree.findall('.//motifs/motif')
        self.motif_names = [motif.attrib['alt'] for motif in motifs]
        self.motif_lengths = [int(motif.attrib['length']) for motif in motifs]

        df = list()
        seqs = tree.findall('.//sequence')
        for seq in seqs:
            hits = seq.findall('./seg/hit')
            for hit in hits:
                row = {'id':seq.attrib['name']}
                row['motif_idx'] = int(hit.attrib['idx'])
                row['start'] = int(hit.attrib['pos'])
                row['stop'] = row['start'] + self.motif_lengths[row['motif_idx']]
                row['length'] = self.motif_lengths[row['motif_idx']]
                row['p_value'] = float(hit.attrib['pvalue'])
                row['strand'] = '+' if (hit.attrib['rc'] == 'n') else '-'
                row['motif_name'] = self.motif_names[row['motif_idx']]
                df.append(row)

        df = pd.DataFrame(df)
        self.df = self.get_significance(df)

    def get_significance(self, df:pd.DataFrame):

        idx_to_significance_map = dict()
        for motif, df_ in df.groupby('motif_name'):
            best_p_value, worst_p_value = df_.p_value.min(), df_.p_value.max()
            delta = np.abs(best_p_value - worst_p_value)
            f = lambda p_value : 1 if (delta == 0) else min(abs(worst_p_value + delta/100 - p_value) / delta, 1)
            idx_to_significance_map.update(df_.p_value.apply(f).to_dict())
        df['significance'] = df.index.map(idx_to_significance_map)
        return df

    def to_df(self):
        return self.df


# Left over from when the IDs were too long and were being truncated by MEME. 

    # def correct_ids(self, correct_ids:list):
    #     # correct_ids = np.array(correct_ids)[~np.isin(correct_ids, self.section_2_df['id'])]
    #     df = self.section_2_df[~self.section_2_df['id'].isin(correct_ids)]
    #     print(f'MASTFile.correct_ids: {len(df)} IDs to correct.')
    #     corrected_ids = dict()

    #     for id_ in df['id']:
    #         matches = [re.search(id_.strip(), correct_id) for correct_id in correct_ids]
    #         idxs = [i for i, match_ in enumerate(matches) if (match_ is not None)]
    #         if len(idxs) == 0:
    #             print(f'MASTFile.correct_ids: No match for ID {id_}.')
    #         else:
    #             corrected_ids[id_] = correct_ids[idxs[0]]
    #     self.section_2_df['id'] = np.where(self.section_2_df['id'].isin(corrected_ids), self.section_2_df['id'].map(corrected_ids), self.section_2_df['id'])
        
    #     n_incorrect = (~self.section_2_df['id'].isin(correct_ids)).sum()
    #     assert n_incorrect == 0, f'MASTFile.correct_ids: {n_incorrect} incorrect IDs remaining.'


# def get_motif_data(df:pd.DataFrame, pattern:str=motif_1a_10):
    
#     motif_df = list()
#     # for row in tqdm(list(df.itertuples()), desc='get_motif_data'):
#     for row in tqdm(list(df.itertuples()), desc='get_motif_data'):

#         seq = {'+':row.seq, '-':str(Seq(row.seq).reverse_complement())}
#         matches = {strand:re.finditer(pattern, s) for strand, s in seq.items()}

#         for strand, matches_ in matches.items():
#             for match in matches_:
#                 # print(matches)
#                 row_ = {'id':row.Index}
#                 row_['strand'] = strand
#                 if strand == '-':
#                     row_['stop'] = len(seq[strand]) - match.start()
#                     row_['start'] = len(seq[strand]) - match.end()
#                 else:
#                     row_['start'] = match.start()
#                     row_['stop'] = match.end()
#                 row_['match'] = match.group(0)
#                 row_['seq'] = seq[strand]
#                 motif_df.append(row_)
#     motif_df = pd.DataFrame(motif_df)
#     return motif_df



