import os 
import io 
import pandas as pd
import numpy as np 
import re 
import xml
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch

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

        with open(path, 'r') as f:
            self.content = f.read()

        self.motif_names = np.unique(re.findall(r'MEME-\d+', self.content))
        self.motif_names = [str(motif_name) for motif_name in self.motif_names]

    def get_regular_expression(self, motif_names:str=None):
        motif_names = self.motif_names if (motif_names is None) else motif_names

        lines = self.content.split('\n')
        
        regex = dict()
        for motif_name in motif_names:
            is_regex_section_start = lambda line : re.search(f'{motif_name} regular expression', line) is not None
            regex_section_start_idx = np.where(np.array([is_regex_section_start(line) for line in lines]))[0][0]
            regex_line_idx = regex_section_start_idx + 2
            regex[motif_name] = lines[regex_line_idx].strip()
        return regex


# <letter id="A" symbol="A" complement="T" name="Adenine" colour="CC0000"/>
# <letter id="C" symbol="C" complement="G" name="Cytosine" colour="0000CC"/>
# <letter id="G" symbol="G" complement="C" name="Guanine" colour="FFB300"/>
# <letter id="T" symbol="T" aliases="U" complement="A" name="Thymine" colour="008000"/>
# <letter id="N" symbol="N" aliases="X." equals="ACGT" name="Any base"/>
# <letter id="V" symbol="V" equals="ACG" name="Not T"/>
# <letter id="H" symbol="H" equals="ACT" name="Not G"/>
# <letter id="D" symbol="D" equals="AGT" name="Not C"/>
# <letter id="B" symbol="B" equals="CGT" name="Not A"/>
# <letter id="M" symbol="M" equals="AC" name="Amino"/>
# <letter id="R" symbol="R" equals="AG" name="Purine"/>
# <letter id="W" symbol="W" equals="AT" name="Weak"/>
# <letter id="S" symbol="S" equals="CG" name="Strong"/>
# <letter id="Y" symbol="Y" equals="CT" name="Pyrimidine"/>
# <letter id="K" symbol="K" equals="GT" name="Keto"/>


    
class MASTXMLFile():

    regex_map = dict()
    regex_map['N'] = '[ACGT]'
    regex_map['V'] = '[ACG]' # Not T
    regex_map['H'] = '[ACT]' # Not G
    regex_map['D'] = '[AGT]' # Not C
    regex_map['B'] = '[CGT]' # Not A
    regex_map['M'] = '[AC]' # Amino 
    regex_map['R'] = '[AG]' # Purine
    regex_map['W'] = '[AT]' # Weak
    regex_map['S'] = '[CG]' # Strong
    regex_map['Y'] = '[CT]' # Pyramidine
    regex_map['K'] = '[GT]' # Keto

    description_map = dict()
    description_map['N'] = '[ACGT]'
    description_map['V'] = '[ACG]'
    description_map['H'] = '[ACT]' # Not G
    description_map['D'] = '[AGT]' # Not C
    description_map['B'] = '[CGT]' # Not A
    description_map['M'] = 'amino' # Amino 
    description_map['R'] = 'purine' # Purine
    description_map['W'] = 'weak' # Weak
    description_map['S'] = 'strong' # Strong
    description_map['Y'] = 'pyramidine' # Pyramidine
    description_map['K'] = 'keto' # Keto

    def __init__(self, path:str):

        with open(path, 'r') as f:
            content = f.read()

        tree = xml.etree.ElementTree.fromstring(content)

        motifs = tree.findall('.//motifs/motif')
        self.motif_names = [motif.attrib['alt'] for motif in motifs]
        self.motif_lengths = [int(motif.attrib['length']) for motif in motifs]
        self.motif_codes = [motif.attrib['id'] for motif in motifs]

        df = list()
        seqs = tree.findall('.//sequence')
        # print(f'MASTXMLFile__init__: Found entries for {len(seqs)} sequences.')
        for seq in seqs:
            for seg in seq.findall('./seg'):
                seg_start = int(seg.attrib['start'])
                seg_seq = seg.find('./data').text.replace('\n', '').replace('\t', '')
                for hit in seg.findall('./hit'):
                    row = {'id':seq.attrib['name']}
                    row['motif_idx'] = int(hit.attrib['idx'])
                    row['start'] = int(hit.attrib['pos'])
                    row['stop'] = row['start'] + self.motif_lengths[row['motif_idx']]
                    row['motif_length'] = self.motif_lengths[row['motif_idx']]
                    row['seq_length'] = int(seq.attrib['length'])
                    row['p_value'] = float(hit.attrib['pvalue'])
                    row['strand'] = '+' if (hit.attrib['rc'] == 'n') else '-'
                    row['motif_name'] = self.motif_names[row['motif_idx']]

                    # Get motif stop and start indices relative to the sequence segment. 
                    start = row['start'] - seg_start
                    stop = row['stop'] - seg_start
                    row['seq'] = seg_seq[start:stop]
                    row['seq_upstream'] = seg_seq[max(0, start - 10):start] 
                    row['seq_downstream'] = seg_seq[stop:min(len(seg_seq), stop + 10)]
                    df.append(row)

        df = pd.DataFrame(df)
        self.df = df
        # self.df = self.get_significance(df)

    # def get_significance(self, df:pd.DataFrame):

    #     idx_to_significance_map = dict()
    #     for motif_name, df_ in df.groupby('motif_name'):
    #         best_p_value, worst_p_value = df_.p_value.min(), df_.p_value.max()
    #         delta = np.abs(best_p_value - worst_p_value)
    #         f = lambda p_value : 1 if (delta == 0) else min(abs(worst_p_value + delta/100 - p_value) / delta, 1)
    #         idx_to_significance_map.update(df_.p_value.apply(f).to_dict())
    #     df['significance'] = df.index.map(idx_to_significance_map)
    #     return df

    def to_df(self):
        return self.df
    
    def __len__(self):
        return self.df['id'].nunique()
    
    def get_regular_expression(self):
        regex = dict()
        for motif_name, motif_code in zip(self.motif_names, self.motif_codes):
            motif_code = list(motif_code)
            motif_regex = ''.join([MASTXMLFile.regex_map.get(symbol, symbol) for symbol in motif_code])
            regex[motif_name] = motif_regex
        return regex
    
    def get_description(self):
        description = dict()
        for motif_name, motif_code in zip(self.motif_names, self.motif_codes):
            motif_code = list(motif_code)
            motif_description = ','.join([f'[{i}] {MASTXMLFile.description_map.get(symbol, symbol)}' for i, symbol in enumerate(motif_code)])
            description[motif_name] = motif_description
        return description
    

    def plot(self, print_regex:bool=False, max_p_value:float=None):

        cmap = mpl.colormaps.get_cmap('tab20c')
        cmap = cmap.resampled(len(self.motif_names))
        palette = {f:cmap(i) for i, f in enumerate(self.motif_names)} # Map each category to a color. 

        figure_df = self.df.copy()
        if max_p_value is not None:
            figure_df = figure_df[figure_df.p_value < max_p_value].copy()

        fig, ax = plt.subplots(figsize=(10, max(1, len(self)//6))) 

        handles, labels = zip(*[(Patch(facecolor=color), label) for label, color in palette.items()])
        ax.legend(handles=handles, labels=labels, loc='upper right')

        y = 0
        y_labels = list()

        for id_, df in figure_df.groupby('id', sort=False):
            x_min, x_max = 0, figure_df.seq_length.iloc[0]
            ax.hlines(y, xmax=x_max, xmin=x_min, color='black', lw=0.5)

            for motif_name, df_ in df.groupby('motif_name'):

                for row in df_.itertuples():
                    dx = row.motif_length if (row.strand == '+') else -row.motif_length 
                    x = row.start if (row.strand == '+') else row.stop 

                    lw = (row.significance * 3) if hasattr(row, 'significance') else 1.5
                    lw = 1.5 if np.isnan(lw) else lw
                    lw = 1.5
                    ax.arrow(x, y, dx, 0, length_includes_head=True, color=palette[motif_name], lw=lw, head_width=0.8, alpha=1)
                    
            y_labels.append(id_)
            y += 1

        ax.set_ylim(ymax=y + 1, ymin=-1)
        ax.set_yticks(np.arange(y), labels=y_labels, fontsize='x-small')

        for _, spine in ax.spines.items():
            spine.set_visible(False)

        if print_regex:
            for motif_name, regex in self.get_regular_expression().items():
                print(f'{motif_name}: {regex}')
        # ax.set_title(str(self.get_regular_expression()))

        plt.show()


# class MASTFile():

#     def __init__(self, path:str, motif_length:int=None):

#         with open(path, 'r') as f:
#             lines = f.readlines()

#         section_1_start = np.where(['SECTION I' in line for line in lines])[0][0]
#         section_2_start = np.where(['SECTION II' in line for line in lines])[0][0]
#         section_3_start = np.where(['SECTION III' in line for line in lines])[0][0]

#         section_2 = '\n'.join(lines[section_2_start:section_3_start])
#         section_2 = re.split(r'[*]{2,}', section_2, flags=re.MULTILINE)[2] # The third thing in this list is the actual data. 
#         section_2 = '\n'.join(section_2.split('\n')[7:]) # Remove the first several lines, which are the headers.
#         section_2_df = pd.read_csv(io.StringIO(section_2), names=['id', 'e_value', 'motif_diagram'], sep=r'\s+', header=None)
#         section_2_df['id'] = section_2_df['id'].str.strip()
#         self.section_2_df = section_2_df
#         self.motif_length = motif_length

    
#     def to_df(self, max_e_value:float=10):
#         # max_e_value = self.section_2_df.e_value.mean() 
#         df = list()

#         for row in self.section_2_df.itertuples():
#             motif_diagram = row.motif_diagram
#             row = row._asdict()
#             for match_ in parse_motif_diagram(motif_diagram, motif_length=self.motif_length):
#                 # print(match_)
#                 row.update(match_)
#                 df.append(row.copy())
#         df = pd.DataFrame(df)
#         df['length'] = int(self.motif_length)
#         df = df[df.e_value < max_e_value].copy()
#         return df

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



