import os 
import pandas as pd 
import numpy as np 
from src.files import FASTAFile, HMMerFile, get_seq_from_fasta
import glob 
import re 
import json

has_right_border = lambda df : np.any(df.annotation.str.contains('right_border'))
has_left_border = lambda df : np.any(df.annotation.str.contains('left_border'))
has_t_dna = lambda df : df.target_name.map(df.groupby('target_name').apply(lambda df_ : has_right_border(df_) and has_left_border(df_), include_groups=False))
has_overlap = lambda start_a, stop_a, start_b, stop_b : not ((start_a > stop_b) or (stop_a < start_b))


def load_hmmer_files(hmmer_dir:str='../data/data-1/hmmer', max_e_value:float=0.05):
    hmmer_dir = os.path.join(hmmer_dir, '*')
    hmmer_df = [HMMerFile(path).to_df() for path in glob.glob(hmmer_dir)]
    hmmer_df = pd.concat([df for df in hmmer_df if (len(df) > 0)])
    hmmer_df = hmmer_df.rename(columns={'query_name':'annotation'})
    hmmer_df = hmmer_df[hmmer_df.e_value < max_e_value].copy() # Get only significant hits. 
    # hmmer_df = hmmer_df[hmmer_df.target_description.str.contains('agro|Agro|Rhizo|rhizo', regex=True)].copy()

    for annotation, df in hmmer_df.groupby('annotation'):
        print(f'load_hmmer: Num. hits for query {annotation}:', len(df))

    # I think we probably need to adjust the starts to be zero-indexed, as with NCBI files.
    hmmer_df['start'] = np.min([hmmer_df.target_to.values, hmmer_df.target_from.values], axis=0) - 1
    hmmer_df['stop'] = np.max([hmmer_df.target_to.values, hmmer_df.target_from.values], axis=0)

    return hmmer_df


def _resolve_borders(target_df:pd.DataFrame):

    assert np.all(target_df.target_name == target_df.target_name.iloc[0]), 'get_t_dna_borders: Expected the DataFrame to correspond to only one target sequence.'
    target_df = target_df.sort_values('e_value', ascending=True)
    borders = dict()
    for row in target_df.itertuples():
        found_overlap = False
        for start, stop in borders.keys():
            if has_overlap(start, stop, row.start, row.stop):
                borders[(start, stop)].append({'e_value':row.e_value, 'strand':row.strand, 'annotation':row.annotation, 'source_id':row.source_id})
                found_overlap = True
                break
        if not found_overlap:
            borders[(row.start, row.stop)] = [{'e_value':row.e_value, 'strand':row.strand, 'annotation':row.annotation, 'source_id':row.source_id}]

    return borders


def _filter_borders(target_df:pd.DataFrame):

    get_next = lambda i : 'none' if (i == (len(target_df) - 1)) else target_df.iloc[i + 1].annotation
    get_previous = lambda i : 'none' if (i == 0) else target_df.iloc[i - 1].annotation
    
    target_df = target_df.sort_values('start')
    target_df['previous_annotation'] = [get_previous(i) for i in range(len(target_df))]
    target_df['next_annotation'] = [get_next(i) for i in range(len(target_df))]

    mask = (target_df.annotation.str.contains('right_border'))
    mask = mask & target_df.next_annotation.str.contains('none|left_border')
    mask = mask & target_df.previous_annotation.str.contains('none|left_border')
    return target_df[mask].copy()


def resolve_borders(hmmer_df:pd.DataFrame, filter_:bool=True):

    hmmer_df = hmmer_df[hmmer_df.annotation != 'overdrive'].copy()

    borders = {target_name:_resolve_borders(df) for target_name, df in hmmer_df.groupby('target_name')}

    hmmer_df = list()
    for target_name, loci in borders.items():
        for (start, stop), annotations in loci.items():
            assert len(annotations) <= 2, 'Expected no more than two overlapping border annotations.'
            row = {'start':start, 'stop':stop}
            row['target_name'] = target_name 
            row.update(annotations[0])
            if len(annotations) > 1:
                row['overlapping_e_value'] = annotations[-1]['e_value']
            hmmer_df.append(row)
    hmmer_df = pd.DataFrame(hmmer_df)

    if filter_:
        hmmer_df = pd.concat([_filter_borders(df).assign(target_name=target_name) for target_name, df in hmmer_df.groupby('target_name')])

    return hmmer_df # , borders




def build_overdrive_dataset(overwrite:bool=False, fasta_dir:str='../data/data-1/ncbi/fasta/', hmmer_df:pd.DataFrame=None, length:int=200, data_dir:str='../data/data-1/', min_length:int=50):

    hmmer_df = resolve_borders(hmmer_df)
    hmmer_df.to_csv(os.path.join(data_dir, 'borders.csv'))
    hmmer_df = hmmer_df[hmmer_df.annotation.str.contains('right_border')].copy() # Only care about right border hits for grabbing the coordinates.
    print(f'build_overdrive_dataset: Found {len(hmmer_df)} right borders.')
    
    dataset_path = os.path.join(data_dir, 'overdrive.csv')
    if os.path.exists(dataset_path) and (not overwrite):
        return pd.read_csv(dataset_path, index_col=0)

    dataset_df = list()
    fasta_paths = glob.glob(os.path.join(fasta_dir, '**/*.fn'), recursive=True)
    print(f'build_overdrive_dataset: Found {len(fasta_paths)} FASTA files.')

    for source_id, df in hmmer_df.groupby('source_id'): # The IDs in the HMMer output are the source file names.
        fasta_path = [path for path in fasta_paths if (re.search(source_id, path) is not None)][0]
        fasta_df = FASTAFile(fasta_path).to_df()

        for row in df.itertuples(): # There can be multiple annotated right borders. 

            # These coordinates are relative to the forward strand, but target_to and target_from are reversed if the hit
            # is on the reverse strand. 
            start = row.stop if (row.strand == '+') else row.start # Get the end of the right border hit. 
            stop = start + length if (row.strand == '+') else start - length # This makes sense. 
            start, stop = min(start, stop), max(start, stop)
            # next_stop = stop + length if (row.strand == '+') else stop - length # This makes sense. 
            
            left_of_right_border_stop = row.start if (row.strand == '+') else row.stop 
            left_of_right_border_start = left_of_right_border_stop - 100 if (row.strand == '+') else left_of_right_border_stop + 100
            left_of_right_border_stop, left_of_right_border_start = max(left_of_right_border_stop, left_of_right_border_start), min(left_of_right_border_stop, left_of_right_border_start)
            # next_start, next_stop = min(stop, next_stop), max(stop, next_stop)

            row_ = dict()
            # The get_seq_from_fasta function handles the reverse complement. 
            row_['seq'] = get_seq_from_fasta(fasta_df, id_=row.target_name, coords=(start, stop), strand=row.strand)
            # row_['next_seq'] = get_seq_from_fasta(fasta_df, id_=row.target_name, coords=(next_start, next_stop), strand=row.strand)
            # Also collect the right border sequence as a sanity check. 
            row_['right_border_seq'] = get_seq_from_fasta(fasta_df, id_=row.target_name, coords=(row.start, row.stop), strand=row.strand)
            row_['left_of_right_border_seq'] = get_seq_from_fasta(fasta_df, id_=row.target_name, coords=(left_of_right_border_start, left_of_right_border_stop), strand=row.strand)
            row_['left_of_right_border_start'] = left_of_right_border_start
            row_['left_of_right_border_stop'] = left_of_right_border_stop
            row_['start'] = start
            row_['stop'] = stop
            row_['strand'] = row.strand 
            row_['source_id'] = source_id
            row_['contig_id'] = row.target_name
            dataset_df.append(row_)

    dataset_df = pd.DataFrame(dataset_df) 
    dataset_df = dataset_df[dataset_df.seq.apply(len) > min_length].copy()
    dataset_df['id'] = [f'{row.contig_id}:{row.start}-{row.stop}' for row in dataset_df.itertuples()] # Make helpful IDs for each overdrive. 
    dataset_df.set_index('id').to_csv(dataset_path)
    return dataset_df.set_index('id')