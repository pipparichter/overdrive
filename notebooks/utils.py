import os 
import pandas as pd 
from src.files import GBKFile 
from src.files import FASTAFile, get_seq_from_fasta
from src.files import HMMerFile
from tqdm import tqdm 
import glob
import math 
import numpy as np 
import networkx as nx
import io
import subprocess
import re

# subprocess.run('export PATH=${HOME}/edirect:${PATH}', shell=True, check=True)
os.environ['PATH'] += os.pathsep + '/home/prichter/edirect'

data_dir = '../data/data-1'

def read_csv_chunk(path:str, start_line:int=309, end_line:int=513, sep:str=r'\s+', names:list=None):

    with open(path, 'r') as f:
        lines = f.readlines()
        lines = lines[start_line:end_line]
        content = ''.join(lines)
    return pd.read_csv(io.StringIO(content), sep=sep, names=names, header=None)

def get_graph_colors(graph, df:pd.DataFrame, field:str=None, palette:dict=None):
    values = df.loc[list(graph.nodes)][field] # Get the field values for each node. 
    # colors = colors.map(palette)
    colors = [palette[value] for value in values]
    return colors

def get_alignment_graph(align_df:pd.DataFrame, order:list=None):
    graph = nx.Graph()
    # Allow custom ordering of graph nodes. 
    nodes = align_df['query'].unique() if (order is None) else order

    for id_ in nodes:
        graph.add_node(id_)
    for row in align_df.itertuples():
        graph.add_edge(row.query, row.target, weight=1/row.bits)
    
    return graph

def get_layout(graph, r:int=5):

    subgraphs = [graph.subgraph(subgraph).copy() for subgraph in nx.connected_components(graph)]
    n = len(subgraphs)
    pos = nx.spring_layout(graph, weight='weight', scale=5)
    # pos = dict()
    dx, dy = 0, 0
    # r = np.linspace(0, 1, n)
    for i, subgraph in enumerate(subgraphs):
        # Get the position of each subgraph individually
        # subgraph_pos = nx.spring_layout(subgraph, k=0.5)
        subgraph_pos = {i:pos[i] for i in subgraph}
        angle = 2 * math.pi * (i / n)
        # dx, dy = r[i] * math.cos(angle), r[i] * math.sin(angle)
        dx, dy = r * math.cos(angle), r * math.sin(angle)
        
        # Shift x coordinates so the subgraphs don't overlap
        subgraph_pos = {n:(x + dx, y + dy) for n, (x, y) in subgraph_pos.items()}
        # dx += 1
        pos.update(subgraph_pos)

    return pos 


# def build_virc_dataset(product:str='VirC1', overwrite:bool=False, gbk_dir:str='../data/data-1/gbk/assemblies', plasmid_contig_ids:list=None):

#     dataset_path = f'{data_dir}/{product.lower()}.csv'
#     if os.path.exists(dataset_path) and (not overwrite):
#         return pd.read_csv(dataset_path, index_col=0)

#     dataset_df = list()
#     for path in tqdm(glob.glob(f'{gbk_dir}/*'), desc='build_virc_dataset'):
#         strain_id = os.path.basename(path).replace('.gbk', '')
#         df = GBKFile(path).to_df()
#         df = df[df['product'].str.contains(product, case=False)].copy()
#         for row in df.itertuples():
#             row_ = dict()
#             row_['contig_id'] = row.contig_id
#             row_['strain_id'] = strain_id 
#             row_['seq'] = row.seq 
#             row_['product'] = row.product
#             row_['locus_tag'] = row.locus_tag
#             row_['gene'] = row.gene
#             dataset_df.append(row_)

#     dataset_df = pd.DataFrame(dataset_df)

#     # There is VirC annotated elsewhere in a handful of the assemblies, but will probably just focus on the VirC found on the plasmid. 
#     if plasmid_contig_ids is not None:
#         mask = dataset_df.contig_id.isin(plasmid_contig_ids)
#         print(f'build_virc_dataset: Removing {(~mask).sum()} {product} proteins not found on the plasmids.')
#         dataset_df = dataset_df[mask].copy()

#     dataset_df.to_csv(dataset_path)
#     return dataset_df





def load_weisberg_2020_metadata():
    metadata_df = pd.read_csv('../data/data-1/weisberg_2020_s1.csv', skiprows=1, comment='#')
    cols = {'Strain.ID':'strain_id', '"Taxonomic" classification*':'taxonomic_classification', '3-letter code*':'biovar', 'Genomospecies':'genomospecies', 'Type  of oncogenic plasmid':'plasmid_type', 'Class of oncogenic plasmid':'plasmid_class'}
    metadata_df = metadata_df[list(cols.keys())].rename(columns=cols)

    # Clean up the strain IDs so they agree with the strain names derived from the plasmid files. 
    strain_id_map = {'Bo542 (pTiBo542)':'Bo542', 'MAFF301001 (pTiSAKURA)':'MAFF301001', 'NCPPB 3554 (K309)':'NCPPB 3554', 'MAFF03-01724 (pRi1724)':'MAFF03-01724', 'B6 (pTiB6)':'Atu B6', 'pTiOctopine':'Octopine'}
    metadata_df['strain_id'] = metadata_df.strain_id.str.strip() # Make sure there's no trailing whitespace.
    metadata_df['strain_id'] = metadata_df.strain_id.replace(strain_id_map)
    metadata_df['strain_id'] = metadata_df.strain_id.str.replace(' ', '_').str.replace('/', '_')
    metadata_df = metadata_df.fillna('none')
    metadata_df['plasmid_type'] = [f'{row.plasmid_type} ({row.plasmid_class})' for row in metadata_df.itertuples()]

    # Create plasmid IDs to match to the plasmid files, this handles the Di1411 strain, which has both an Ri and Ti plasmid.
    metadata_df['source_id'] = [f'{row.strain_id}_{row.plasmid_class.lower()}' for row in metadata_df.itertuples()]
    metadata_df = metadata_df.drop_duplicates('source_id')
    return metadata_df