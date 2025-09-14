import os 
import re 
from tqdm import tqdm 
import pandas as pd 
import numpy as np 
# from Bio.Seq import Seq
# # from Bio.Data import CodonTable
# import Bio.Seq

# TODO: I should automatically detect the feature labels instead of specifying beforehand.
# TODO: Organize the fields a little better. 
# TODO: Still getting a not divisible by 3 warning when translating pseudogenes, I think because things are running off the end of a contig. 

# If a gene as at the edge of a contig, it may be marked partial, but can will still be translated. The only CDS features which are not translated
# are those marked as pseudo, which have a frameshift, internal stop, etc. There are also non-CDS features marked as pseudo, which may not be translated.

class GBKFile():

    fields = ['contig_id', 'product', 'note', 'pseudo', 'gene', 'locus_tag'] # Basic fields. 
    dtypes = {'start':int, 'stop':int, 'strand':int, 'pseudo':bool, 'continuous':bool}
    for field in fields: # Set all other fields to the string datatype.
        dtypes[field] = dtypes.get(field, str) 

    qualifier_pattern = re.compile(r'/([a-zA-Z_]+)="([^"]+)"') # Two capturing groups so that re.findall gets a list of tuples. 
    # The question mark is for lazy matching, so it only matches until the first /translation qualifier. 
    # cds_feature_pattern = r'CDS\s+.*?/translation=".*?"'
    cds_feature_pattern = r'CDS\s+.*?(?=\s{2,}(?:gene|sig_peptide|rRNA|tRNA|CDS)\s{2,}|\Z)'

    # The order of these coordinate patterns is important, as want to prioritize the outer match; coordinates can be of the form complement(join(..)),
    # for example, and I want to make sure I don't just match the substring. It is also important to enable multiline, as sometimes the coordinates span multiple lines.
    coordinate_pattern = re.compile('|'.join([r'(complement\(.+\))', r'(join\(.+\))', r'(order\(.+\))', r'([\<\d]+\.\.[\>\d]+)']), flags=re.MULTILINE)

    @staticmethod
    def parse_coordinate(coordinate:str):
        '''For information on location operators, see section 3.4.2.2 of https://www.insdc.org/submitting-standards/feature-table/#3.4.2.2'''
        # In bacteria, the joins can be used in cases like programmed frameshifts (e.g. if there is ribosomal slippage)
        parsed_coordinate = list()

        strand = 1
        if re.match(r'complement\((.+)\)', coordinate) is not None:
            strand = -1 
            coordinate = re.match(r'complement\((.+)\)', coordinate).group(1) 
        
        continuous = True
        if re.match(r'join\((.+)\)', coordinate) is not None:
            continuous = False 
            coordinate = re.match(r'join\((.+)\)', coordinate).group(1)

        if re.match(r'order\((.+)\)', coordinate) is not None:
            continuous = False 
            coordinate = re.match(r'order\((.+)\)', coordinate).group(1)  

        for range_ in coordinate.split(','): # Handles the case of a potential join. 
            try: # Edge case where coordinate looks like join(876461,1..1493)
                start, stop = range_.split('..')
            except ValueError:
                start, stop = range_, range_
            partial = ('1' if ('<' in start) else '0') + ('1' if ('>' in stop) else '0')
            start = int(start.replace('<', '').replace('>', ''))
            stop = int(stop.replace('<', '').replace('>', ''))
            parsed_coordinate.append({'start':start, 'stop':stop, 'continuous':continuous, 'strand':strand, 'partial':partial})

        return parsed_coordinate

    @staticmethod
    def parse_qualifiers(qualifiers:list, delimiter:str=';'):
        # Need to account for the fact that a single entry can have multiples of the same qualifier.
        parsed_qualifiers = dict()
        for field, value in qualifiers:
            if (field not in parsed_qualifiers):
                parsed_qualifiers[field] = value
            else:
                value = value.replace(delimiter, ' ')
                # assert (delimiter not in value), f"GBKFile.parse_qualifiers: There is already a '{delimiter}' in \"{value}\" for field {field}. Need to use a different delimiter."
                parsed_qualifiers[field] += delimiter + value
        return parsed_qualifiers

    @staticmethod
    def parse_feature(qualifiers:str) -> dict:
        # Extract the gene coordinates, which do not follow the typical field pattern. 
        coordinate = re.search(GBKFile.coordinate_pattern, qualifiers).group(0)
        pseudo = ('/pseudo' in qualifiers)
        # Remove all newlines or any more than one consecutive whitespace character.
        # This accomodates the fact that some of the fields are multi-line. 
        qualifiers = re.sub(r'[\s]{2,}|\n', ' ', qualifiers) 
        qualifiers = re.findall(GBKFile.qualifier_pattern, qualifiers) # Returns a list of matches. 
        qualifiers = GBKFile.parse_qualifiers(qualifiers) # Convert the qualifiers to a dictionary. 

        parsed_feature = list()
        for parsed_coordinate in GBKFile.parse_coordinate(coordinate): # There can be multiple coordinates when there are joins involved. 
            parsed_feature_ = dict()
            parsed_feature_['coordinate'] = coordinate
            parsed_feature_['pseudo'] = pseudo
            parsed_feature_['seq'] = re.sub(r'[\s]+', '', qualifiers.get('translation', 'none')) # Remove leftover whitespace in sequence. 
            parsed_feature_['locus_tag'] = qualifiers.get('locus_tag', 'none')
            # if parsed_feature['locus_tag'] == 'AS1F9_05626':
            parsed_feature_['product'] = qualifiers.get('product', 'none')
            parsed_feature_['gene'] = qualifiers.get('gene', 'none')
            parsed_feature_['note'] = qualifiers.get('note', 'none')
            parsed_feature_.update(parsed_coordinate)
            parsed_feature.append(parsed_feature_)

        return parsed_feature 

    @staticmethod
    def parse_contig(contig:str) -> pd.DataFrame:
        seq = re.search(r'ORIGIN(.*?)(?=(//)|$)', contig, flags=re.DOTALL).group(1)
        contig_id = contig.split()[0].strip()

        contig = contig.split('ORIGIN')[0] # Remove the nucleotide sequence. 
        contig = contig.split('FEATURES')[-1] # Remove the header information. 
        cds_features = re.findall(GBKFile.cds_feature_pattern, contig, flags=re.DOTALL)

        if len(cds_features) == 0: # Catches the case where the contig is not associated with any gene features. 
            return contig_id, seq, None

        cds_features = [feature.replace('CDS', '').strip() for feature in cds_features] # Tuples of form (feature, data). 
        df = pd.DataFrame([parsed_feature for feature in cds_features for parsed_feature in GBKFile.parse_feature(feature)])
        df['contig_id'] = contig_id

        return contig_id, seq, df 


    def __init__(self, path:str):
        
        self.path, self.df, self.contigs = path, list(), dict() 
        # self.plasmid_id = os.path.basename(path).replace('.gbk', )

        with open(path, 'r') as f:
            content = f.read()

        # If there are multiple contigs in the file, the set of features corresponding to the contig is marked by a "contig" feature.
        # I used the lookahead match here because it is not treated as a capturing group (so I don't get LOCUS twice). 
        contigs = re.findall(r'LOCUS(.*?)(?=LOCUS|$)', content, flags=re.DOTALL) # DOTALL flag means the dot character also matches newlines.
 
        for contig in contigs:
            contig_id, contig_seq, contig_df = GBKFile.parse_contig(contig)
            self.contigs[contig_id] = re.sub(r'[\n\s0-9]', '', contig_seq).upper() # Remove line numbers and newlines from the contig sequence
            if (contig_df is not None):
                self.df.append(contig_df)

        # It's important to reset the index after concatenating so every feature has a unique label for the subsequent evidence merging. 
        self.df = pd.concat(self.df).reset_index(drop=True)


    def to_fasta(self, path:str=None, fmt:str='nt', add_plasmid_id:bool=True):

        content = ''
        if fmt == 'nt':
            for contig_id, contig in self.contigs.items():
                content += f'>{contig_id}\n'
                content += f'{contig}\n'

        elif fmt == 'aa':
            for row in self.df.itertuples():
                content += f'>{row.locus_tag}|{row.gene}|{row.product}\n'
                content += row.seq + '\n'
        
        if path is None:
            return content 
        else:
            with open(path, 'w') as f:
                f.write(content)

    def to_df(self):
        df = self.df.copy()
        df['file_name'] = os.path.basename(self.path)
        return df



