import numpy as np
import pandas as pd

def get_nnmd(df, score_column, ess_column='essentiality'):
    """Given essentiality scores or similar for a selection of gene targets in a pooled CRISPR screen with bulk readout, null-normalized mean difference (NNMD) measures the separation between scores for essential and non-essential genes.
    It is defined as the difference in means between scores of both groups divided by the standard deviation of the non-essential genes.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the essentiality scores and essentiality annotation
    score_column : str
        Column name of the essentiality scores
    ess_column : str
        Column name of the essentiality annotation
    
    Returns
    -------
    float
        The NNMD score
    """
    tdf = df[df[ess_column].isin(['essential', 'non-essential'])]
    means = tdf.groupby(ess_column).apply(lambda x: x[score_column].mean())
    diff = means['essential'] - means['non-essential']
    std_ne = tdf[tdf[ess_column]=='non-essential'][score_column].std()
    return np.abs(diff / std_ne)

def annotate_essentiality(gene_list):
    essential = pd.read_csv("../resources/bagel_supplement/CEGv2.txt", sep='\t', index_col=0)
    non_essential = pd.read_csv("../resources/bagel_supplement/NEGv1.txt", sep='\t', index_col=0)
    assert np.intersect1d(gene_list, essential.index).size + np.intersect1d(gene_list, non_essential.index).size > 0, "Gene lists not overlapping, maybe your gene symbols are wrong?"
    return ['essential' if gene in essential.index else 'non-essential' if gene in non_essential.index else 'unknown' for gene in gene_list]

def parse_sam_file(sam_file):
    data = []
    with open(sam_file, 'r') as f:
        for line in f:
            # Skip header lines if any
            if line.startswith('@'):
                continue
            
            fields = line.strip().split('\t')
            # Extract relevant fields
            query_name = fields[0]
            flag = int(fields[1])
            reference_name = fields[2]
            position = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            rnext = fields[6]
            pnext = int(fields[7])
            tlen = int(fields[8])
            seq = fields[9]
            qual = fields[10]
            # Additional optional fields
            optional_fields = fields[11:]
            # Store the information
            data.append({'QueryName': query_name, 'Flag': flag, 'ReferenceName': reference_name, 'Position': position,
                         'MappingQuality': mapq, 'CIGAR': cigar, 'RNext': rnext, 'PNext': pnext, 'TLEN': tlen,
                         'Sequence': seq, 'Quality': qual, 'OptionalFields': optional_fields})
    return data

# Define a dictionary mapping tag abbreviations to descriptive names
TAG_DESCRIPTIONS = {
    'AS': 'Alignment score',
    'XN': 'Number of ambiguous bases',
    'XM': 'Number of mismatches',
    'XO': 'Number of gap opens',
    'XG': 'Number of gap extensions',
    'NM': 'Edit distance to the reference',
    'MD': 'String for mismatching positions',
    'YT': 'Type of read',
    'NH': 'Number of reported alignments'
    # Add other tag descriptions as needed
}

def parse_optional_fields(optional_fields):
    """
    Parse a list of optional fields from a SAM file into a dictionary.

    Args:
        optional_fields (list): A list of strings representing the optional fields in a SAM file.

    Returns:
        dict: A dictionary where the keys are the tag names and the values are the tag values.
    """
    parsed_fields = {}
    
    for field in optional_fields:
        tag, type_and_value = field.split(':', 1)
        tag_name = TAG_DESCRIPTIONS[tag] if tag in TAG_DESCRIPTIONS else tag
        type_, value = type_and_value.split(':', 1)
        
        # Convert the value to the appropriate type
        if type_ == 'i':
            value = int(value)
        elif type_ == 'f':
            value = float(value)
        # For 'Z' (string) and 'H' (hex string), keep as string
        # so no conversion needed
        
        parsed_fields[tag_name] = value
    
    return parsed_fields

def process_alignment(sam_path, keep_duplicates=False):
    alignment = pd.DataFrame(parse_sam_file(sam_path)).set_index('QueryName')

    # Annotate strand
    alignment['strand'] = None
    alignment['strand'][(alignment['Flag'] == 16) | (alignment['Flag'] == 272)] = '-'
    alignment['strand'][(alignment['Flag'] == 0) | (alignment['Flag'] == 256)] = '+'
    # Parse optional fields and add as columns
    optFields = pd.DataFrame(alignment['OptionalFields'].apply(parse_optional_fields).values.tolist(), index=alignment.index)
    alignment = pd.concat([alignment, optFields], axis=1)
    
    # Remove and rename columns
    alignment.drop(['RNext', 'PNext', 'TLEN', 'OptionalFields'], inplace=True, axis=1)
    alignment.rename(columns={'ReferenceName': 'chromosome', 'Position': 'Start_position'}, inplace=True)
    alignment.index.name = 'guide_id'
    # Calculate end position
    alignment['End_position'] = alignment['Start_position'] + alignment['Sequence'].str.len()
    # Remove unmapped reads
    alignment = alignment[alignment['Flag'] != 4]

    alignment = alignment.sort_values('MappingQuality', ascending=False)
    salignment = alignment if keep_duplicates else alignment[~alignment.index.duplicated(keep='first')]
    guide_df = pd.read_csv('../resources/all_gRNA_data.csv', index_col='unique_guide_name')
    tab = pd.merge(guide_df, salignment, left_index=True, right_index=True, how='left')
    
    return tab