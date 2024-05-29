"""
Annotate exon number for each guide RNA in the SAM file

- Parse the alignment results (SAM/BAM files) to extract the alignment positions (chromosome, start, end).
- Load the transcript annotation file (GTF/GFF) into memory.
- For each alignment position, check which transcripts it overlaps with based on genomic coordinates.
- Determine the exon(s) within the transcript that the alignment overlaps with.
"""


import pandas as pd
import numpy as np
from tqdm.auto import tqdm
from gtfparse import read_gtf

# Process CLI
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_sam", help="Input SAM file", required=True)
parser.add_argument("-g", "--gtf", help="GTF file", required=True)
parser.add_argument("-x", "--extra_output", help="Exon number series", required=True)
parser.add_argument("-o", "--annotated_csv", help="Output csv SAM-like file", required=True)
args = parser.parse_args()

# Functions
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

def process_alignment(sam_path):
    alignment = pd.DataFrame(parse_sam_file(sam_path)).set_index('QueryName')

    # Annotate strand
    alignment['strand'] = None
    alignment['strand'][(alignment['Flag'] == 16) | (alignment['Flag'] == 272)] = '-'
    alignment['strand'][(alignment['Flag'] == 0) | (alignment['Flag'] == 256)] = '+'
    # Remove and rename columns
    alignment.drop(['RNext', 'PNext', 'TLEN', 'OptionalFields'], inplace=True, axis=1)
    alignment.rename(columns={'ReferenceName': 'chromosome', 'Position': 'Start_position'}, inplace=True)
    alignment.index.name = 'guide_id'
    # Calculate end position
    alignment['End_position'] = alignment['Start_position'] + alignment['Sequence'].str.len()
    # Remove unmapped reads
    alignment = alignment[alignment['Flag'] != 4]

    alignment = alignment.sort_values('MappingQuality', ascending=False)
    salignment = alignment[~alignment.index.duplicated(keep='first')]
    guide_df = pd.read_csv('../resources/all_gRNA_data.csv', index_col='unique_guide_name')
    tab = pd.merge(guide_df, salignment, left_index=True, right_index=True, how='left')
    
    return tab

# Read
alignment = pd.DataFrame(parse_sam_file(args.input_sam)).set_index('QueryName')
alignment['End'] = alignment['Position'] + alignment['Sequence'].str.len()

# read gtf file (is this the right one??)
gtf = read_gtf(args.gtf)
gtf = pd.DataFrame(gtf, columns=gtf.columns)

# Annotate exon number
# Problem: I find multiple exons at the same location
best_hits = {}
for i in tqdm(range(len(alignment))):
    if alignment.iloc[i, 3] == 0: # 'MappingQuality'
        best_hits[i] = None
        continue
    ref, start, end = alignment.iloc[i, [1, 2, -1]]
    sgtf = gtf[(gtf.seqname==ref) & (gtf.feature == 'exon')]
    is_overlapping = (((sgtf.start < start) & (sgtf.end > start)) | ((sgtf.start < end) & (sgtf.end > end)))
    sgtf = sgtf[is_overlapping]
    exons = sgtf.exon_number.unique()
    if len(exons) == 0:
        # no matching exon
        exon_number = None
    elif len(exons) == 1:
        # exactly one exon matching
        exon_number = exons[0]
    elif len(exons) > 1:
        # Multiple different exon numbers for position were found
        if len(sgtf.exon_version) > 1:
            # Take most recent exon version
            sgtf_ordered = sgtf.sort_values('exon_version', ascending=False)
        # Now we just take the first exon number
        exon_number = sgtf_ordered.exon_number.iloc[0]
    best_hits[i] = exon_number
    if i % 1000 == 0:
        print(i, exon_number)
print('Unique exon numbers: ', len(pd.unique(list(best_hits.values()))))
s = pd.Series(best_hits, name='exon_number')
s.index.name = 'guide_id'
s.to_csv(args.extra_output)

print("writing to file")
alignment['exon_number'] = s.values
alignment.to_csv(args.annotated_csv)