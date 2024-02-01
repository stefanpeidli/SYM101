#!/usr/bin/env python
"""\
This script trims reads to the region behind the specified primer sequence.

Usage:
    variable_trimming.py <input> <output> <primer> <guidelength> <primerrange> <threads>
Example:
    variable_trimming.py input.fastq output.fastq GGGGACAGGATGCGGAGAAACAGC 20 0,10 4
"""
__author__ = "Stefan Peidli"
__date__ = "29.01.2023"
__email__ = "stefan.peidli@embl.de"

import argparse
from pathlib import Path
from joblib import Parallel, delayed

import pandas as pd
import numpy as np
from tqdm.auto import tqdm
from scipy.spatial.distance import hamming

### parse arguments ###
parser = argparse.ArgumentParser(description='Optional app description')
parser.add_argument('-i','--input', type=str, help='Input file')
parser.add_argument('-o','--output', type=str, help='Output file')
parser.add_argument('-q','--qc', type=str, help='QC output log file')
parser.add_argument('-p','--primer', type=str, help='Primer sequence')
parser.add_argument('-g','--guidelength', type=int, help='Guide length')
parser.add_argument('-r','--primerrange', type=str, help='Range where start of primer is found. Give like this: 0,10')
parser.add_argument('-t','--threads', type=int, help='Number of threads', default=-1)
args = parser.parse_args()

input_file = args.input
output_file = args.output
qc_file = args.qc
seq_primer = args.primer
pmin, pmax = args.primerrange.split(',')
pmin = int(pmin)
pmax = int(pmax)
guide_length = args.guidelength
threads = args.threads

### functions ###
def qual_to_float(qual):
    return list(map(ord, qual))

def make_substrings(s, K):
    return [s[i: j] for i in range(len(s)) for j in range(i + 1, len(s) + 1) if len(s[i:j]) == K]

def trim_sequence(read, quality):
    """Trim sequence to region behind primer and return trimmed sequence and quality.
    
    Parameters
    ----------
    read : str
        Read sequence
    quality : str
        Read quality
    
    Returns
    -------
    trimmed_seq : str
        Trimmed sequence
    trimmed_qual : str
        Trimmed quality
    primer_idx : int
        Index of primer start
    primer_hit : str
        Primer sequence
    primer_distance : float
        Hamming distance between primer and primer_hit
    primer_certainty : float
        Difference between first and second best primer hit primer_distance
    """
    K = len(seq_primer)
    
    # search primer in expected region
    subs = make_substrings(read[pmin:pmax+K], K)
    weights = make_substrings(quality[pmin:pmax+K], K)
    primer_alignment = pd.DataFrame(index=np.arange(len(subs)), columns=['idx', 'sequence', 'hamming_distance'])
    for j, sub, w in zip(range(len(subs)), subs, weights):
        ham = hamming(list(seq_primer), list(sub), qual_to_float(w))
        primer_alignment.iloc[j] = [j+pmin, sub, ham]
    primer_alignment = primer_alignment.sort_values('hamming_distance')
    primer_idx, primer_hit, primer_distance = primer_alignment.iloc[0]
    primer_certainty = primer_alignment.hamming_distance.iloc[1] - primer_alignment.hamming_distance.iloc[0]
    
    # select guide position after primer
    trimmed_seq = read[primer_idx+K : primer_idx+K+guide_length]
    trimmed_qual = quality[primer_idx+K : primer_idx+K+guide_length]
    
    return trimmed_seq, trimmed_qual, primer_idx, primer_hit, primer_distance, primer_certainty

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

### parallel processing ###
def process(line1, line2, line4):
    index = line1
    sequence = line2.replace('\n', '')
    qual = line4.replace('\n', '')
    trimmed_seq, trimmed_qual, _, _, primer_distance, primer_certainty = trim_sequence(sequence, qual)
    return index, trimmed_seq, trimmed_qual, primer_distance, primer_certainty

# Parallel processing
print('Reading input file...')
with open(input_file, 'r') as INPUT:
    lines = INPUT.readlines()
    n_reads = int(len(lines) / 4)
    print('Trimming reads...')
    results = Parallel(n_jobs=threads)(delayed(process)(line1, line2, line4) for line1, line2, line3, line4 in tqdm(chunker(lines, size=4), total=n_reads))

# write to file (both trimmed fastq and QC log)
print('Writing to file...')
Path(output_file).parent.mkdir(parents=True, exist_ok=True)
Path(qc_file).parent.mkdir(parents=True, exist_ok=True)
with open(output_file, mode='w') as OUTPUT, open(qc_file, 'w') as QC:
    # Track QC metrics
    primer_distance_sum = 0
    primer_certainty_sum = 0
    n_removed = 0
    for index, trimmed_seq, trimmed_qual, primer_distance, primer_certainty in tqdm(results):
        primer_certainty_sum += primer_certainty
        primer_distance_sum += primer_distance
        if (primer_distance > 0.3) or (primer_certainty < 0.25):
            # skip read if primer is not found or too close to second best primer
            n_removed += 1
            continue
        # write trimmed read to file
        OUTPUT.write(index)
        OUTPUT.write(trimmed_seq+'\n')
        OUTPUT.write('+\n')
        OUTPUT.write(trimmed_qual+'\n')
    primer_certainty_sum /= n_reads
    primer_distance_sum /= n_reads
    QC.write(f'Mean_primer_distance\tMean_primer_certainty\tFiltered_reads\tTotal_reads\n')
    QC.write(f'{primer_certainty_sum}\t{primer_distance_sum}\t{n_removed}\t{n_reads}\n')

print('Done!')