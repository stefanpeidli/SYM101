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
