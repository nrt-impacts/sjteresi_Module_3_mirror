#!/usr/bin/env python3

__author__ = "Anna Haber"

import logging
import coloredlogs
import os
import sys
import argparse
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

def load_tassel_results(file_name):
    """
    Loads in the .txt file of TASSEL results as a pandas.DataFrame.
    Args:
        path_to_file (str)
    """
    # Load in the file, keeping only relevant columns
    results = pd.read_csv(file_name, sep='\t', header='infer',
                          usecols=['Trait','Marker','Chr','Pos','df','F','p'])
    # Drop rows containing NA values.
    results = results.dropna(axis=0)

    return results

def extract_pval(df):
    """
    Extracts p-values as a 1-D numpy array for FDR correction.
    Args:
        pandas.DataFrame
    """
    # Pull out p-value column and convert to 1-D array.
    pvals = df['p'].values

    return pvals

def correct_pval_fdr(p_array, in_alpha):
    """
    Correct the p-values and output an alpha cutoff value.
    Args:
        1-D array of p-values
        input alpha (float)
    """
    _, pvals_fdr, _, _ = multipletests(p_array, in_alpha, method='fdr_bh')

    return _, pvals_fdr, _, _

def filter_snps(df, corr_pvals, alpha):
    """
    Output a new df of significant SNPs based on original alpha
    and corrected p-values from FDR.
    Args:
        pandas.DataFrame of marker statistics
        corr_pvals = array of corrected p-values
        alpha (float)
    """
    # Add corrected p value as a new column to df.
    corr_pvals = pd.Series(corr_pvals)
    df['p_FDR'] = corr_pvals

    # Now select only the rows containing corrected p-values
    # less than the alpha value.
    filtered = df[df['p_FDR']<alpha]

    return filtered

def save_filtered(df, filename):
    """Save the filtered SNPs to a tab-delimited file."""
    df.to_csv(filename, sep='\t', index=False, header=True)

def validate_args(args, logger):
    """Raise if an input argument is invalid."""

    if not os.path.isfile(args.results_file):
        logger.critical("Argument 'results_file' is not file")
        raise ValueError("%s is not a directory" % (args.results_file))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load input data')
    path_main = os.path.abspath(__file__)

    parser.add_argument('--results_file', '-r', type=str, 
                        help='parent path of TASSEL results file')
    parser.add_argument('--alpha_level', '-a', type=float,
                        help='input alpha level (float) for p-value correction')
    parser.add_argument('--out_filename', '-o', type=str,
                        help='name of output file')
    args = parser.parse_args()
    args.results_file = os.path.abspath(args.results_file)
    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    logger.info("Loading in TASSEL results as df...")
    results_df = load_tassel_results(args.results_file)

    logger.info("Extracting p-values for correction...")
    pvals = extract_pval(results_df)

    logger.info("Correcting p-values via FDR...")
    _, pvals_fdr, _, _ = correct_pval_fdr(pvals, args.alpha_level)

    logger.info("Filtering SNPs for significance...")
    filtered_snps = filter_snps(results_df, pvals_fdr, args.alpha_level)
    save_filtered(filtered_snps, args.out_filename)
