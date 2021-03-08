#!/usr/bin/env python3

__author__ = "Anna Haber"

import logging
import coloredlogs
import os
import sys
import argparse
import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

def load_results(file_name):
    """
    Loads in the .txt file of TASSEL results as a pandas.DataFrame.
    Args:
        path_to_file (str)
    """
    # Load in the file, keeping only relevant columns
    results = pd.read_csv(file_name, sep=',', header='infer',
        usecols=['SNP','Chromosome','Position ','P.value','FDR_Adjusted_P-values'])
    # Drop rows containing NA values.
    results = results.dropna(axis=0)

    return results

def filter_snps(df, alpha):
    """
    Output a new df of significant SNPs based on original alpha
    and corrected p-values from FDR.
    Args:
        pandas.DataFrame of marker statistics
        alpha (float)
    """
    # Select only the rows containing corrected p-values
    # less than the alpha value.
    filtered = df[df['FDR_Adjusted_P-values']<alpha]

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
                        help='parent path of GAPIT results file')
    parser.add_argument('--alpha_level', '-a', type=float,
                        help='input alpha level (float) for p-value correction')
    parser.add_argument('--out_filename', '-o', type=str,
                        help='name of output file')
    args = parser.parse_args()
    args.results_file = os.path.abspath(args.results_file)
    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    logger.info("Loading in GAPIT results as df...")
    results_df = load_results(args.results_file)

    logger.info("Filtering SNPs for significance...")
    filtered_snps = filter_snps(results_df, args.alpha_level)
    logger.info("Saving file of filtered SNPs...")
    save_filtered(filtered_snps, args.out_filename)
