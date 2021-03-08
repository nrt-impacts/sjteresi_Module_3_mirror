#!/usr/bin/env python3

__author__ = "Anna Haber"

import logging
import coloredlogs
import os
import sys
import argparse
import pandas as pd
from datetime import datetime
from dateutil.parser import parse

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

def import_merge(human_data, multispec_data):
    """
    Import and merge the ground_data_2019.csv and 
    multispec_sept_mean.csv data.
    Args:
        2 instances of path_to_file (str)
    """
    # Import the ground data (relevant columns only)
    ground_data = pd.read_csv(human_data, sep=',', header='infer',
                              usecols=['plot','Anthesis','Silking','Sowing'])
    # Import the multispectral means data (relevant columns only).
    multispec = pd.read_csv(multispec_data, sep=',', header='infer',
                            usecols = ['id','NIR_9_2_19','REDEDGE_9_2_19'])
    # Rename the multispec 'id' column to match the ground_data 'plot' col.
    multispec['plot'] = multispec['id']
    multispec = multispec.drop('id', axis=1)

    # Merge the dataframes on column 'plot'.
    merged_df = ground_data.merge(right=multispec, how='inner', on='plot')
    # Drop rows containing NAs.
    merged_df = merged_df.dropna(axis=0)

    # Covert anthesis, silking, and sowing dates to datetime objects
    # and append as new columns.
    merged_df['Anthesis_Date'] = [datetime.strptime(x, '%m/%d/%Y') for x in 
                                  merged_df['Anthesis']]
    merged_df['Silking_Date'] = [datetime.strptime(x, '%m/%d/%Y') for x in 
                                 merged_df['Silking']]
    merged_df['Sowing_Date'] = [datetime.strptime(x, '%m/%d/%Y') for x in 
                                merged_df['Sowing']]
    # calculate days to anthesis
    # subtract sowing date from anthesis date iteratively for each row
    days_to_anth = []
    for i in range(len(merged_df)):
        curr_anth = merged_df.iloc[i,6]
        curr_sow = merged_df.iloc[i,8]
        dta = (curr_anth - curr_sow).days
        days_to_anth.append(dta)
    # append days_to_anth list as new column
    merged_df['Days_to_Anthesis'] = days_to_anth

    # now do the same for days to silking
    days_to_silk = []
    for i in range(len(merged_df)):
        curr_silk = merged_df.iloc[i,7]
        curr_sow = merged_df.iloc[i,8]
        dts = (curr_silk - curr_sow).days
        days_to_silk.append(dts)
    merged_df['Days_to_Silking'] = days_to_silk

    # Keep only the necessary columns.
    merged_data = merged_df[['plot', 'NIR_9_2_19', 'REDEDGE_9_2_19',
                             'Days_to_Anthesis', 'Days_to_Silking']]

    return merged_data

def obs_key(file_name):
    """
    Import the obs_2019_key.csv data.
    Args:
        path_to_file (str)
    """

    obs_key = pd.read_csv(file_name, sep=',', header='infer')

    return obs_key

def geno_plot_dict(df):
    """
    Extract a dictionary of plot ID vs. genotype from a pandas.DataFrame
    Args:
        pandas.DataFrame (i.e. obs_key generated by previous function)
    """
    # Converts dataframe columns into lists.
    geno = df['Genotype'].tolist()
    plot = df['Plot'].tolist()
    # Creates the dictionary with plotID as key and genotype as value.
    geno_plot_dict = dict(zip(plot, geno))

    return geno_plot_dict

def add_genotype(df, dictionary, plotcol):
    """
    Add a column of genotypes corresponding to the plot ID
    in the same row.
    Args:
        pandas.DataFrame
        Dictionary with plot ID as key and genotype as value
        Column of the df containing plot ID information (str)
    """
    # Creates a new column 'Genotype' using the dictionary.
    df['Genotype'] = df[plotcol].replace(to_replace=dictionary, value=None)
    # Reorders columns so 'Genotype' is in the front (TASSEL format).
    df = df[['Genotype', 'plot', 'Days_to_Anthesis', 'Days_to_Silking',
             'NIR_9_2_19', 'REDEDGE_9_2_19']]

    return df

def validate_args(args, logger):
    """Raise if an input argument is invalid."""

    if not os.path.isfile(args.human_data):
        logger.critical("Argument 'human_data' is not file")
        raise ValueError("%s is not a directory" % (args.human_data))
    if not os.path.isfile(args.multispec_data):
        logger.critical("Argument 'multispec_data' is not file")
        raise ValueError("%s is not a directory" % (args.multispec_data))
    if not os.path.isfile(args.obs_key):
        logger.critical("Argument 'obs_key' is not file")
        raise ValueError("%s is not a directory" % (args.obs_key))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load input data')
    path_main = os.path.abspath(__file__)

    parser.add_argument('--human_data', '-d', type=str, default=os.path.join(
                        path_main, '../../', 'data/ground_data_2019.csv'),
                        help='parent path of human data file')
    parser.add_argument('--multispec_data', '-m', type=str, default=os.path.join(
                        path_main, '../../', 'data/multispec_sept_mean.csv'),
                        help='parent path of multispectral means file')
    parser.add_argument('--obs_key', '-o', type=str, default=os.path.join(
                        path_main, '../../', 'data/obs_2019_key.csv'),
                        help='parent path of observation key')
    args = parser.parse_args()
    args.human_data = os.path.abspath(args.human_data)
    args.multispec_data = os.path.abspath(args.multispec_data)
    args.obs_key = os.path.abspath(args.obs_key)
    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    logger.info('Importing observation key...')
    ObsKey = obs_key(args.obs_key)
    logger.info('Constructing plot ID vs. genotype dictionary...')
    PlotGeno = geno_plot_dict(ObsKey)

    logger.info('Merging and formatting ground and multispectral data...')
    Merged_Data = import_merge(args.human_data, args.multispec_data)
    
    logger.info("Adding 'Genotype' column to data...")
    Merged_Geno = add_genotype(Merged_Data, PlotGeno, 'plot')
    
    logger.info('Saving phenotype data...')
    Merged_Geno.to_csv('../data/merged_pheno_data.csv', index=False, header=True)