import pandas as pd
import numpy as np
import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
# from sklearn.datasets import load_iris

import argparse
import os
import timeit

import matplotlib.pyplot as plt

# import logging
# import coloredlogs

# adapted from
# https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60


def readData(input_file):
    """
    Reads in snp matrix and drops rows with nan

    args:
        input_file (str): path to input file

    returns:
        data (Pandas DataFrame): Interpolatated dataframe, snps columns, each
        row a genotype, but that information is contained in first_col

        first col (Pandas Series): 1 column, the genotypes, string type
    """
    range1 = [i for i in range(0,200000)]
    data = pd.read_csv(input_file,
                        sep='\t',
                        header=None,
                        names=range1,
                        low_memory=False)
    data = data.iloc[1:]

    first_col = data.iloc[:,0]  # first column (names)
    first_col.columns = first_col.iloc[0]  # make first row as the header
    first_col = first_col.drop(1, axis=0)  # drop the first row now
    first_col = first_col.drop(2, axis=0)  # drop the second row, contains NA
    first_col = first_col.astype(str)

    data = data.drop(columns = [0])  # drop the first column
    data.columns = data.iloc[0]  # make first row as the header
    data = data.drop(1, axis=0)  # drop the first row now
    data = data.apply(pd.to_numeric)

    #print(f"Before: {data[data.isna().any(axis=1)]}")
    data = data.interpolate()
    #print(f"After: {data[data.isna().any(axis=1)]}")

    data = data.dropna(axis=0, how='any')  # drop rows with NA, drops 1 row
    return data, first_col

    # We have now interpolated values

def standardizeData(snpsNumeric):
    """
    Use StandardScaler to standardize features onto unit scale

    Args:
        snpsNumeric: numeric SNP matrix as pandas df

    Returns:
        (numpy.ndarray): The normalized data
    """
    features = snpsNumeric.columns.values.tolist()
    x = snpsNumeric.loc[:, features].values
    x = StandardScaler().fit_transform(x)
    return x

def getPrincipalComponents(standardized_data, id_col):
    """
    Get principal components for the data

    Args:
        data (numpy.ndarray): the standardized data received from standardizeData

    Returns:
        finalDf (Pandas DataFrame): Dataframe of the principal components
    """
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(standardized_data)
    principalDf = pd.DataFrame(data=principalComponents)
    finalDf = pd.concat([principalDf, id_col], axis=1)
    return finalDf, pca.explained_variance_ratio_


def plotPCs(principalDf, pca_variance):
    """
    TODO add stuff
    """
    x = principalDf.iloc[:,0]
    y = principalDf.iloc[:,1]
    plt.scatter(x,y)
    plt.title('Maize Population Structure PCA')
    plt.xlabel("Principal Component 1: {:.2f}%".format(pca_variance[0]))
    plt.ylabel("Principal Component 2: {:.2f}%".format(pca_variance[1]))
    plt.show()


    #TODO: save figure
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Population density PCA')
    path_main = os.path.abspath(__file__)
    parser.add_argument('input_file', type=str, help='numeric SNP matrix')
    #parser.add_argument('id_col', type=str, help='String of the ID col')


    args = parser.parse_args()
    args.input_file = os.path.abspath(args.input_file)

    data, first_col = readData(args.input_file)
    standardized_data = standardizeData(data)
    final_df, pca_variance = getPrincipalComponents(standardized_data, first_col)

    plotPCs(final_df, pca_variance)