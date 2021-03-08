"""
module to read in large numerical 2D SNP matrices.

Author: Serena G. Lotreck, lotrecks@msu.edu
Date: 04/16/2020
"""
import numpy as np
import pandas as pd
import csv

def readData(input_file, na_conversion, shape):
    """
    Reads input file and populates emptyDf. Converts all values to int,
    with NA being coded as the int provided in na_conversion. Keeps first row
    to use as column names.

    parameters:
        input_file: str, path to input file
        na_converstion: int, number greater than 1 to use as code for NA
        shape: tuple of two int, the shape (dimensions) of the data

    returns: a populated df representing the data
    """
    print('Initial shape of the data is {}'.format(shape))
    if headers != 0:
        shape[0] = shape[0] - (headers+1)
    print('Shape after accounting for headers is {}'.format(shape))

    print('Initializing empty dataframe...')
    df = initializeDf(shape)

    print('Reading in the data...')
    with open(input_file) as csvfile:
        reader = csv.DictReader(csvfile)
        rowNum = 0
        columns = []
        for row in reader:
            if rowNum == 0:
                columns.append(list(row.keys()))
                print('Naming columns. Example of column names: {}'.format(columns[:4]))
                df.columns = columns
            for key, value in row.items(): #TODO: see if there's a better way to do this by replacing whole row at once with values in dict
                df.loc[rowNum,key] = value
            print('Just read row number {}'.format(rowNum))
            rowNum += 1

    print('Data has been read! Head: {}'.format(df.head()))

    return df 


def initializeDf(shape):
    """
    Initializes a pandas df with the shape of the data to be read.

    parameters:
        shape, tuple of int: the shape (dimensions) of the data

    returns: emptyDf, an empty array to be used for population,
    datatype of all columns is int
    """
    emptyArray = np.empty(shape,dtype=int)
    emptyDf = pd.DataFrame(emptyArray)

    return emptyDf


def strToInt(na_conversion, convert):
    """
    Converts string numbers to ints, and NA to the int represented by
    na_conversion.

    parameters:
        na_converstion: int, number used to represent NA in the new dataframe
        convert: the str to convert
    """
    if convert == 'NA':
        converted = na_conversion
    else:
        converted = int(convert)

    return converted
