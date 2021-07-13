""" Builds a bipartite network of the SureChEMBL database

Reads in cpd & patent data (ids & dates) from SureChEMBL data 
in Data\SureChemBLMAP directory. Builds a bipartite network of cpds & patents,
dates associated, with cpds connected to patents if a cpd appears in a patent.
The earliest date of entry (both cpd & patent) is used.

"""

import igraph as ig
import pickle
import pandas as pd
import numpy as np
import time
from tqdm import tqdm
import os


def read_data(fp):
    """ Read in SureChemBL data into a dataframe.

    Generates a dataframe with SureChEMBL cpd ID, patent ID, and associated date.

    Args:
        fp: filepath to .txt file containing SureChemBL data
    
    Returns:
        A pandas dataframe with three columns and N rows. All rows are represented as strings.
    """
    return pd.read_csv(fp,
                       delimiter="\t",
                       usecols=[0, 4, 5],
                       names=["cpdID", "patentID", "Date"])


def get_ids_dates(df, c, id_date_dict, unique_ids):
    """ Gets data (ids & dates) from SureChemBL dataframe

    Builds a list of unique ids, and a dictionary of ids associated with the earliest date of entry. 
    Applied to both compounds or patents.

    Args:
        df: individual dataframe of SureChemBL data (from read_data())
        c: "cpdID" or "patentID", depending on which data is used
        id_date_dict: existing dictionary linking ids & dates
        unique_ids: existing list of unique ids

    Returns:
        list of all unique ids, dictionary of all unique ids with earliest date of entry
    """
    #Find unique compounds
    unique_ids = list(set(df[c].tolist() + unique_ids))

    for index, row in tqdm(df.iterrows(), total=df.shape[0]):
        #Add ID if not present in database
        if row[c] not in id_date_dict:
            id_date_dict[row[c]] = row["Date"]
        else:
            #If the id is there, check if the patent date is earlier than before
            if time.strptime(row["Date"], "%Y-%m-%d") < time.strptime(
                    id_date_dict[row[c]], "%Y-%m-%d"):
                #If so, replace date with earlier time
                id_date_dict[row[c]] = row["Date"]

    return unique_ids, id_date_dict


def build_network(cpds, patents, cpd_date_dict, patent_date_dict):
    """ Builds the igraph network of cpds & patents.

    Takes compounds, patents, and associated dates to build a bipartite igraph network. Compounds are linked
    to patents if a compound appears in a specific patent. Type (cpd vs patent) is specified through 
    the "type" variable.

    Args:
        cpds: list of all unique cpd ids
        patents: list of all unique patent ids
        cpd_date_dict: dictionary associating all cpd ids with the earliest date of entry
        patent_date_dict: dictionary associating all patent ids with the earliest date of entry

    Returns:
        igraph network G
"""
    G = ig.Graph()

    ### Add nodes ###

    ### Link nodes ###

    return G


def main():
    ### Read in data ###
    #Build list of all unique compounds & patents, as well as dictionaries with dates
    cpds = []
    cpd_dates = {}
    patents = []
    patent_dates = {}

    for f in os.listdir("Data/SureChemblMAP/"):
        print("---- Analzying", f, "----")
        df = read_data("Data/SureChemblMAP/" + f)

        cpds, cpd_dates = get_ids_dates(df, "cpdID", cpd_dates, cpds)
        print("Current unique cpds:", len(cpds))

        patents, patent_dates = get_ids_dates(df, "patentID", patent_dates,
                                              patents)
        print("Current unique patents:", len(patents))

    #Save cpds & patents (including dictionaries)
    pickle.dump(cpds, file=open("Data/CpdPatentIdsDates/nique_cpds.p", "wb"))
    pickle.dump(cpd_dates,
                file=open("Data/CpdPatentIdsDates/cpd_date_dict.p", "wb"))
    pickle.dump(patents,
                file=open("Data/CpdPatentIdsDates/unique_patents.p", "wb"))
    pickle.dump(patent_dates,
                file=open("Data/CpdPatentIdsDates/patent_date_dict.p", "wb"))

    ### Create graph ###

    # TODO(JFM): Do this
    G = build_network()


if __name__ == "__main__":
    main()
