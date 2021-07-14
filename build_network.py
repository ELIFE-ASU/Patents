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


def build_network(cpds, patents, cpd_date_dict, patent_date_dict, edges):
    """ Builds the igraph network of cpds & patents.

    Takes compounds, patents, and associated dates to build a bipartite igraph network. Compounds are linked
    to patents if a compound appears in a specific patent. Type (cpd vs patent) is specified through 
    the "type" variable.

    Args:
        cpds: list of all unique cpd ids
        patents: list of all unique patent ids
        cpd_date_dict: dictionary associating all cpd ids with the earliest date of entry
        patent_date_dict: dictionary associating all patent ids with the earliest date of entry
        edges: list of tuples in (cpd_id, patent_id) form

    Returns:
        igraph network G
"""
    G = ig.Graph()

    ### Add nodes ###
    G.add_vertices(len(cpds) + len(patents))
    G.vs["name"] = cpds + patents  #name field is the respective cpd/patent ID

    #Add dates to each node
    all_dates = []
    for cpd in cpds:
        all_dates.append(cpd_date_dict[cpd])
    for patent in patents:
        all_dates.append(patent_date_dict[patent])
    G.vs["date"] = all_dates

    #Types - cpd:0, patent:1
    G.vs["type"] = [0] * len(cpds) + [1] * len(patents)

    ### Add edges ###
    indexed_edges = []
    #Find the index of each cpd & patent id
    for edge in edges:
        #Use G.vs.find to take advantage of igraph dictionary
        indexed_edges.append(
            (G.vs.find(edge[0]).index, G.vs.find(edge[1]).index))
    G.add_edges(indexed_edges)

    print(ig.summary(G))
    return G


def get_cpd_patent_info(data_fp):
    """ Saves compound and patent information from SureChemBL mapping

    Reads in all SureChemBL mapping data and creates lists of all unique compounds & patents, as well as 
    dictionaries which link compounds and patents to the earliest data of entry.

    Args:
        data_fp: filepath to SureChemBL mapping data(pre-downloaded)

    Returns:
        None, but saves all data to pickle files to "Data/CpdPatentIdsDates" directory
    
    """

    for f in os.listdir(data_fp):  #full dataset
        # f = "SureChEMBL_map_20210101.txt"  #test case
        ## Note: variable declarations should be outside the for loop for full dataset analysis
        cpds = []
        cpd_dates = {}
        patents = []
        patent_dates = {}

        print("---- Analzying", f, "----")
        df = read_data(data_fp + f)

        cpds, cpd_dates = get_ids_dates(df, "cpdID", cpd_dates, cpds)
        print("Current unique cpds:", len(cpds))

        patents, patent_dates = get_ids_dates(df, "patentID", patent_dates,
                                              patents)
        print("Current unique patents:", len(patents))

        #Save cpds & patents (including dictionaries)
        ## Note: pickle dumps *should* be outside the for loop for full dataset analysis
        pickle.dump(cpds,
                    file=open(
                        "Data/CpdPatentIdsDates/unique_cpds" + f[15:-4] + ".p",
                        "wb"))
        pickle.dump(
            cpd_dates,
            file=open("Data/CpdPatentIdsDates/cpd_date_dict" + f[15:-4] + ".p",
                      "wb"))
        pickle.dump(
            patents,
            file=open("Data/CpdPatentIdsDates/unique_patents" + f[15:-4] + ".p",
                      "wb"))
        pickle.dump(patent_dates,
                    file=open(
                        "Data/CpdPatentIdsDates/patent_date_dict" + f[15:-4] +
                        ".p", "wb"))


def get_cpd_patent_relations(data_fp):
    """ Builds a relation between compounds and patents

    Relates compounds to patents where they appear. The data structure used is a list of tuples
    in (cpd, patent) form, with the ids of cpds & patents used as unique identifiers.

    Args:
        data_fp: filepath to SureChemBL mapping

    Returns:
        None, but saves all data to pickle files to "Data/CpdPatentIdsDates" directory
    
    """

    for f in os.listdir(data_fp):
        #f = "SureChEMBL_map_20210101.txt"  #test case
        ## Note: list definition should be outside for loop for full dataset
        cpd_patent_edges = []

        print("---- Analzying", f, "----")
        df = read_data(data_fp + f)

        #Extent an existing list of tuples, each tuple is a relation between (cpd, patent)
        cpd_patent_edges.extend(
            list(set([t for t in list(zip(df["cpdID"], df["patentID"]))])))

        ## Note: pickle dump should be outside for loop for full dataset
        pickle.dump(cpd_patent_edges,
                    file=open(
                        "Data/CpdPatentIdsDates/cpd_patent_edges" + f[15:-4] +
                        ".p", "wb"))


def main():
    ### Read in data ###

    # #Build list of all unique compounds & patents, as well as dictionaries with dates
    get_cpd_patent_info("Data/SureChemblMAP/")

    # #Build edge list (cpd, patent)
    print("\n\n--- Edges --- \n")
    get_cpd_patent_relations("Data/SureChemblMAP/")

    # # ### Create graph ###
    # # #Note - takes ~90GB and ~20 minutes to build the full network
    # G = build_network(
    #     pickle.load(
    #         file=open("Data/CpdPatentIdsDates/unique_cpds_2021.p", "rb")),
    #     pickle.load(
    #         file=open("Data/CpdPatentIdsDates/unique_patents_2021.p", "rb")),
    #     pickle.load(
    #         file=open("Data/CpdPatentIdsDates/cpd_date_dict_2021.p", "rb")),
    #     pickle.load(
    #         file=open("Data/CpdPatentIdsDates/patent_date_dict_2021.p", "rb")),
    #     pickle.load(
    #         file=open("Data/CpdPatentIdsDates/cpd_patent_edges_2021.p", "rb")))

    # pickle.dump(G, file=open("Data/Graphs/G_full_2021.p", "wb"))

    # ### Cpd & Patent subgraphs ###
    # #Test with 2021 network first
    # G = pickle.load(file=open("Data/Graphs/G_full_2021.p", "rb"))
    # G_cpd, G_patent = G.bipartite_projection()

    # print("--- Compound Projection --- ")
    # print(ig.summary(G_cpd), end="\n\n")
    # print("--- Patent projection ---")
    # print(ig.summary(G_patent))


if __name__ == "__main__":
    main()
