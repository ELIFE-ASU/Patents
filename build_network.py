""" Builds a bipartite network of the SureChEMBL database

Reads in cpd & patent data (ids & dates) from SureChEMBL data
in Data\SureChemBLMAP directory. Builds a bipartite network of cpds & patents,
dates associated, with cpds connected to patents if a cpd appears in a patent.
The earliest date of entry (both cpd & patent) is used.

"""

from numpy.lib.shape_base import split
import igraph as ig
import pickle
import pandas as pd
import numpy as np
import time
import datetime
from tqdm import tqdm
import os
from collections import defaultdict
from itertools import combinations
from itertools import islice
from itertools import repeat
from multiprocessing import Pool
from functools import partial
import calendar


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


def build_bipartite_network(cpds, patents, cpd_date_dict, patent_date_dict,
                            edges):
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

    # #List of all quarterly updates (avoids initial data dump)
    #TODO: code skipped 20150401 for some reason
    updates = [
        "20150401", "20150701", "20151001", "20160101", "20160401", "20160701",
        "20161001", "20170101", "20170401", "20170701", "20171001", "20180101",
        "20180401", "20180701", "20181001", "20190101", "20190401", "20190701",
        "20191001", "20200101", "20200401", "20200701", "20201001", "20210101"
    ]
    test = ["20150401"]
    for update in updates:  # in os.listdir(data_fp):  #full dataset
        f = "SureChEMBL_map_" + update + ".txt"

        ## Note: variable declarations should be outside the for loop for full dataset analysis
        cpds = []
        cpd_dates = {}
        patents = []
        patent_dates = {}

        print("---- Analzying", f, "----")
        df = read_data(data_fp + f)

        months = get_all_months(df["Date"].iloc[0], df["Date"].iloc[-1])

        for month in months:
            criterion = df["Date"].map(lambda x: x.startswith(month[:-2]))
            split = df[criterion]
            print("\n----- Building", month[:-3], "-----")

            print("-- Unique Cpds --")
            cpds, cpd_dates = get_ids_dates(split, "cpdID", cpd_dates, cpds)

            print("-- Unique Patents --")
            patents, patent_dates = get_ids_dates(split, "patentID",
                                                  patent_dates, patents)

            #Build patent-cpd relationships for each month
            get_cpd_patent_relations(split, "_" + month[:-3])

            #Save cpds & patents (including dictionaries)
            ## Note: pickle dumps *should* be outside the for loop for full dataset analysis
            pickle.dump(cpds,
                        file=open(
                            "Data/CpdPatentIdsDates/unique_cpds_" + month[:-3] +
                            ".p", "wb"))
            pickle.dump(cpd_dates,
                        file=open(
                            "Data/CpdPatentIdsDates/cpd_date_dict_" +
                            month[:-3] + ".p", "wb"))
            pickle.dump(patents,
                        file=open(
                            "Data/CpdPatentIdsDates/unique_patents_" +
                            month[:-3] + ".p", "wb"))
            pickle.dump(patent_dates,
                        file=open(
                            "Data/CpdPatentIdsDates/patent_date_dict_" +
                            month[:-3] + ".p", "wb"))


def get_cpd_patent_relations(df, label):
    """ Builds a relation between compounds and patents

    Relates compounds to patents where they appear. The data structure used is a list of tuples
    in (cpd, patent) form, with the ids of cpds & patents used as unique identifiers.
    Also relates patents to all compounds which appear in them. This is a dictionary in
    {patent: [cpds]} form, also with ids as unique identifiers.

    Args:
        df: dataframe containing SureChemBL map data

    Returns:
        None, but saves all data to pickle files to "Data/CpdPatentIdsDates" directory

    """

    cpd_patent_edges = []
    patent_cpd_edges = defaultdict(list)

    #Extend an existing list of tuples, each tuple is a relation between (cpd, patent)
    cpd_patent_edges.extend(
        list(set([t for t in list(zip(df["cpdID"], df["patentID"]))])))

    #Builds a dictionary of {patent: [cpd]} relations
    print("-- Cpd-patent Edges --")
    for index, row in tqdm(df.iterrows(), total=df.shape[0]):
        patent_cpd_edges[row["patentID"]].append(row["cpdID"])

    ## Note: pickle dump should be outside for loop for full dataset
    pickle.dump(cpd_patent_edges,
                file=open(
                    "Data/CpdPatentIdsDates/cpd_patent_edges" + label + ".p",
                    "wb"))

    pickle.dump(patent_cpd_edges,
                file=open(
                    "Data/CpdPatentIdsDates/patent_cpd_edges" + label + ".p",
                    "wb"))


def build_cpd_network(cpds, cpd_date_dict, patent_cpd_links):
    """ Builds a network of compounds, connected by occurrence within the same patent.

    Builds a igraph network with SureChemBL compounds as nodes, and edges between compounds
    are created when two compounds appear in the same patent together.

    Args:
        cpds: list of all unique compounds in SureChemBL
        cpd_date_dict: all compounds associated with date of first entry
        patent_cpd_links: finds all compounds associated with each patent

    Returns:
        an igraph network of SureChemBL compounds

    """
    G = ig.Graph()

    ### Add nodes ###
    G.add_vertices(len(cpds))
    G.vs["name"] = cpds  #name field is the respective cpd/patent ID

    #Add dates to each node
    all_dates = []
    for cpd in cpds:
        all_dates.append(cpd_date_dict[cpd])
    G.vs["date"] = all_dates

    ### Add vertices ###
    # pool = Pool(processes=10)
    # patent_cpd_slice = dict(islice(patent_cpd_links.items(), 1000))
    # print(patent_cpd_slice)

    ## Note: serial attempt
    c = 0  #loop counter
    es = []
    for s in tqdm(patent_cpd_links.values()):
        if c < 10000:  #Every 10k loops, add edges to the graph (otherwise store them)
            es.extend(find_cpd_cpd_edges(G, s))
            c += 1
        else:
            G.add_edges(es)
            es = []  #Reset the counter and edge list
            c = 0

    #Add any extra edges4
    G.add_edges(es)

    ## Note: parallel attempts
    # start = time.time()
    # graphs = pool.map(partial(find_cpd_cpd_edges, G), patent_cpd_slice.values())

    # for graph in tqdm(graphs):
    #     G.add_edges(graph.get_edgelist())

    # print("Took:", time.time() - start)

    # print("Adding edges...")
    # for es in tqdm(edges):
    #     G.add_edges(list(es))

    print(ig.summary(G))
    return G


def find_cpd_cpd_edges(graph, patent_cpds):
    """" Adds edges to an igraph network G

    Takes a dictionary of patent-cpd links in {patent: [cpds]} form, adds edges between
    all compounds associated with a single patent. Meant to be called in parallel from
    build_cpd_network().

    Args:
        G: igraph network, with nodes as compounds
        patent_cpd_links: dictionary of all compounds associated with patents

    Returns:
        igraph network with edges between compounds

    """
    #try:  #TODO: assert that patent_cpd is a list
    #Take advantage of O(1) lookup of graphs
    linked_cpds_indicies = [
        graph.vs.find(c).index for c in list(set(patent_cpds))
    ]
    return list(combinations(linked_cpds_indicies, 2))

    # ## Note: parallel attempt to add edges directly to G
    # # print(linked_cpds_indicies)
    # #Add edges between all possible combinations of nodes found within a patent
    # graph.add_edges(list(combinations(linked_cpds_indicies, 2)))
    # print(ig.summary(graph))

    # return graph

    # else:
    #     print("Not a dictionary", type(patent_cpd_links))
    #     print(patent_cpd_links)


def subtract_month(date):
    """ Subtracts one month from a given date

    Args:
        date (string, in form YYYY-MM-DD): date to be subtracted from

    Returns:
        (string, in form YYYY-MM-DD): date after subtracting one month. The
        number of days is the maximum number of days in that particular month.
    """
    date = datetime.datetime.strptime(date, "%Y-%m-%d")
    month = date.month - 2
    month = month % 12 + 1
    year = date.year - month // 12
    day = calendar.monthrange(year, month)[-1]
    return datetime.date(year, month, day).strftime("%Y-%m-%d")


def get_all_months(start, end):
    """Build a list of all months within a given range.

    Takes a range of dates from a dataframe, then builds a list of all months
    from that range in the form YYYY-MM-DD, where DD is the last day of each month

    Args:
        start (string, in form YYYY-MM-DD): start date
        end ((string, in form YYYY-MM-DD)): end date

    Returns:
        range (list): list of dates in YYYY-MM-DD format
    """
    months = []
    while end > start:
        months.append(end)
        end = subtract_month(end)
    return months


def main():
    ### Read in data ###

    # #Build list of all unique compounds & patents, as well as dictionaries with dates
    get_cpd_patent_info("Data/SureChemblMAP/")

    ### Create cpd-patent graph ###
    #Note - takes ~90GB and ~20 minutes to build the full network

    # for update in updates:
    #     print("--- Building:", update, "---")
    #     G = build_bipartite_network(
    #         pickle.load(file=open(
    #             "Data/CpdPatentIdsDates/unique_cpds" + update +
    #             ".p", "rb")),
    #         pickle.load(file=open(
    #             "Data/CpdPatentIdsDates/unique_patents" + update +
    #             ".p", "rb")),
    #         pickle.load(file=open(
    #             "Data/CpdPatentIdsDates/cpd_date_dict" + update +
    #             ".p", "rb")),
    #         pickle.load(file=open(
    #             "Data/CpdPatentIdsDates/patent_date_dict" + update +
    #             ".p", "rb")),
    #         pickle.load(file=open(
    #             "Data/CpdPatentIdsDates/cpd_patent_edges" + update +
    #             ".p", "rb")))

    #     pickle.dump(G, file=open("Data/Graphs/G_" + update + ".p", "wb"))

    ### Cpd & Patent subgraphs ###
    # G_cpd, G_patent = G.bipartite_projection(multiplicity=False)
    # pickle.dump(G_cpd,
    #             file=open("Data/Graphs/G_cpd_" + update + ".p", "wb"))
    # pickle.dump(G_patent,
    #             file=open("Data/Graphs/G_patent_" + update + ".p", "wb"))
    # #Test size (for possible later projections)
    # sizes = G.bipartite_projection_size()
    # print(sizes)
    # print(update + "," + str(sizes[0]) + "," + str(sizes[1]) + "," +
    #         str(sizes[2]) + "," + str(sizes[3]),
    #         file=f)

    # ### Build cpd-cpd graph ###
    # test = ["20141231"]
    # print("\n\n --- Building Graphs --- \n")
    # for update in test:
    #     print("--- Building:", update, "---")
    #     for c in range(17, 38):  # len(pre-2014 data) / 5.5 million = 38
    #         G = build_cpd_network(
    #             pickle.load(file=open(
    #                 "Data/CpdPatentIdsDates/unique_cpds" + update +
    #                 ".p", "rb")),
    #             pickle.load(file=open(
    #                 "Data/CpdPatentIdsDates/cpd_date_dict" + update +
    #                 ".p", "rb")),
    #             pickle.load(file=open(
    #                 "Data/CpdPatentIdsDates/patent_cpd_edges" + update + "_" +
    #                 str(c) + ".p", "rb")))

    #         # pickle.dump(G, file=open("Data/Graphs/G_cpd_" + update + ".p", "wb")) ## Too much memory
    #         G.save("Data/Graphs/G_cpd_" + update + "_" + str(c) + ".gmlz",
    #                format="graphmlz")  #save in zipped-gml format to save memory
    #         del (G)  #remove G from memory to free up space

    #     #TODO: rebuild 20210101 cpd-patent graph (overwrote it accidentally)


if __name__ == "__main__":
    main()
