""" Find a variety of network data from SureChemBL graphs.

Catch-all code for getting data from SureChemBL graphs. Instead of downloading
full graphs (~4 GB each), this code will find & store necessary data for
analyses on local machines.

"""

import igraph as ig
import numpy as np
import pickle
from itertools import islice
import time


def get_degrees(G):
    """ Finds the degree distribution of a igraph network

    Args:
        G (igraph object): igraph network, contains a variety of data including degrees between nodes

    Returns:
        G.degree (list): list of degrees (in order of igraph vertex index)
    """
    return G.degree()


def get_id_degree(G):
    """ Creates a dictionary of compound ids & associated degrees

    Args:
        G (igraph object): igraph network, contains a variety of data including degrees between nodes

    Returns:
        id_degree_dict (dictionary): Associates SurechemBL compound ids with degree value
    """
    id_degree_dict = dict(zip(G.vs["name"], G.degree()))

    return id_degree_dict


def read_cpdcpd_graph(update):
    """ Reads a cpd-cpd graph in .gmlz format

    Takes a stored igraph network in .gmlz format and reads it into igraph form.
    Assumes that the data are stored in the Data/Graphs/ directory and in the format
    G_cpd_<update>.gmlz.

    Args:
        update (string): date (YYYYMMDD) of SureChemBL data dump

    Returns:
        G (igraph object): cpd-cpd network pertaining to a specific update
    """
    G = ig.Graph.Read_GraphMLz("Data/Graphs/G_cpd_" + update + ".gmlz")
    # print("Loaded gml")
    # print(ig.summary(G))

    # #One time function - to save all gmlz files as pickles to save time
    # pickle.dump(G, file=open("Data/Graphs/G_cpd_" + update + ".p", "wb"))
    # print("Dumped G")

    # G = pickle.load(file=open("Data/Graphs/G_cpd_" + update + ".p", "rb"))
    # print(ig.summary(G))

    return G


def main():
    #List of all quarterly updates (avoids initial data dump)
    updates = [
        "20150401", "20150701", "20151001", "20160101", "20160401", "20160701",
        "20161001", "20170101", "20170401", "20170701", "20171001", "20180101",
        "20180401", "20180701", "20181001", "20190101", "20190401", "20190701",
        "20191001", "20200101", "20200401", "20200701", "20201001", "20210101"
    ]
    test = ["20150401"]
    for update in updates:
        start = time.time()
        print("--- Analyzing:", update, "---")
        G = read_cpdcpd_graph(update)

        degrees = get_degrees(G)
        print("Avg Degree:", np.mean(degrees))
        pickle.dump(degrees,
                    file=open("Data/Degrees/degrees_" + update + ".p", "wb"))
        id_degrees = get_id_degree(G)
        pickle.dump(id_degrees,
                    file=open("Data/Degrees/id_degrees_" + update + ".p", "wb"))

        print("Time:", time.time() - start)

        print()


if __name__ == "__main__":
    main()
