""" Find a variety of network data from SureChemBL graphs.

Catch-all code for getting data from SureChemBL graphs. Instead of downloading
full graphs (~4 GB each), this code will find & store necessary data for
analyses on local machines.

"""

import igraph as ig
import numpy as np
import pickle
from itertools import islice
from itertools import zip_longest
import time
import os
from tqdm import tqdm


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
    print("Loaded gml")
    #print(ig.summary(G))

    # #One time function - to save all gmlz files as pickles to save time
    # pickle.dump(G, file=open("Data/Graphs/G_cpd_" + update + ".p", "wb"))
    # print("Dumped G")

    # G = pickle.load(file=open("Data/Graphs/G_cpd_" + update + ".p", "rb"))
    # print(ig.summary(G))

    return G


def get_degrees_from_graphs():
    """Finds all degrees and id:degree pairs from SureChemBL cpd-cpd graphs

    Assumes all quarterly graphs are in 'Data/Graphs/', saves all degree lists
    and id:degree pairs to pickle files in 'Data/Degrees/' directory
    """
    #List of all quarterly updates (avoids initial data dump)
    updates = [
        "20150401", "20150701", "20151001", "20160101", "20160401", "20160701",
        "20161001", "20170101", "20170401", "20170701", "20171001", "20180101",
        "20180401", "20180701", "20181001", "20190101", "20190401", "20190701",
        "20191001", "20200101", "20200401", "20200701", "20201001", "20210101"
    ]
    test = ["20141231"]
    for update in test:
        for i in range(39):  #All pre-20141231 updates
            start = time.time()
            print("--- Analyzing:", update, i, "---")
            G = read_cpdcpd_graph(update + "_" + str(i))

            degrees = get_degrees(G)
            id_degrees = get_id_degree(G)
            del (G)
            print("Avg Degree:", np.mean(degrees))
            pickle.dump(degrees,
                        file=open(
                            "Data/Degrees/degrees_" + update + "_" + str(i) +
                            ".p", "wb"))

            pickle.dump(id_degrees,
                        file=open(
                            "Data/Degrees/id_degrees_" + update + "_" + str(i) +
                            ".p", "wb"))

            print("Time:", time.time() - start)
            del (degrees)
            del (id_degrees)

        print()


def get_degree_distributions():
    """ Finds and saves degrees across SureChemBL update graphs

    Collects degrees in one place (all_degrees, list of degree distributions)
    and finds the average degree distribution (avg, single list). Saves both
    to the 'Data/Degrees/' directory

    """
    #Load all degrees associated with quarterly updates
    print("\n----- Loading Degree Distributions -----\n")
    all_degrees = []
    for f in tqdm(os.listdir("Data/Degrees/")):
        if f.startswith("degrees_"):
            degrees = pickle.load(file=open("Data/Degrees/" + f, "rb"))
            degrees.sort(reverse=True)
            all_degrees.append(degrees)
    pickle.dump(all_degrees, file=open("Data/Degrees/all_degrees.p", "wb"))

    #Calculate & plot average of all degrees (avg over each column)
    #Note: takes ~20 minutes
    print("\n----- Calculating average degree distribution -----\n")
    avg = [
        np.ma.average(np.ma.masked_values(temp_list, None))
        for temp_list in zip_longest(*all_degrees)
    ]
    pickle.dump(avg, file=open("Data/Degrees/avg_degree_list.p", "wb"))


def calculate_preferential_attachment():
    """ Calculates preferential attachment index (see Rednar 2004) for SureChemBL degrees

    Uses data stored in the 'Data/Degrees/id_degrees_*' files to build a preferential
    attachment index. Saves both the full list of id:degree pairs (summed over time)
    and the preferential attachment indicies of each degree to pickle files in
    'Data/Degrees' directory
    """
    print("\n----- Building id-degree dictionary -----\n")
    full_id_degrees = {}
    size = 63  #there are 63 unique update files

    i = 0
    for f in tqdm(os.listdir("Data/Degrees/")):
        if f.startswith("id_degrees_"):
            #Load id_degree dictionary
            id_degrees = pickle.load(file=open("Data/Degrees/" + f, "rb"))

            # TODO Add / update to full dictionary.

            for key, value in id_degrees.items():
                #If a id is in the dictionary, update the appropriate degree value
                if key in full_id_degrees:
                    full_id_degrees[key][i] = sum(full_id_degrees[key]) + value
                # If a id is not in the dictionary, add it with a list of len(updates) 0s and update the appropriate degree value
                else:
                    full_id_degrees[key] = [0] * size
                    full_id_degrees[key][i] = value
            i += 1

    print(list(islice(full_id_degrees.items(), 10)))

    #Replace all zero values with the last non-zero value
    print("\n----- Replacing zero values -----\n")
    for key, value in tqdm(full_id_degrees.items()):
        arr = np.array(value)
        prev = np.arange(len(arr))
        prev[arr == 0] = 0
        prev = np.maximum.accumulate(prev)
        full_id_degrees[key] = arr[prev]
    pickle.dump(full_id_degrees,
                file=open("Data/Degrees/full_id_degrees.p", "wb"))

    #Calculate preferential attachment
    print("\n----- Calculating preferential attachment -----\n")
    pref_attach_dict = {}
    for key, value in tqdm(full_id_degrees.items()):
        c = 0
        attachments = []
        while c < len(value) - 1:
            attachments.append(value[c + 1] - value[c])
            c += 1
        pref_attach_dict[key] = np.mean(attachments)
    pickle.dump(pref_attach_dict,
                file=open("Data/Degrees/pref_attach_dict.p", "wb"))


def main():
    #Calculate and save degrees & id:degree pairs from SureChemBL updates
    #get_degrees_from_graphs()

    #Store all degree distributions in a single list
    get_degree_distributions()

    #Calculate preferential attachment index for entire SureChemBL dataset
    calculate_preferential_attachment()


if __name__ == "__main__":
    main()
