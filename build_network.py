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
from collections import defaultdict
from itertools import combinations
import calendar
import subprocess
import os
from tqdm import tqdm


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

    for index, row in df.iterrows():  #, total=df.shape[0]:
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

    #print(ig.summary(G))
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
        "20141231", "20150401", "20150701", "20151001", "20160101", "20160401",
        "20160701", "20161001", "20170101", "20170401", "20170701", "20171001",
        "20180101", "20180401", "20180701", "20181001", "20190101", "20190401",
        "20190701", "20191001", "20200101", "20200401", "20200701", "20201001",
        "20210101"
    ]
    recent_updates = [
        "20210401", "20210701", "20211001", "20220101", "20220401", "20220401",
        "20220701", "20221001", "20230101"
    ]
    test = ["20141231"]
    test = ["20150401"]  #Testing Agave
    for update in recent_updates[0:2]:  # in os.listdir(data_fp):  #full dataset
        f = "SureChEMBL_map_" + update + ".txt"

        print("---- Analzying", f, "----")
        df = read_data(data_fp + f)
        df = df.sort_values(by=["Date"])
        print(df.head())

        months = get_all_months(df["Date"].iloc[0], df["Date"].iloc[-1])

        for month in months:
            ## Note: variable declarations should be outside the for loop for full dataset analysis
            cpds = []
            cpd_dates = {}
            patents = []
            patent_dates = {}

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
            volume_fp = "/Volumes/Macintosh HD 4/SureChemBL/CpdPatentIdsDates/"
            pickle.dump(cpds,
                        file=open(
                            volume_fp + "unique_cpds_" + month[:-3] + ".p",
                            "wb"))
            pickle.dump(cpd_dates,
                        file=open(
                            volume_fp + "cpd_date_dict_" + month[:-3] + ".p",
                            "wb"))
            pickle.dump(patents,
                        file=open(
                            volume_fp + "unique_patents_" + month[:-3] + ".p",
                            "wb"))
            pickle.dump(patent_dates,
                        file=open(
                            volume_fp + "patent_date_dict_" + month[:-3] + ".p",
                            "wb"))

            # # Move all files to Google Drive
            # subprocess.run([
            #     "rclone", "moveto", "/scratch/jmalloy3/CpdPatentIdsDates",
            #     "SureChemBL_Patents:CpdPatentIdsDates"
            # ])


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
    print("-- Cpd-Patent Edges --")
    for index, row in df.iterrows():  #, total=df.shape[0]:
        patent_cpd_edges[row["patentID"]].append(row["cpdID"])

    ## Note: pickle dump should be outside for loop for full dataset
    volume_fp = "/Volumes/Macintosh HD 4/SureChemBL/CpdPatentIdsDates/"
    pickle.dump(cpd_patent_edges,
                file=open(volume_fp + "cpd_patent_edges" + label + ".p", "wb"))

    pickle.dump(patent_cpd_edges,
                file=open(volume_fp + "patent_cpd_edges" + label + ".p", "wb"))


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
    for s in patent_cpd_links.values():
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

    #print(ig.summary(G))
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


def build_month_list(start, end):
    """ Builds a list of all months in a given range

    Args:
        start (int): year describing the start of the data
        end (int): year describing the end of the data (inclusive)

    Returns:
        list: list of all update months in format "YYYY-MM"
    """
    updates = []
    for year in range(start, end + 1):  # all years through the given end
        for month in range(1, 13):  #include 12 months
            if month < 10:
                updates.append(str(year) + "-0" + str(month))
            else:
                updates.append(str(year) + "-" + str(month))

    return updates


def build_cpd_ID_mapping(fp):
    """ Builds a dictionary mapping SureChemBL IDs to numerical indicies,
    for ease of building an igraph network

    Args:
        fp (string): filepath to location of compound data

    Returns:
        none: does save dictionary to cp_ID_index_dict.p
    """
    allcpds = pd.read_pickle(fp + "SureChemBL_allCpds.p")
    unique_cpds = list(set(allcpds["SureChEMBL_ID"].tolist()))
    print("Cpds with IDs:", len(allcpds["SureChEMBL_ID"].tolist()))
    print("All compounds:", len(allcpds))
    print("Unique cpds:", len(list(set(allcpds["SureChEMBL_ID"].tolist()))))

    cpd_dict = dict(zip(unique_cpds, np.arange(0, len(unique_cpds), 1)))

    pickle.dump(cpd_dict, file=open("../../../mnt/Archive/Shared/PatentData/SureChemBL/Cpd_Data/cpd_ID_index_dict.p", "wb"))


def build_patent_ID_mapping(updates, fp):
    """ Build a mapping between patents and igraph indicies (0-len(patents))

    Args:
        updates (list): list of months
        fp (string): filepath to GDrive patent info

    Returns:
        None: saves patent indicies to fp+"patent_ID_index_dict.p"
    """
    patents = []
    for update in tqdm(updates):
        if os.path.isfile("../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Patent_Cpd_Edges/patent_cpd_edges_" + update + ".p"):
            patent_edges = pickle.load(
                file=open("../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Patent_Cpd_Edges/patent_cpd_edges_" + update + ".p", "rb"))
            patents.extend(patent_edges.keys())

    unique_patents = list(set(patents))
    print("Patents from patent edges:", len(patents))
    print("Unique patents:", len(unique_patents))

    patent_ids = dict(zip(unique_patents, np.arange(0, len(unique_patents), 1)))

    pickle.dump(patent_ids, file=open("../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/patent_ID_index_dict.p", "wb"))


def replaceIds(updates, fp, cpd_id_dict, patent_id_dict):
    """ Replace SureChemBL ids with igraph indicies - will save igraph memory

    Args:
        updates (list): all months in a certain range (YYYY-MM)
        fp (string): filepath to compound data
        cpd_id_dict (dictionary): links SureChemBL ids to igraph indicies
        patent_id_dict (dictionary): links patent ids to igraph indicies

    Returns:
        none: saves each update to a pickle file
    """
    #number to add to patents to avoid duplicate igraph indicies
    num_cpds = len(cpd_id_dict)
    print("Num Cpds:", num_cpds)
    print("Num patents:", len(patent_id_dict))

    print("Max cpd value:", max(cpd_id_dict.values()))
    print("Max patent value:", max(patent_id_dict.values()))

    for update in tqdm(updates):
        patent_cpd_edges = pickle.load(
            file=open("../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Patent_Cpd_Edges/patent_cpd_edges_" + update + ".p", "rb"))

        patent_id_edges = {}  #New dictionary to hold patent/id relations

        #Track number of compounds that do not appear in SureChemBL compound list
        count = 0
        failed = 0

        #Replace SureChemBL cpd ids with igraph indicies
        for patent, cpds in patent_cpd_edges.items():
            indicies = []
            count += len(cpds)
            for cpd in cpds:
                try:
                    indicies.append(cpd_id_dict[cpd])
                except:
                    failed += 1  #Count failures if compound doesn't appear in cp_id_dict

            #Link patent index with all compound indicies associated with it
            patent_id_edges[patent_id_dict[patent] + num_cpds] = indicies

        #Save each month's edges
        pickle.dump(patent_id_edges,
                    file=open("../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Patent_ID_Edges/patent_id_edges" + update + ".p", "wb"))


def build_cpd_edgelist(updates, fp):
    """ Build unique edgelist for igraph network

    Args:
        updates (list): list of months (YYYY-MM format)
        fp (string): filepath to CpdPatentIdsDates directory

    Returns:
        None: saves edgelist to pickle file (index_edgelist.p)
    """
    edgelist = {}
    for update in updates:
        #Load all patent edges
        patent_index_edges = pickle.load(
            file=open(fp + "patent_id_edges" + update + ".p", "rb"))

        #Loop through all compounds in a particular month
        for cpds in patent_index_edges.values():
            #Only consider patents with more than two compounds
            if len(cpds) > 1:
                #Add combinations to growing edgelist using dictionaries and sets
                for cpd in cpds:
                    if cpd not in edgelist.keys():
                        edgelist[cpd] = set()
                        edgelist[cpd].update([x for x in cpds if x != cpd])
                    else:
                        edgelist[cpd].update([x for x in cpds if x != cpd])

    pickle.dump(edgelist, file=open(fp + "index_edgelist.p", "wb"))

    del (edgelist)


def build_bipartite_edgelist(updates, fp):
    """ Builds patent-cpd edges using igraph indicies

    Args:
        updates (list): list of months (YYYY-MM format)
        fp (string): filepath to CpdPatentIdsDates directory

    Returns:
        None: saves edges to pickle file (index_edgelist_bipartite.p in CpdPatentIdsDates)
    """
    edges = []
    max_value = 0

    for update in updates:
        patent_index_edges = pickle.load(
            file=open("../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Patent_ID_Edges/patent_id_edges" + update + ".p", "rb"))

        for patent, cpds in patent_index_edges.items():
            if patent > max_value:
                max_value = patent

            for cpd in cpds:
                edges.append((patent, cpd))

    pickle.dump(edges,
                file=open(
                    "../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/index_edgelist_bipartite.p",
                    "wb"))
    print("Max Value is:", max_value)

    return edges


def build_full_bipartite_network(edgelist, cpd_id_dict, patent_id_dict):
    """ Builds full igraph network containing patents and compounds

    Args:
        edgelist (list of sets): list of all edges between patent & compound indicies
        cpd_id_dict (dict): links SureChemBL cpd ids with igraph indicies
        patent_id_dict (dict): links patent ids with igraph indicies
    """
    print("Sum of cpd & patent id dicts is:",
          len(cpd_id_dict) + len(patent_id_dict))
    G = ig.Graph()

    #Add nodes
    G.add_vertices(len(cpd_id_dict) + len(patent_id_dict))
    G.vs["name"] = [*cpd_id_dict] + [*patent_id_dict]
    #Type is cpd/patent to distinguish bipartite nature of nodes
    G.vs["type"] = [0] * len(cpd_id_dict) + [1] * len(patent_id_dict)

    #Add edges
    G.add_edges(edgelist)

    print(ig.summary(G))

    del (edgelist)
    del (cpd_id_dict)
    del (patent_id_dict)

    pickle.dump(G, file=open("../../../mnt/Archive/Shared/PatentData/SureChemBL/Graphs/cpd_patent_G.p", "wb"))


def main():
    ### Read in data ###

    # # # #Build list of all unique compounds & patents, as well as dictionaries with dates
    # get_cpd_patent_info(
    #     "~/../../../../Volumes/Macintosh HD 4/SureChemBL/Cpd_Data/")

    ### Create cpd-patent graph ###
    #Note - takes ~90GB and ~20 minutes to build the full network

    ### Test rclone from GDrive ###
    # Move all files to Google Drive

    updates = build_month_list(1976, 2022)

    # for update in ["2021-06"]: #updates:
    #     # ## Moves data to scratch on Agave

    #     # for label in [
    #     #         "unique_cpds", "unique_patents", "cpd_date_dict",
    #     #         "patent_date_dict", "cpd_patent_edges", "patent_cpd_edges"
    #     # ]:
    #     # subprocess.run([
    #     #     "rclone",
    #     #     "copy",
    #     #     "SureChemBL_Patents:CpdPatentIdsDates/" + label + "_" + update +
    #     #     ".p",
    #     #     "/scratch/jmalloy3/CpdPatentIdsDates/",
    #     # ])

    #     print("--- Building:", update, "---")

    #     #Build bipartite network
    #     G = build_bipartite_network(
    #         pickle.load(file=open(
    #             "../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Unique_Cpds/unique_cpds_"
    #             + update + ".p", "rb")),
    #         pickle.load(file=open(
    #             "../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Unique_Patents/unique_patents_"
    #             + update + ".p", "rb")),
    #         pickle.load(file=open(
    #             "../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Cpd_Date_Dict/cpd_date_dict_"
    #             + update + ".p", "rb")),
    #         pickle.load(file=open(
    #             "../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Patent_Date_Dict/patent_date_dict_"
    #             + update + ".p", "rb")),
    #         pickle.load(file=open(
    #             "../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Cpd_Patent_Edges/cpd_patent_edges_"
    #             + update + ".p", "rb")))

    #     pickle.dump(G,
    #                 file=open(
    #                     "../../mnt/Archive/Shared/PatentData/SureChemBL/Graphs/G_" +
    #                     update + ".p", "wb"))

    #     del (G)

    #     #Build cpd-cpd network
    #     G = build_cpd_network(
    #         pickle.load(file=open(
    #             "../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Unique_Cpds/unique_cpds_"
    #             + update + ".p", "rb")),
    #         pickle.load(file=open(
    #             "../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Cpd_Date_Dict/cpd_date_dict_"
    #             + update + ".p", "rb")),
    #         pickle.load(file=open(
    #             "../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Patent_Cpd_Edges/patent_cpd_edges_"
    #             + update + ".p", "rb")))

    #     # pickle.dump(G, file=open("/scratch/jmalloy3/Graphs/G_cpd_" + update + ".p", "wb")) ## Too much memory
    #     G.save("../../mnt/Archive/Shared/PatentData/SureChemBL/Graphs/G_cpd_" + update +
    #            ".gmlz",
    #            format="graphmlz")  #save in zipped-gml format to save memory

    #     del (G)  #remove G from memory to free up space

    # #Move all Graphs to GDrive
    # subprocess.run([
    #     "rclone", "moveto", "/scratch/jmalloy3/Graphs",
    #     "SureChemBL_Patents:Graphs"
    # ])

    # #Delete all files from scratch (already backed up in GDrive)
    # subprocess.run(["rm", "-r", "/scratch/jmalloy3/CpdPatentIdsDates/"])

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

    ### Build Full igraph Network  ###
    
    # Step 0.5: Map SureChemBL cpd ids to igraph vertex labels - SHOULD ONLY BE RUN ONCE
    fp = "../../../mnt/Archive/Shared/PatentData/SureChemBL/Cpd_Data/"
    build_cpd_ID_mapping(fp)

    # # Load cpd-id dictionary
    # cpd_id_dict = pickle.load(file=open(fp + "cpd_ID_index_dict.p", "rb"))

    # # Make dictionary of patent-igraph indicies - SHOULD ONLY BE RUN ONCE
    # replaceIds(updates,
    #            "../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Patent_Cpd_Edges/",
    #            cpd_id_dict)

    # build_edgelist(updates, fp)  #NOTE: runs out of memory on Agave :(

    ## Build Full igraph cpd-patent network ###
    fp = "../../mnt/Archive/Shared/PatentData/SureChemBL/Cpd_Data/"

    #Step 1: build patent-id dictionary - should only be run once
    build_patent_ID_mapping(updates, fp)

    #Step 2: Update patent-cpd-id files to include patent ids
    #Load cpd-id dictionary - should only be run once
    cpd_id_dict = pickle.load(file=open(
        "../../mnt/Archive/Shared/PatentData/SureChemBL/Cpd_Data/cpd_ID_index_dict.p",
        "rb"))
    patent_id_dict = pickle.load(file=open(
        "../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/patent_ID_index_dict.p",
        "rb"))

    print("Num Cpds:", len(cpd_id_dict))
    print("Patents", len(patent_id_dict))

    replaceIds(updates, fp, cpd_id_dict, patent_id_dict)

    # #Step 3: Make edgelist of patent-cpd edges, using igraph ids - should only be run once
    edgelist = build_bipartite_edgelist(updates, fp)

    # Step 4: Build & save full igraph network
    # edgelist = pickle.load(
    #     file=open("/scratch/jmalloy3/Patents/index_edgelist_bipartite.p", "rb"))

    # #GDrive filepath
    # cpd_id_dict = pickle.load(file=open("G:/Shared drives/SureChemBL_Patents/Cpd_Data/cpd_ID_index_dict.p", "rb"))
    # patent_id_dict = pickle.load(file=open("G:/Shared drives/SureChemBL_Patents/CpdPatentIdsDates/patent_ID_index_dict.p", "rb"))

    print("Num cpds:", len(cpd_id_dict))
    print("Num patents:", len(patent_id_dict))

    build_full_bipartite_network(edgelist, cpd_id_dict, patent_id_dict)
    #
    # Step 5: Add cpd names & patent ids


if __name__ == "__main__":
    main()
