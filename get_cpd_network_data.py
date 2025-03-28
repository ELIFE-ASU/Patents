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
from itertools import accumulate
import os
import subprocess
import pandas as pd


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
    Assumes that the data are stored in the /Volumes/Macintosh HD 4/SureChemBL/Graphs/G_Cpd directory and in the format
    G_cpd_<update>.gmlz.

    Args:
        fp (string): filepath to a specific cpd-cpd graph

    Returns:
        G (igraph object): cpd-cpd network pertaining to a specific update
    """
    fp = "../../mnt/Archive/Shared/PatentData/SureChemBL/Graphs/G_Cpd/G_cpd_" + update + ".gmlz"
    G = ig.Graph.Read_GraphMLz(fp)
    print("Loaded graph:", update)
    print(ig.summary(G))

    # #One time function - to save all gmlz files as pickles to save time
    # pickle.dump(G, file=open("Data/Graphs/G_cpd_" + update + ".p", "wb"))
    # print("Dumped G")

    # G = pickle.load(file=open("Data/Graphs/G_cpd_" + update + ".p", "rb"))
    # print(ig.summary(G))

    return G


def get_network_stats(start, stop):
    """Finds basic network statistics SureChemBL cpd-cpd graphs in a given range

    Calculates num nodes, num edges, avg degree, max degree,
    avg clustering coefficient, largest connected component size

    Args:
        start (int): year of starting point for analysis
        end (int): year of ending point (inclusive)

    Returns:
        (none): writes a file containing the basic network statistics for each month
                in a given myself
    """
    updates = build_month_list(start, stop)
    data = []

    for update in updates:
        # subprocess.run([
        #     "rclone",
        #     "copy",
        #     "SureChemBL_Patents:Graphs/G_cpd_" + update + ".p",
        #     "/scratch/jmalloy3/Graphs/",
        # ])

        G = read_cpdcpd_graph(update)
        network_stats = {}

        degrees = get_degrees(G)
        id_degrees = get_id_degree(G)
        # del (G)
        network_stats["Nodes"] = G.vcount()
        network_stats["Edges"] = G.ecount()
        network_stats["Avg Degree"] = np.mean(degrees)
        network_stats["Max Degree"] = max(degrees)
        network_stats["LCC Size"] = G.clusters().giant().vcount()
        network_stats["Clustering coefficient"] = G.transitivity_undirected()
        print(network_stats)
        print()
        data.append(network_stats)

        #Writing every month of stastistics, degrees, and ID degrees
        df = pd.DataFrame(data)
        pickle.dump(df,
                file=open(
                    "../../mnt/Archive/Shared/PatentData/SureChemBL/NetworkStats/stats_" + str(start) + "_" +
                    str(stop) + ".p", "wb"))

        pickle.dump(degrees,
                    file=open(
                        "../../mnt/Archive/Shared/PatentData/SureChemBL/Degrees/Months/Degrees/degrees_" + update +
                        ".p", "wb"))

        pickle.dump(id_degrees,
                    file=open(
                        "../../mnt/Archive/Shared/PatentData/SureChemBL/Degrees/Months/IdDegrees/id_degrees_" +
                        update + ".p", "wb"))

    df = pd.DataFrame(data)
    pickle.dump(df,
                file=open(
                    "../../mnt/Archive/Shared/PatentData/SureChemBL/NetworkStats/stats_" + str(start) + "_" +
                    str(stop) + ".p", "wb"))


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


def get_degree_distributions():
    """ Finds and saves degrees across SureChemBL update graphs

    Collects degrees in one place (all_degrees, list of degree distributions)
    and finds the average degree distribution (avg, single list). Saves both
    to the 'Data/Degrees/' directory

    """
    #Load all degrees associated with quarterly updates
    print("\n----- Loading Degree Distributions -----\n")
    all_degrees = []
    for f in os.listdir("Data/Degrees/"):
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


def link_id_degrees(full_id_degrees, id_degrees, i, bins):
    """ Puts the ids of specific updates/times into a full dictionary of id-[degrees]

    Takes in a specific update (id_degrees), and adds the degrees present in that
    update to the full id:degree dictionary. The degrees are a list of len(bins),
    and each item is the sum of the previous item and the total provided in id_degrees.

    Args:
        full_id_degrees (dictionary): links a SureChemBL id with a list of degrees
            over time. The degrees are added together as they increase.
        id_degrees (dictionary): links a SureChemBL id with a specific degree from
            a specific update.
        i (int): the increment of the specific update, provides the necessary place
            for the full_id_degrees dictionary
        bins (int): number of updates

    Returns:
        full_id_degrees: updated full dictionary of id:[degrees]
    """
    for key, value in id_degrees.items():
        #If a id is in the dictionary, update the appropriate degree value
        if key in full_id_degrees:
            full_id_degrees[key][i] = value
        # If a id is not in the dictionary, add it with a list of len(updates) 0s and update the appropriate degree value
        else:
            full_id_degrees[key] = [0] * bins
            full_id_degrees[key][i] = value

    return full_id_degrees


def replace_zeroes(full_id_degrees):
    """Replaces all zero values using itertools.accumulate()

    Takea the full dictionary of id:[degree], then finds the cumulative sum
    across the entire list of degrees and replaces the original degree list with
    this cumulative list (for preferential attachment purposes)

    Args:
        full_id_degrees (dictionary): links a SureChemBL id with a list of degrees
            over time. The degrees are added together as they increase

    Returns:
        full_id_degrees: full dictionary with cumulative sums as degree list
    """
    #Finds cumulative sum of values
    print("\n----- Replacing zero values -----\n")
    for key, value in full_id_degrees.items():
        cum_value = list(accumulate(value))
        full_id_degrees[key] = cum_value

    return full_id_degrees


def pref_attachment_calculation(full_id_degrees):
    """ Calculates preferential attachment over a list of degrees

    Given a list of increasing degrees associated with each SureChemBL id, the
    preferential attachment index sum(value[n+1] - value[n]), over all n is calculated

    Args:
        full_id_degrees (dictionary): links a SureChemBL id with a list of degrees
            over time. The degrees are added together as they increase

    Returns:
        pref_attach_dict: associated a preferential attachment index with each
            SureChemBL ID
    """

    print("\n----- Calculating preferential attachment -----\n")
    pref_attach_dict = {}
    for key, value in full_id_degrees.items():
        c = 0
        attachments = []
        while c < len(value) - 1:
            attachments.append(value[c + 1] - value[c])
            c += 1
        pref_attach_dict[key] = np.mean(attachments)

    return pref_attach_dict


def calculate_preferential_attachment(start, stop):
    """ Calculates preferential attachment index (see Rednar 2004) for SureChemBL degrees
    across all patents

    Uses data stored in the 'Data/Degrees/id_degrees_*' files to build a preferential
    attachment index. Saves both the full list of id:degree pairs (summed over time)
    and the preferential attachment indicies of each degree to pickle files in
    'Data/Degrees' directory
    """
    print("\n----- Building id-degree dictionary -----\n")
    updates = build_month_list(start, stop)

    full_id_degrees = {}
    bins = len(updates)

    i = 0  #update number (to link id_degree lists with updates)
    for update in updates:
        # subprocess.run([
        #     "rclone",
        #     "copy",
        #     "SureChemBL_Patents:Degrees/Months/id_degrees_" + update + ".p",
        #     "/scratch/jmalloy3/Degrees/Months/",
        # ])

        #Load id_degree dictionary - attempting through GDrive
        id_degrees = pickle.load(file=open(
            "SureChemBL_Patents:Degrees/Months/id_degrees_" + update + ".p",
            "rb"))

        # #Load id_degree dictionary
        # id_degrees = pickle.load(file=open(
        #     "/scratch/jmalloy3/Degrees/Months/id_degrees_" + update +
        #     ".p", "rb"))

        full_id_degrees = link_id_degrees(full_id_degrees, id_degrees, i, bins)
        i += 1  #increment position

        # #Remove file from scratch
        # subprocess.run([
        #     "rclone", "moveto",
        #     "/scratch/jmalloy3/Degrees/Months/id_degrees_" + update + ".p",
        #     "SureChemBL_Patents:Degrees/Months/id_degrees_" + update + ".p",
        # ])

    #print(list(islice(full_id_degrees.items(), 10)))

    full_id_degrees = replace_zeroes(full_id_degrees)

    #print(list(islice(full_id_degrees.items(), 10)))

    pickle.dump(full_id_degrees,
                file=open(
                    "/scratch/jmalloy3/Degrees/full_id_degrees_" + str(start) +
                    "_" + str(stop) + ".p", "wb"))

    pref_attach_dict = pref_attachment_calculation(full_id_degrees)

    #print(list(islice(pref_attach_dict.items(), 10)))

    pickle.dump(pref_attach_dict,
                file=open(
                    "/scratch/jmalloy3/pref_attach_dict_" + str(start) + "_" +
                    str(stop) + ".p", "wb"))


def clear_scratch(start, stop):
    """ Move files from scratch to GDrive

    Move: id_degrees (in Degrees/Months) - worked
        : degrees (in Degrees/Months) - worked
        : network_stats (in NetworkStats)
        : full_id_degrees (in Degrees)
        : pref_attach_dict (in scratch/jmalloy3)

    Delete: G_cpd_XXXX-MM.p, in Graphs/

    """
    fps_months = [
        "Degrees/Months/id_degrees_", "Degrees/Months/degrees_", "Graphs/C_cpd_"
    ]
    fps_years = [
        "NetworkStats/stats_", "Degrees/full_id_degrees_",
        "pref_attach_dict_"
    ]

    updates = build_month_list(start, stop)
    for update in updates:
        for f in fps_months:
            #Move all files to GDrive
            subprocess.run([
                "rclone", "moveto", "/scratch/jmalloy3/" + f + update + ".p",
                "SureChemBL_Patents:" + f + update + ".p"
            ])

    for f in fps_years:
        #Move all files to GDrive
        subprocess.run([
            "rclone", "moveto",
            "/scratch/jmalloy3/" + f + str(start) + "_" + str(stop) + ".p",
            "SureChemBL_Patents:" + f + str(start) + "_" + str(stop) +
            ".p"
        ])


def build_increments(start, stop, increment):
    """ Build a range of years from start to stop, given a specific increment

    Args:
        start (int): starting year
        stop (int): ending year
        increment (int): number of years in one time series

    Returns:
        year_increments (list): list of tuples, each containing one time series
    """
    year_increments = []
    start_increment = start
    stop_increment = start + increment - 1 #keep increment-range (not increment + 1)
    while stop_increment <= stop:
        year_increments.append((start_increment, stop_increment))
        start_increment += increment
        stop_increment += increment

    return year_increments

def main():
    start = 2021
    stop = 2022

    ## NOTE: building preferential attachement across entire network

    # five_year_increments = build_increments(start, stop, 5)
    # ten_year_increments = build_increments(start, stop, 10)
    # #twenty_year_increments = build_increments(start, stop, 20)

    # for inc in five_year_increments + ten_year_increments: # + twenty_year_increments:
    #     start, stop = inc[0], inc[1]

    #Calculate basic high-level network stats from SureChemBL updates
    get_network_stats(start, stop)

    start = 1976
    stop = 1979
    get_network_stats(start, stop)


    # # #Store all degree distributions in a single list
    # # #get_degree_distributions()

    # #Calculate preferential attachment index for a range of the SureChemBL dataset
    # calculate_preferential_attachment(start, stop)

    #clear_scratch(start, stop)


if __name__ == "__main__":
    main()
