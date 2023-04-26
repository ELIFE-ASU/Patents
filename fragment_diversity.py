import igraph as ig
import ast
import re
import itertools
from tqdm import tqdm
import os
import pickle
import pandas as pd
from random import sample


def get_list(line):
    """ Evaluation a line to find the description and list containing data
        for eventual addition to an igraph graph

    Args:
        line (str): specific line to evaluate

    Returns:
        list: 2-element list of [0]: description (vertices, edges, vertexcolors, or edgecolors) 
            and [1]: list of values matching that description
    """
    #Split each line on the end of the first word
    line = line.split("s ")
    #values = json.loads(line[1].replace(" ", ", "))
    try:
        #Lists of ints (vertices & edges) can be treated with literal_eval
        values = ast.literal_eval(line[1].replace(" ", ", "))
    except ValueError:
        #Lists of strings/chars need to be treated with regex (for some reason)
        values = re.sub('[\[\]]', '', line[1])
        values = values.split()

    return [line[0] + "s", values]


def build_fragment(start, end, lines, vscolor_map, escolor_map):
    """ Builds an igraph graph of a single fragment

    Args:
        start (int): line number of the vertices entry
        end (int): line number of the edgecolours entry
        lines (list): list of all text file lines

    Returns:
        igraph object: graph containing vertex, edges, and vs/es colors
    """
    g = ig.Graph()

    #Label map from MA node labels to igraph node labels
    node_map = {}

    for line in lines[start:end]:
        #Isolate description
        desc = get_list(line)

        #Iterativley(?) build igraph fragment
        #Add vertices to igraph (guarenteed to be first)
        if desc[0] == "Vertices":
            g.add_vertices(len(desc[1]))

            #Create a map for vertex numbers to igraph values (for edges)
            for i in range(len(desc[1])):
                node_map[desc[1][i]] = i

        #Add edges based on label map built in vetex
        elif desc[0] == "Edges":
            for edge in desc[1]:
                g.add_edges([(node_map[edge[0]], node_map[edge[1]])])

        #Add vertex & edge colors
        elif desc[0] == "VertexColours":
            for element in desc[1]:
                if element not in vscolor_map:
                    vscolor_map[element] = len(vscolor_map) + 1

            g.vs["color"] = [vscolor_map[x] for x in desc[1]]

        elif desc[0] == "EdgeColours":
            for element in desc[1]:
                if element not in escolor_map:
                    escolor_map[element] = len(escolor_map) + 1

            g.es["color"] = [escolor_map[x] for x in desc[1]]

    return g


def check_iso(candidate_fragments, all_frags, frag_count):
    """ Finds the unique fragments within a set of candidates

    Args:
        candidate_fragments (list): list of igraph graphs from a single text file
        all_frags (list): list of all unique igraph graphs from previous files

    Returns:
        list: updated list of unique igraph graphs
    """
    #Check all candidate fragments...
    for possible_f in candidate_fragments:
        is_unique = True
        #...against all unique fragments
        for i in range(len(all_frags)):
            if possible_f.isomorphic_vf2(all_frags[i],
                                         color1=possible_f.vs["color"],
                                         color2=all_frags[i].vs["color"],
                                         edge_color1=possible_f.es["color"],
                                         edge_color2=all_frags[i].es["color"]):
                #if candidate is isomorphic, it is not unique, therefore break off check
                is_unique = False

            if not is_unique:
                frag_count[i] += 1
                break

        #If it makes it through the full list of fragments, it is unique!
        if is_unique:
            all_frags.append(possible_f)
            frag_count[len(all_frags) - 1] = 1  #Start count at 1

    #Return updated full list
    return all_frags, frag_count


def build_month_increments(start, stop):
    """ Build all monthly increments from the start year to stop year in the
    format YEAR-MONTH

    Args:
        start (int): start year of increments
        stop (int): end year of increments

    Returns:
        list: list of strings holding the YEAR-MONTH increments
    """
    months = []
    while start <= stop:
        for month in [
                "01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
                "11", "12"
        ]:
            months.append(str(start) + "-" + month)
        start += 1

    return months


def get_files(fp, labels):
    """ Find all file paths to specific log files to analyze

    These all correspond to fully completed MAs

    Args:
        fp (str): filepath to overall AssemblyValue directory
        labels (list): all labels of completed MAs

    Returns:
        (list): list of all filepaths
    """
    possible_dataSources = ["AuthorCpds_Done/", "AssigneeCpds_Done/"]

    fps = []
    for label in labels:
        if "SCHEMBL" in label:
            for dataSource in possible_dataSources:
                if os.path.exists(fp + dataSource + label + ".txt"):
                    fps.append(fp + dataSource + label + ".txt")
                    break
        elif label.endswith("_new"):
            fps.append(fp + "NewDatabase_Done/"+ label[:-4] +
                       ".txt")
        elif label.endswith("_full"):
            fps.append(fp + "FullDatabase_Done/" + label[:-4] +
                       ".txt")

    return fps


def main():
    ### Using Author Cpds Done to build fragments

    # fp = "Data/AssemblyValues/NewDatabase_Done/" #old newdatabase code

    MA_df_completed = pd.read_csv(
        "Data/AssemblyValues/ALLSAMPLEDcpds_AssemblyGo_COMPLETED.csv")

    print(MA_df_completed.head())

    print(type(list(MA_df_completed["earliest_date"])[0]))

    month = "1980-01"
    month_df = MA_df_completed[MA_df_completed["earliest_date"].str.startswith(
        month)].sample(1000)

    print(month_df)

    files = get_files("Data/AssemblyValues/", list(month_df["label"]))
    print(files)

    # labels = []
    # for index, row in tqdm(MA_df_completed.iterrows(),
    #                        total=len(MA_df_completed)):
    #     labels.append(get_labels(row["label"], row["earliest_date"]))

    # print(labels[0:10])

    # fp = "Data/AssemblyValues/"

    # # months = build_month_increments(1980, 2020) #old newdatabase code

    # df = pd.read_csv("Data/ID_months.csv")

    # all_IDs = {}

    # # for month in months:
    # for index, row in df.iterrows(): #, total=df.shape[0]):
    #     # print("--- Analyzing", row["month"], "---")

    #     # #Old newdatabase code
    #     # files = [x for x in os.listdir(fp) if x.startswith(month)]
    #     # files = [x for x in files if x.endswith(".txt")]
    #     # print(len(files))

    #     fragments = []
    #     all_frags = []
    #     frag_count = {}
    #     vscolor_map = {}
    #     escolor_map = {}

    #     #Get IDs
    #     IDs = ast.literal_eval(row["IDs"])

    #     #Sample 1000 per month (to make this computationally tractable)
    #     if len(IDs) > 1000:
    #         IDs = sample(IDs, 1000)

    #     #Store ids (for agave parallelization)
    #     all_IDs[row["month"]] = IDs

    #     for f in tqdm(IDs):
    #         if os.path.exists(fp + f + ".txt"):
    #             with open(fp + f + ".txt") as f:
    #                 lines = f.readlines()
    #                 lines = [l.strip() for l in lines]

    #         else:
    #             pass

    #         #Isolate fragments lines
    #         for i in range(len(lines)):
    #             if lines[i] == "======":
    #                 if lines[i + 1].startswith("Vertices"):
    #                     #i+1 is the start of the fragment definition, i+5 is the end
    #                     fragments.append(
    #                         build_fragment(i + 1, i + 5, lines, vscolor_map,
    #                                     escolor_map))

    #         ## More testing - check isomorphism within a single output file
    #         all_frags, frag_count = check_iso(fragments, all_frags, frag_count)

    #     pickle.dump(all_frags, file=open("Data/AssemblyValues/Fragments/authorFrags_" + row["month"] + "_updated.p", "wb"))
    #     pickle.dump(frag_count, file=open("Data/AssemblyValues/Fragments/authorFragsCount_" + row["month"] + "_updated.p", "wb"))

    # pickle.dump(all_IDs, file=open("Data/AssemblyValues/Fragments/all_ids_updated.p", "wb"))
    # pickle.dump(vscolor_map, file=open("Data/AssemblyValues/Fragments/vscolor_map.p", "rb"))
    # pickle.dump(escolor_map, file=open("Data/AssemblyValues/Fragments/vscolor_map.p", "rb"))

    # print(vscolor_map)
    # print(escolor_map)
    # print("Unique fragment size:", len(all_frags))
    # # for frag in all_frags:
    # #     print(frag.vs["color"], frag.es["color"])
    # #     print("----")


if __name__ == "__main__":
    main()
