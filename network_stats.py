import time
import igraph as ig
import pickle
import numpy as np

#TODO (JFM): add toplevel comments


def nodes_before_date(G, cutoff_date):
    """
    Finds all nodes that exist in graph G before a given date

    Analyzes the "date" parameter of igraph network G, and returns the vertex 
    labels of all nodes which were introduced to the network before a specific 
    date

    Args:
        G: igraph network G, has "id", "date", "type" field
        date: a specific date, in YYYY-MM-DD format

    Returns:
        A list of node labels
    """
    G.vs["dt"] = [time.strptime(d, "%Y-%m-%d") for d in G.vs["date"]]
    print(ig.summary(G))

    early_nodes = G.vs.select(dt_lt=cutoff_date)
    return [v["name"] for v in early_nodes if v["type"] == 0]


def get_attachment_index(G, G_new, name_list):
    """ Calculate the attachment index for a given list of cpd names

    The degrees of the given nodes are calculated for G & G_new (G + 1 timestep),
    and the attachment index for each nodes is calcualted from the Barabasi textbook 
    formula (degree at timestep t+1 - degree at timestep t)

    Args:
        G: igraph network at timestep t
        G_new: igraph network at timestep t+1
        name_list: list of cpd names to be calculated

    Returns:
        Attachment index for all nodes
    """
    attachment_indexes = []
    for n in name_list:
        try:
            attachment_indexes.append(G_new.degree(n) - G.degree(n))
        except:
            print()

    print(attachment_indexes[0:100])
    print(len(attachment_indexes))
    print(np.mean(attachment_indexes))

    #TODO (JFM): finish calculating attachment index


def main():
    ### Find degrees of compounds before a specific date ###

    #Test graph
    G = pickle.load(file=open("Data/Graphs/G_20201001.p", "rb"))
    d = "2020-11-01"
    cpd_names = nodes_before_date(G, time.strptime(d, "%Y-%m-%d"))

    G_new = pickle.load(file=open("Data/Graphs/G_20210101.p", "rb"))

    ### Calculate attachment index ###
    get_attachment_index(G, G_new, cpd_names)

    #TODO (JFM): analysis of attachment index (linear scaling?, etc...)


if __name__ == "__main__":
    main()
