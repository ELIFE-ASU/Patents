import igraph as ig
import pickle
import subprocess


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


def load_graph(month):
    """ Loads a igraph graph (stored in pickle format)

    Args:
        fp (string): filepath to igraph graph

    Returns:
        igraph graph: Graph of 1 month of SureChemBL
    """
    subprocess.run([
        "rclone",
        "copy",
        "SureChemBL_Patents:Graphs/G_cpd_" + month + ".p",
        "/scratch/jmalloy3/",
    ])
    return pickle.load(file=open("/scratch/jmalloy3/G_cpd_" + month +
                                 ".p", "rb"))


def get_pageRank(G, month):
    G.vs["pagerank"] = G.pagerank()

    #Zip ChEMBL ids and pagerank values (don't care about the graph structure)
    pickle.dump(list(zip(G.vs["name"], G.vs["pagerank"])),
                file=open("/scratch/jmalloy3/pagerank_" + month + ".p", "wb"))
    subprocess.run([
        "rclone", "copy", "/scratch/jmalloy3/pagerank_" + month + ".p",
        "SureChemBL_Patents:PageRank/"
    ])
    subprocess.run(["rm", "/scratch/jmalloy3/pagerank_" + month + ".p"])
    subprocess.run(["rm", "/scratch/jmalloy3/G_cpd_" + month + ".p"])


def main():
    #Find pagerank for all months in dataset
    months = build_month_list(1980, 2019)

    for month in months:
        G = load_graph(month)

        get_pageRank(G, month)


if __name__ == "__main__":
    main()
