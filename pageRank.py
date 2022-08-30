import igraph as ig
import pickle
import subprocess
from functools import reduce
import pandas as pd


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
    """ Calculates the pagerank of the graph of a given month

    Args:
        G (igraph Graph): graph of a specific month
        month (str): month (to be used for file IO labeling)
    """
    #Calculate pagerank
    G.vs["pagerank"] = G.pagerank()

    #Zip ChEMBL ids and pagerank values (don't care about the graph structure)
    pickle.dump(list(zip(G.vs["name"], G.vs["pagerank"])),
                file=open("/scratch/jmalloy3/pagerank_" + month + ".p", "wb"))

    #Ensure pagerank analysis ends up in Google Drive
    subprocess.run([
        "rclone", "copy", "/scratch/jmalloy3/pagerank_" + month + ".p",
        "SureChemBL_Patents:PageRank/"
    ])

    #Remove files from scratch to keep it empty-ish
    subprocess.run(["rm", "/scratch/jmalloy3/pagerank_" + month + ".p"])
    subprocess.run(["rm", "/scratch/jmalloy3/G_cpd_" + month + ".p"])


def tidy_pageRank(month):
    """ Tidyverse pagerank data

    Args:
        month (str): month of specific analysis

    Returns:
        pandas dataframe: dataframe containing SureChemBL id and pagerank value, for given month
    """
    fp = "SureChemBL_Patents:PageRank/"

    #Load pagerank results from GDrive
    subprocess.run([
        "rclone",
        "copy",
        fp + "pagerank_" + month + ".p",
        "~/Patents/PageRank/",
    ])

    #create dataframe from pagerank results, with columns "ID" and month
    return pd.DataFrame(pickle.load(
        file=open("~/Patents/PageRank/pagerank_" + month + ".p", "rb")),
                        columns=["ID", str(month)])


def create_dataframe(dataframes):
    """ Build an overarching dataframe of all pagerank values

    In the form ID: value1, value2, value3... where valueN is each month value

    Args:
        dataframes (list): list of dataframes containing ID and monthly pagerank values

    Returns:    
        None (writes full dataframe & dropped nan dataframe to /scratch)
    """
    #Merges dataframe
    df = reduce(
        lambda left, right: pd.merge(left, right, on=["ID"], how="outer"),
        dataframes)
    df.to_csv("/scratch/jmalloy3/pagerank_analysis.csv")

    #Drops nans, so only compounds present in every month appear (may not be useful, but could be?)
    df = df.dropna()
    df.to_csv("scratch/jmalloy3/pagerank_lastingCpds.csv")


def main():
    #Calculate pagerank for all months in dataset
    months = build_month_list(1980, 2019)

    # for month in months:
    #     G = load_graph(month)

    #     get_pageRank(G, month)

    #Tidy data from pagerank
    dataframes = []
    for month in months:
        dataframes.append(tidy_pageRank(month))

    create_dataframe(dataframes)


if __name__ == "__main__":
    main()
