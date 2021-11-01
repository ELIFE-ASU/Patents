import igraph as ig
import numpy as np
import pickle
from random import sample


def read_month(fp):
    """ Read in a specific month's igraph network

    Args:
        fp (string): filepath to igraph network

    Returns:
        igraph network: contains all compounds before or in a specific month, along with all patents
    """
    return pickle.load(file=open(fp, "rb"))


def main():
    #Testing - read in a compound in 1980-12, and track degrees over 1980
    fp = "/scratch/jmalloy3/Patents/Graphs/cpd_patent_1980-12.p"
    G = read_month(fp)

    num_cpds = 21641384  #Numbers from cpd/patent id dicts used to build igraph network
    num_patents = 4578946

    seq = G.vs.select(lambda v: v.index < num_cpds)
    print([v.index for v in seq])
    print(sample([v.index for v in seq], 1))


if __name__ == "__main__":
    main()