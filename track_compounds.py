import numpy as np
import pickle
import igraph as ig


def main():
    #Testing - load id-degrees from a given month
    # id_deg = pickle.load(file=open(
    #     "G:\\Shared drives\\SureChemBL_Patents\\Degrees\\Months\\id_degrees_1980-01.p",
    #     "rb"))
    # id_deg2 = pickle.load(file=open(
    #     "G:\\Shared drives\\SureChemBL_Patents\\Degrees\\Months\\id_degrees_1980-02.p",
    #     "rb"))
    # print(id_deg)

    # full_id_deg = pickle.load(file=open(
    #     "G:\\Shared drives\\SureChemBL_Patents\\Degrees\\full_id_degrees_1980_1984.p",
    #     "rb"))

    # full_id_deg2 = pickle.load(file=open(
    #     "G:\\Shared drives\\SureChemBL_Patents\\Degrees\\full_id_degrees_1985_1989.p",
    #     "rb"))

    # test_id = "SCHEMBL21799"

    # print(full_id_deg[test_id])
    # print(len(full_id_deg[test_id]))
    # print()
    # print(full_id_deg2[test_id])
    # print(len(full_id_deg2[test_id]))

    ##BASIC TESTING NOW
    cpd_patent_edges = pickle.load(file=open(
        "G:\\Shared drives\\SureChemBL_Patents\\CpdPatentIdsDates\\cpd_patent_edges_1974-01.p",
        "rb"))
    print(cpd_patent_edges)
    print()
    patent_cpd_edges = pickle.load(file=open(
        "G:\\Shared drives\\SureChemBL_Patents\\CpdPatentIdsDates\\patent_cpd_edges_1974-01.p",
        "rb"))
    print(patent_cpd_edges)

    # G1 = pickle.load(file=open(
    #     "G:\\Shared drives\\SureChemBL_Patents\\Graphs\\G_cpd_1973-01.p",
    #     "rb"))

    # G2 = pickle.load(file=open(
    #     "G:\\Shared drives\\SureChemBL_Patents\\Graphs\\G_cpd_1973-02.p",
    #     "rb"))

    # G1_edges = G1.vs["name"]
    # G2_edges = G2.vs["name"]

    # all_cpds = list(set(G1_edges + G2_edges))

    # edge_dict = {}
    # for cpd in G1_edges:
    #     print(G1.vs.find(cpd))


if __name__ == "__main__":
    main()
