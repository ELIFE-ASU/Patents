import igraph as ig
import pickle
import numpy as np
import os
import re
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


def get_pickled_files(month):
    try:
        patent_date_dict = pickle.load(
            open(
                "../../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Patent_Date_Dict/patent_date_dict_"
                + month + ".p", "rb"))
        patent_cpd_edges = pickle.load(
            open(
                "../../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Patent_Cpd_Edges/patent_cpd_edges_"
                + month + ".p", "rb"))
    except EOFError:
        print(f"{month} has corrupted files")
        patent_date_dict = None
        patent_cpd_edges = None

    return patent_date_dict, patent_cpd_edges


def get_patent_stats():
    total_patents = []
    new_patents = []
    avg_degree = []
    months = []

    for month in build_month_list(1976, 2022):
        print(f"Analyzing {month}...")

        months.append(month)

        patent_date_dict, patent_cpd_edges = get_pickled_files(month)

        if patent_date_dict:
            n_patents = len(patent_date_dict.keys())
            new_patents.append(n_patents)

            if total_patents:
                total_patents.append(n_patents + total_patents[-1])
            else:
                total_patents.append(n_patents)
        else:
            new_patents.append(np.nan)
            total_patents.append(np.nan)

        if patent_cpd_edges:
            degrees = []
            for _, cpds in patent_cpd_edges.items():
                degrees.append(len(cpds))
            avg_degree.append(np.mean(degrees))
        else:
            degrees.append(np.nan)
            avg_degree.append(np.nan)

    return months, total_patents, new_patents, avg_degree


def save_dataframe(months, total_patents, new_patents, avg_degree):
    df = pd.DataFrame({
        'month': months,
        'total_patents': total_patents,
        'new_patents': new_patents,
        'avg_degree': avg_degree
    })
    df.to_csv(
        "../../../../mnt/Archive/Shared/PatentData/SureChemBL/NetworkStats/patent_stats_updatedv2.csv"
    )


def main():
    months, total_patents, new_patents, avg_degree = get_patent_stats()

    save_dataframe(months, total_patents, new_patents, avg_degree)


if __name__ == "__main__":
    main()
