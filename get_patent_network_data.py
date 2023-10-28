import igraph as ig
import pickle
import numpy as np
import os
import re
import pandas as pd

def get_graph_fps():
    return ["../../../../mnt/Archive/Shared/PatentData/SureChemBL/Graphs/"+fp for fp in os.listdir("../../../../mnt/Archive/Shared/PatentData/SureChemBL/Graphs") if re.match(r'^cpd_patent_[0-9]', fp)]


def get_patent_stats(fps):
    total_patents = []
    new_patents = []
    avg_degree = []
    months = []

    for f in fps:
        print(f"Analyzing {f}...")

        months.append(f.split("_")[-1].split(".")[0])

        G = pickle.load(open(f, "rb"))
        G_patents = G.vs(type=1)

        total_patents.append(len(G_patents.select(_degree_gt=0)))

        if new_patents:
            new_patents.append(total_patents[-1] - total_patents[-2])
        else:
            new_patents.append(total_patents[-1])

        avg_degree.append(np.mean(G.degree(G_patents.select(_degree_gt=0))))

    return months, total_patents, new_patents, avg_degree

def save_dataframe(months, total_patents, new_patents, avg_degree):
    df = pd.DataFrame(
    {'month': months,
     'total_patents': total_patents,
     'new_patents': new_patents,
     'avg_degree': avg_degree
    })
    df.to_csv("../../../../mnt/Archive/Shared/PatentData/SureChemBL/NetworkStats/patent_stats_updated.csv")

def main():
    graphs_fps = get_graph_fps()

    months, total_patents, new_patents, avg_degree = get_patent_stats(graph_fps)

    save_dataframe(months, total_patents, new_patents, avg_degree)

if __name__=="__main__":
    main()