import pandas as pd
import os
from tqdm import tqdm


def get_MA(fp):
    """ Get AssemblyGo MA value

    Args:
        fp (str): filepath to a particular .txt AssemblyGo output file

    Returns:
        label (str): label of the compound which was analzyed (empty if a failure)
        MA (int): assemblyGo MA value (-1 if a failure)
    """
    with open(fp) as f:
        lines = f.readlines()

    try:
        #molfile will be the last element in 0th line
        label = lines[0].split()[-1].split("/")[-1].split(".")[0]

        #MA will be last elemnt in -2nd line (will be an int)
        MA = int(lines[-2].split()[-1])

        return label, MA
    except:
        return "", -1


def main():
    ### Read in .txt files, return the label & MA value

    #Create & save one dataframe per database, eventually to merge with data validation csvs in analyze_ma_dist.ipynb
    database_fps = {
        "NewDatabase": "NewDatabase_Done",
        "FullDatabase": "FullDatabase_Done",
    }

    for database, name in database_fps.items():
        print("----- Analzying", database, "-----")
        MA_values = []
        fp = "Data/AssemblyValues/" + name + "/"
        for file in os.listdir(fp):
            if file.endswith(".txt"):
                label, MA = get_MA(fp + file)
                MA_values.append({
                    "label": label,
                    "MA_assemblyGo": int(MA)
                })

        MA_df = pd.DataFrame(MA_values)
        MA_df.to_csv("Data/AssemblyValues/" + database + "_AssemblyGo.csv")


if __name__ == "__main__":
    main()
