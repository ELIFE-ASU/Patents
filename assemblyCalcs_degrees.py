# import assemblycalculator as ac
# import multiprocessing as mp
import pickle
import pandas as pd
from itertools import islice
import numpy as np
from random import sample


def get_data(fp):
    """ Read in id/degree dictionary

    Contains a SureChemBL ID and the corresponding degree for a specific month

    Args:
        fp (str): filepath to id_degree pickle file

    Returns:
        pandas Dataframe: Dataframe with id as the index and a "Degree" column
    """
    data = pickle.load(file=open(fp, "rb"))

    return pd.DataFrame.from_dict(data, orient="index", columns=["Degree"])


def sample_quantiles(df, q, n):
    """ Sample n compounds from each of q quantiles within a id/degree dataframe

    Args:
        df (dataframe): SureChemBL id is the index, "Degree" holds the degree of that particular id
        q (int): number of quantiles to sample from (anything over 75 seems to break, not sure why)
        n (int): sample size from each quantile

    Returns:
        dataframe: sampled slice of original dataframe with max q compounds per quantile
    """
    df["quantile"] = pd.qcut(df.Degree,
                             q=50,
                             duplicates="drop",
                             labels=list(np.arange(0, 50, 1)))

    out = pd.concat([
        df[df["quantile"].eq(label)].sample(n, replace=True)
        for label in list(np.arange(0, 50, 1))
    ])
    #Remove duplicate indices (if possible)
    return out[~out.index.duplicated(keep="last")]


def main():
    ### Sample compounds across a range of degrees
    # (use quantiles to give a uniform sample across degrees)
    df = get_data("../../Downloads/id_degrees_2019-12.p")

    output_df = sample_quantiles(df, q=50, n=200)

    ### Link SurechemBL ids with inchi values
    SBL_df = pickle.load(
        file=open("../../Downloads/SureChemBL_allCpds.p", "rb"))
    print(type(SBL_df))

    sample = SBL_df[SBL_df["SureChEMBL_ID"].isin(list(output_df.index.values))]
    print(sample)
    sample.to_csv("sample_TESTFULLIDS_2019-12.csv")

    print(output_df)
    output_df.to_csv("sample_byDegreeQuantile_2019-12.csv")


if __name__ == "__main__":
    main()
