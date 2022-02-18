import assemblycalculator as ac
import multiprocessing as mp
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


def calculate_statistics(inchi):
    """ Uses the assembly calculator package to find the assembly index of a compound

    Uses try/except blocks to account for various assembly calculator/rdkit errors,
    giving 240 seconds as a timeout (only ~9500 compounds)

    Args:
        inchi (string): compound inchi descriptor

    Returns:
        set: dictionary of MA values (mc and fragment methods) and inchi
    """
    try:
        #Use these (small) parameters for assembly index sampling -
        #  goal is to determine ranges, not necessarily exact values
        mc_ai = ac.calculate_ma(inchi,
                                240,
                                "monte-carlo",
                                num_frags_hist=5000,
                                path_samples=10000)
    except:
        mc_ai = -1

    try:
        frag_ai = ac.calculate_ma(inchi, 240, "fragment")
    except:
        frag_ai = -1

    return {"mc_ai": mc_ai, "frag_ai": frag_ai, "inchi": inchi}


def get_stats(inchis):
    """ Wrapper for parallel assembly value computations

    Args:
        inchis (list): list of all inchis to calculate MA values

    Returns:
        list: list of dictionaries, all continaing MC/Fragment MA values & corresponding inchi
    """
    #Set up parallelization
    pool = mp.Pool(64)

    assembly_results = [pool.map(calculate_statistics, inchis)]

    pool.close()
    pool.join()

    return assembly_results


def main():
    ### Sample compounds across a range of degrees
    # (use quantiles to give a uniform sample across degrees)
    # df = get_data("../../Downloads/id_degrees_2019-12.p")

    # output_df = sample_quantiles(df, q=50, n=200)

    ### Link SurechemBL ids with inchi values
    # SBL_df = pickle.load(
    #     file=open("../../Downloads/SureChemBL_allCpds.p", "rb"))
    # print(type(SBL_df))

    # sample = SBL_df[SBL_df["SureChEMBL_ID"].isin(list(output_df.index.values))]
    # print(sample)
    # sample.to_csv("sample_TESTFULLIDS_2019-12.csv")

    # print(output_df)
    # output_df.to_csv("sample_byDegreeQuantile_2019-12.csv")

    ### Assembly stats over degree quantiles
    df_inchis = pd.read_csv(
        "Data/AssemblyDegreeCorr/sample_TESTFULLIDS_2019-12.csv")

    #Calculate and save assembly results (will worry about correlation later)
    assembly_results = get_stats(list(df_inchis["InChI"]))
    pickle.dump(assembly_results,
                file=open(
                    "Data/AssemblyDegreeCorr/TESTFULLIDS_assemblyValues.p",
                    "wb"))


if __name__ == "__main__":
    main()
