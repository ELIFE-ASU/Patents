# import igraph as ig
# import numpy as np
import pickle
import pandas as pd
from tqdm import tqdm
import os
import heapq
import scipy.stats as stats
from random import sample


def build_cpd_df(fp):
    """ Takes 29 separate compound data files and combines them into a single pandas dataframe for ease of access

    Args:
        fp (string): Filepath to SureChemBL data files (assuming G drive goes to jmalloy3 Google Account)

    Returns:
        None - but does write a pickled dataframe to SureChemBL_Patents/Cpd_Data/ directory
    """
    dfs = []
    for f in tqdm(os.listdir(fp)):
        if f.endswith(".txt"):
            dfs.append(pd.read_csv(fp + f, sep="\t", header=0))

    df = pd.concat(dfs, ignore_index=True)
    print(df)
    pickle.dump(df, file=open(fp + "SureChemBL_allCpds.p", "wb"))

    del df


def find_highest_degrees(df, n, start, stop):
    """ Finds the n highest-degree compounds within a specific date range

    Saves various data associated with those n comopunds - smiles, inchi,
    inchikey, degree, preferential attachment value

    Args:
        df (pandas dataframe): dataframe containing all SureChemBL compounds
        n (int): the number of highest-degree compounds to select
        start (int): 1st year of the range
        stop (int): last year of the range
    """
    print("----------", start, stop, "----------")

    #Finding the top 10 preferential attachment compounds (from 1980-1984 as a test)
    full_id_degrees = pickle.load(file=open(
        "G:\\Shared drives\\SureChemBL_Patents\\Degrees\\full_id_degrees_" +
        str(start) + "_" + str(stop) + ".p", "rb"))
    pref_attach_dict = pickle.load(file=open(
        "G:\\Shared drives\\SureChemBL_Patents\\pref_attach_dict_" +
        str(start) + "_" + str(stop) + ".p", "rb"))

    #Find n compounds with largest degree
    highest_degree_cpds = heapq.nlargest(n,
                                         full_id_degrees,
                                         key=full_id_degrees.get)

    highest_degree_cpds_df = df[df["SureChEMBL_ID"].isin(highest_degree_cpds)]

    pref_attach_values = list(pref_attach_dict.values())

    #Extra information to be added to the csv output file
    degrees = []
    pref_attach_highestCpd_values = []
    pref_attach_percentiles = []

    for cpd in tqdm(highest_degree_cpds_df["SureChEMBL_ID"]):
        #Degree of compound
        degrees.append(full_id_degrees[cpd][-1])

        #Preferential attachment value
        pref_attach_highestCpd_values.append(pref_attach_dict[cpd])

        #Percentile of preferential attachment value
        pref_attach_percentiles.append(
            stats.percentileofscore(pref_attach_values, pref_attach_dict[cpd]))

    highest_degree_cpds_df["degree"] = degrees
    highest_degree_cpds_df["pref_attach_value"] = pref_attach_highestCpd_values
    highest_degree_cpds_df["pref_attach_percentile"] = pref_attach_percentiles

    highest_degree_cpds_df.to_csv(
        "G:\\Shared drives\\SureChemBL_Patents\\Cpd_Data/highest_degree_data_" +
        str(start) + "_" + str(stop) + "_1000.csv")

    print()


def find_llanos_cpds(fp, df):
    """ Tests various compounds found in Llanos et al (2019) in SureChemBL data

    Llanos et al used Reaxys data to find the most popular compounds. This checks
    where those compounds appear, if at all, in SureChembL patent data

    Args:
        df (pandas dataframe): dataframe of all SureChemBL chemistry
    """

    cpds_1980_2015_inchi = {
        "acetic anhydride":
            "InChI=1S/C4H6O3/c1-3(5)7-4(2)6/h1-2H3",
        "methanol":
            "InChI=1S/CH4O/c1-2/h2H,1H3",
        "methyl iodide":
            "InChI=1S/CH3I/c1-2/h1H3",
        "diazomethane":
            "InChI=1S/CH2N2/c1-3-2/h1H2",
        "formaldehyde":
            "InChI=1S/CH2O/c1-2/h1H2",
        "benzaldehyde":
            "InChI=1S/C7H6O/c8-6-7-4-2-1-3-5-7/h1-6H",
        "copper(II) oxide":
            "InChI=1S/Cu.O",
        "ethanol":
            "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "benzoyl chloride":
            "InChI=1S/C7H5ClO/c8-7(9)6-4-2-1-3-5-6/h1-5H",
        "carbon monoxide":
            "InChI=1S/CO/c1-2",
        "water (2000)":
            "InChI=1S/H2O/h1H2",
        "Trifluoroacetic acid (2000)":
            "InChI=1S/C2HF3O2/c3-2(4,5)1(6)7/h(H,6,7)",
        "Phenylacetylene (2000)":
            "InChI=1S/C8H6/c1-2-8-6-4-3-5-7-8/h1,3-7H",
        "benzyl bromide (2000)":
            "InChI=1S/C7H7Br/c8-6-7-4-2-1-3-5-7/h1-5H,6H2"
    }

    #Find stats for Llanos compounds - use 2015 data for stats (I really need to make a consensus graph)
    full_id_degrees = pickle.load(file=open(
        "G:\\Shared drives\\SureChemBL_Patents\\Degrees\\full_id_degrees_2015_2019.p",
        "rb"))
    pref_attach_dict = pickle.load(file=open(
        "G:\\Shared drives\\SureChemBL_Patents\\pref_attach_dict_2015_2019.p",
        "rb"))
    pref_attach_values = list(pref_attach_dict.values())

    #Loop through Llanos compounds
    with open(fp + "llanos_cpds.csv", "a") as f:
        f.write(
            "name,inchi,SureChemBL_ID,degree,pref_attach_value,pref_attach_percentile\n"
        )
        for name, inchi in cpds_1980_2015_inchi.items():
            s = df[df["InChI"] == inchi]
            if not s.empty:  #if SureChemBL holds that compound, save id & stats
                #Degree of compound
                degree = full_id_degrees[s.iloc[0]["SureChEMBL_ID"]][-1]

                #Preferential attachment value
                pref_attach_value = pref_attach_dict[s.iloc[0]["SureChEMBL_ID"]]

                #Percentile of preferential attachment value
                pref_attach_percentile = stats.percentileofscore(
                    pref_attach_values,
                    pref_attach_dict[s.iloc[0]["SureChEMBL_ID"]])

                f.write(name + ",\"" + inchi + "\"," +
                        s.iloc[0]["SureChEMBL_ID"] + "," + str(degree) + "," +
                        str(pref_attach_value) + "," +
                        str(pref_attach_percentile) + "\n")

            else:  #if not, no name nor stats
                f.write(name + ",\"" + inchi + "\",na,na,na,na\n")


def build_month_increments(start, stop):
    """ Build all monthly increments from the start year to stop year in the
    format YEAR-MONTH

    Args:
        start (int): start year of increments
        stop (int): end year of increments

    Returns:
        list: list of strings holding the YEAR-MONTH increments
    """
    months = []
    while start <= stop:
        for month in [
                "01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
                "11", "12"
        ]:
            months.append(str(start) + "-" + month)
        start += 1

    return months


def sample_compounds_unique(n, months, cpd_df):
    """ Sample compounds which are uniquely added in a specific month

    This uniquess is determined by determing when a compound is added in a month
    and has not been present in the patent record before that month.

    Args:
        n (int): Number of compounds to sample every month
        months (list): list of months to sample from
        cpds (list): all SureChemBL IDs of compounds added in a specific month
        cpd_df (pandas dataframe): Master dataframe of all compounds
    """
    sample_inchis = {}

    print("----- Sampling Novel Compounds -----")
    fp = "/Volumes/Macintosh HD 4/SureChemBL/CpdPatentIdsDates/New_Ids/"
    for month in tqdm(months):
        new_cpds = pickle.load(file=open(fp + "newIds_updated_" + month +
                                         ".p", "rb"))
        #Only sample if there are more than 1000 compounds
        if len(new_cpds) > n:
            sample_cpds = sample(new_cpds, n)
        else:
            sample_cpds = new_cpds

        sub_df = cpd_df[cpd_df["SureChEMBL_ID"].isin(sample_cpds)]
        sample_inchis[month] = list(sub_df["InChI"])

    print("\n----- Saving compounds -----")
    pickle.dump(sample_inchis,
                file=open("Data/sample_inchi_1000_NEW_1976-1979.p", "wb"))


def sample_compounds(n, months, cpd_df):
    """ Sample n compounds from each month, initially with overlap allowed

    Args:
        n (int): number of compounds to sample
        months (string): description of month, e.g. 1980-01
        cpd_df (pandas dataframe): contains information for each compound in SureChemBL, including InChIKey

    Returns:
        list: list of all randomly sampled compounds (in inchi?)
    """

    #Inchis for all sampled compounds
    sample_inchis = {}
    sample_ids = {}

    fp = "/Volumes/Macintosh HD 4/SureChemBL/CpdPatentIdsDates/Unique_Cpds/"
    print("----- Sampling Full Compounds ------\n")
    for month in tqdm(months):
        cpds = pickle.load(file=open(fp + "unique_cpds_" + month + ".p", "rb"))

        sample_cpds = sample(cpds, n)

        sub_df = cpd_df[cpd_df["SureChEMBL_ID"].isin(sample_cpds)]
        sample_inchis[month] = list(sub_df["InChI"])
        sample_ids[month] = list(sub_df["SureChEMBL_ID"])

    #Save sampled inchis to pickle files
    print("\n----- Saving Data -----")
    pickle.dump(sample_inchis,
                file=open("Data/sample_inchi_" + str(n) + "_1976-1979.p", "wb"))
    

def count_cpds(months):
    all_cpds = []
    for month in months:
        unique_cpds = pickle.load(open(f"/../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Unique_Cpds/unique_cpds_{month}.p", "rb"))
        all_cpds = list(set(all_cpds + unique_cpds))
        
    print(f"Total cpds, 1962-1979: f{len(all_cpds)}")

def main():
    # # ### Highest Degree compounds ###
    # data_fp = "/Volumes/Macintosh HD 4/SureChemBL/Cpd_Data/"
    # # build_cpd_df(data_fp)  #NOTE: only needs to be run once

    # cpd_df = pickle.load(file=open(data_fp + "SureChemBL_allCpds.p", "rb"))
    # print(cpd_df.columns)

    # ### Statistics over highest degree compounds ###
    # n = 1000  #Number of compounds to find
    # for range in [(1980, 1984), (1985, 1989), (1990, 1994), (1995, 1999),
    #               (2000, 2004), (2005, 2009), (2010, 2014), (2015, 2019)]:
    #     find_highest_degrees(cpd_df, n, range[0], range[1])

    # ### Testing Llanos et al (2019) compounds ###
    # find_llanos_cpds(data_fp, cpd_df)

    # ### Sampling compounds for MA analysis ###
    # sample_compounds_unique(1000, build_month_increments(1976, 1979), cpd_df)
    # sample_compounds(1000, build_month_increments(1976, 1979), cpd_df)

    ### Counting Stats (for SI) ###

    ## Get all cpds before 1980
    count_cpds(build_month_increments(1962,1979))


if __name__ == "__main__":
    main()
