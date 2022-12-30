import pickle
import pandas as pd
import rdkit.Chem as Chem
from tqdm import tqdm
from rdkit import RDLogger


def link_structure_IDs():
    """ Link randomly sampled structures to SureChemBL IDs

    Returns:
        None: but saves data to random_cost_50k_MA_IDs.csv on Macintosh 4 volume
    """
    random_cpds = pd.read_csv(
        "/Volumes/Macintosh HD 4/SureChemBL/Cost/random_cost_50k_MA.csv")
    print(random_cpds.head())

    surechembl_df = pickle.load(file=open(
        "/Volumes/Macintosh HD 4/SureChemBL/Cpd_Data/SureChemBL_allCpds.p",
        "rb"))
    print(surechembl_df.head())

    # Add surechembl ids to random cpd df through merging
    full_df = pd.merge(random_cpds,
                       surechembl_df,
                       how="left",
                       left_on=["smiles"],
                       right_on=["SMILES"])

    print(full_df)

    #Save merged df
    full_df.to_csv(
        "/Volumes/Macintosh HD 4/SureChemBL/Cost/random_cost_50k_MA_IDs.csv")


def build_month_increments(start, stop):
    """ Build month increments in the form YYYY-MM

    Args:
        start (int): Starting year
        stop (int): Ending year

    Returns:
        list: list of strings in the form YYYY-MM (e.g., "1980-01")
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


def find_earliest_date(ids):
    #Loop through every cpd date dictionary, filtering out by ids for each one

    months = build_month_increments(1962, 2020)
    #Reverse the month list to ensure the last date added is the earliest
    months.reverse()

    # #Testing
    id_date_dict = {}
    for month in tqdm(months):
        cpd_date_dict = pickle.load(file=open(
            "/Volumes/Macintosh HD 4/SureChemBL/CpdPatentIdsDates/Cpd_Date_Dict/cpd_date_dict_"
            + month + ".p", "rb"))

        for id in ids:
            if id in cpd_date_dict:
                id_date_dict[id] = cpd_date_dict[id]

        #print([cpd_date_dict[id] for id in ids if id in cpd_date_dict])

    #     print(results)

    return id_date_dict


def main():
    RDLogger.DisableLog('rdApp.*')
    tqdm.pandas()

    ### Cost/timewise analysis, starting with randomly sampled compounds

    # #Step 1: link structures with SureChemBL IDs
    # link_structure_IDs()

    # Step 2: link ids with starting dates (on only 2901 cpds - eventually, figure out why this is - probably something to do with smiles?)
    full_df = pd.read_csv(
        "/Volumes/Macintosh HD 4/SureChemBL/Cost/random_cost_50k_MA_IDs.csv")

    full_df = full_df.dropna()

    print(full_df)

    id_date_dict = find_earliest_date(list(full_df["SureChEMBL_ID"]))

    full_df["date"] = full_df["SureChEMBL_ID"].map(id_date_dict)

    print(full_df)

    full_df.to_csv(
        "/Volumes/Macintosh HD 4/SureChemBL/Cost/random_cost_50k_MA_IDs_dates.csv"
    )


if __name__ == "__main__":
    main()
