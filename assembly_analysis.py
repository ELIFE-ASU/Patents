import pickle
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
from datetime import datetime as dt
from tqdm import tqdm


def load_files(fp):
    """ Loads a given file

    Args:
        fp (str): filepath to a particular datafile

    Returns:
        data (usually pandas dataframe): data file
    """
    if fp.endswith(".csv"):
        return pd.read_csv(fp)
    elif fp.endswith(".p"):
        return pickle.load(file=open(fp, "rb"))
    else:
        return None


def get_dates_inchis():
    """ Find all sampled months & inchis across the patent project

    Returns:
        lists: list of months & inchis, all in order of sampling
    """
    print("--- Loading patent cpds ---")
    patents_df = load_files(
        "Data/AssemblyValues/ALLSAMPLEDcpds_AssemblyGo_MW.csv")
    patents_months = list(patents_df["month"])
    patents_inchis = list(patents_df["InChI"])

    print("--- Loading full/new cpds --- ")
    full_df = load_files("Data/AssemblyValues/fullCpds_AssemblyValues.csv")
    full_months = list(full_df["month"])
    full_inchis = list(full_df["inchi"])

    new_df = load_files("Data/AssemblyValues/newCpds_AssemblyValues.csv")
    new_months = list(new_df["month"])
    new_inchis = list(new_df["inchi"])

    print("--- Loading 1970s & 2020s cpds --- ")
    new_70s = load_files("Data/sample_inchi_1000_NEW_1976-1979.p")
    for key, inchis in new_70s.items():
        new_months.extend([key] * len(inchis))
        new_inchis.extend(inchis)

    new_20s = load_files("Data/sample_inchi_1000_NEW_2020-2022.p")
    for key, inchis in new_20s.items():
        new_months.extend([key] * len(inchis))
        new_inchis.extend(inchis)

    full_70s = load_files("Data/sample_inchi_1000_1976-1979.p")
    for key, inchis in full_70s.items():
        full_months.extend([key] * len(inchis))
        full_inchis.extend(inchis)

    full_20s = load_files("Data/sample_inchi_1000_2020-2022.p")
    for key, inchis in full_20s.items():
        full_months.extend([key] * len(inchis))
        full_inchis.extend(inchis)

    months = patents_months + full_months + new_months
    inchis = patents_inchis + full_inchis + new_inchis

    return months, inchis


def get_mw(inchis):
    """ Finds the exact molecular weight of a list of inchis

    Args:
        inchis (list): inchi descriptors of a set of molecules

    Returns:
        list: list of molecular weights
    """
    mws = []
    fail_count = 0
    for inchi in tqdm(inchis):
        try:
            mw = Descriptors.ExactMolWt(Chem.MolFromInchi(inchi))
        except:
            mw = -1
            fail_count += 1

        mws.append(mw)

    print("Failures:", fail_count, end="\n\n")

    return mws


def get_bonds(inchis):
    """ Get the number of bonds from a set of molecules

    Args:
        inchis (list): inchi descriptions of molecules

    Returns:
        list: list of number of bonds
    """

    bonds = []
    fail_count = 0
    for inchi in tqdm(inchis):
        try:
            bond = len(Chem.MolFromInchi(inchi).GetBonds())
        except:
            bond = -1
            fail_count += 1

        bonds.append(bond)

    print("Failures:", fail_count, end="\n\n")

    return bonds


def main():
    ### End goal - list of months with molecular weights

    months, inchis = get_dates_inchis()

    # mws = get_mw(inchis)
    bonds = get_bonds(inchis)

    # ## Save all MW data into a dataframe

    # mw_df = pd.DataFrame({"Inchi": inchis, "month": months, "mw": mws})
    # print(mw_df.head())
    # mw_df.to_csv("Data/AssemblyValues/ALLCPDS_mw.csv")

    ## Save all bonds data into a dataframe

    bond_df = pd.DataFrame({"Inchi": inchis, "month": months, "bonds": bonds})
    print(bond_df.head())
    bond_df.to_csv("Data/AssemblyValues/ALLCPDS_bonds.csv")


if __name__ == "__main__":
    main()
