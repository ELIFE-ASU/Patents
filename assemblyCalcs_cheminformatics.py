import rdkit.Chem as Chem
import rdkit.Chem.rdMolDescriptors as Desc
from rdkit import RDLogger
import pickle
from tqdm import tqdm
import pandas as pd
import numpy
import os


def get_mw(inchi):
    """ Finds the molecular weight for a given inchi string

    Args:
        inchi (str): inchi string of a SureChemBL compound

    Returns:
        float: molecular weight (or NaN if rdkit failure)
    """
    try:
        m = Chem.MolFromInchi(inchi)
        return Desc.CalcExactMolWt(m)
    except:
        return numpy.nan


def main():
    RDLogger.DisableLog('rdApp.*')

    #Test: load an assembly subset, find MWs for all cpds
    new_df = pd.DataFrame()
    full_df = pd.DataFrame()

    for f in tqdm(os.listdir("Data/AssemblyValues/")):
        if f.startswith("assembly_values_1000_"):
            data = pickle.load(file=open("Data/AssemblyValues/" + f, "rb"))
            df = pd.DataFrame(data)

            df["mw"] = df["inchi"].apply(get_mw)

            if "FULL" in f:
                full_df = full_df.append(df)
            else:
                new_df = new_df.append(df)

    new_df.to_csv("Data/AssemblyValues/newCpds_AssemblyValues.csv")
    full_df.to_csv("Data/AssemblyValues/fullCpds_AssemblyValues.csv")


if __name__ == "__main__":
    main()
