import rdkit.Chem as Chem
from rdkit import RDLogger
import pandas as pd
from tqdm import tqdm
import pickle
import os


def inchi_to_mol(inchi, index, database):
    """ Translates a smiles string to a mol file and saves it in the appropriate location

    Args:
        smiles (str): smiles description of a molecule
        index (int): identifying number of the molecule (unique to database)
        database (str): name of database (corresponds to a directory in 'Data')
    """
    try:
        mol = Chem.MolFromInchi(inchi)
        print(Chem.MolToMolBlock(mol),
              file=open(
                  "Data/AssemblyValues/" + database + "/" + str(index) + ".mol",
                  "w+"))
    except:
        pass


def smiles_to_mol(smiles, index, database):
    """ Translates a smiles string to a mol file and saves it in the appropriate location

    Args:
        smiles (str): smiles description of a molecule
        index (int): identifying number of the molecule (unique to database)
        database (str): name of database (corresponds to a directory in 'Data')
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        print(Chem.MolToMolBlock(mol),
              file=open(
                  "Data/AssemblyValues/" + database + "/" + str(index) + ".mol",
                  "w+"))
    except:
        pass


def move_completed_files(fp):
    """ Moves completed MA files from analysis directory to completed directory

    Args:
        fp (str): filepath to data directory containing Full & New Databases 
    """
    # print("----- Moving Full Database -----")
    # completed_files = [x for x in os.listdir(fp + "/FullDatabase/") if x.endswith(".txt")]

    # for f in tqdm(completed_files):
    #     os.rename(fp + "/FullDatabase/" + f, fp + "/FullDatabase_Done/" + f)
    #     os.rename(fp + "/FullDatabase/" + f[:-4] + ".mol", fp + "/FullDatabase_Done/" + f[:-4] + ".mol", )

    print("----- Moving New Database -----")
    completed_files = [
        x for x in os.listdir(fp + "/NewDatabase/") if x.endswith(".txt")
    ]

    for f in tqdm(completed_files):
        os.rename(fp + "/NewDatabase/" + f, fp + "/NewDatabase_Done/" + f)
        os.rename(
            fp + "/NewDatabase/" + f[:-4] + ".mol",
            fp + "/NewDatabase_Done/" + f[:-4] + ".mol",
        )


def main():
    #Progress meter
    tqdm.pandas()

    #Disable RDKit logger
    RDLogger.DisableLog("rdApp.*")

    # ## Full database to mol files
    # print("------ Full Database -----")
    # full_df = pd.read_csv("Data/AssemblyValues/fullCpds_AssemblyValues.csv")

    # full_df.progress_apply(lambda x: inchi_to_mol(
    #     x["inchi"], x["month"] + "_" + str(x.name), "FullDatabase"),
    #                        axis=1)

    # ## Full database to mol files
    # print("------ New Database -----")
    # new_df = pd.read_csv("Data/AssemblyValues/newCpds_AssemblyValues.csv")

    # new_df.progress_apply(lambda x: inchi_to_mol(x["inchi"], x["month"] + "_" +
    #                                              str(x.name), "NewDatabase"),
    #                       axis=1)

    # ## Cost database - random & percentiles - to mol
    # print("------ Cost Random  -----")
    # cost_df = pd.read_csv("Data/Cost/random_cost_50k.csv")
    # cost_df = cost_df.rename(columns={"Unnamed: 0": "ID"
    #                                  })  #Make sure ID column is named correctly

    # #Make sure ids are unique
    # # print(len(set(list(cost_df["ID"])))) #Should be 50k - YEP
    # cost_df.progress_apply(
    #     lambda x: smiles_to_mol(x["smiles"], x["ID"], "CostRandom"), axis=1)

    # print("------ Cost Percentiles  -----")
    # cost_percentile_df = pd.read_csv(
    #     "Data/Cost/random_cost_percentiles_50k.csv")
    # cost_percentile_df = cost_percentile_df.rename(
    #     columns={"Unnamed: 0": "ID"})  #Make sure ID column is named correctly)

    # cost_percentile_df.progress_apply(
    #     lambda x: smiles_to_mol(x["smiles"], x["ID"], "CostRandomPercentiles"),
    #     axis=1)

    # print("----- Reaxys Random -----")
    # reaxys_df = pd.read_csv("Data/Reaxys/random_100k.csv")
    # reaxys_df = reaxys_df.rename(
    #     columns={"Unnamed: 0": "ID"})  #Make sure ID column is named correctly

    # reaxys_df.progress_apply(
    #     lambda x: smiles_to_mol(x["smiles"], x["ID"], "ReaxysRandom"), axis=1)

    # print("----- Reaxys Product Percentiles -----")
    # reaxys_df = pd.read_csv(
    #     "Data/Reaxys/random_uniform_productPercentiles_110k.csv")
    # reaxys_df = reaxys_df.rename(
    #     columns={"Unnamed: 0": "ID"})  #Make sure ID column is named correctly

    # reaxys_df.progress_apply(lambda x: smiles_to_mol(
    #     x["smiles"], x["ID"], "ReaxysProductPercentiles"),
    #                          axis=1)

    print("----- Reaxys Reference Percentiles -----")
    reaxys_df = pd.read_csv(
        "Data/Reaxys/random_uniform_referencePercentiles_96k.csv")
    reaxys_df = reaxys_df.rename(
        columns={"Unnamed: 0": "ID"})  #Make sure ID column is named correctly

    reaxys_df.progress_apply(lambda x: smiles_to_mol(
        x["smiles"], x["ID"], "ReaxysReferencePercentiles"),
                             axis=1)

    # move_completed_files("Data/AssemblyValues")


if __name__ == "__main__":
    main()
