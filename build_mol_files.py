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


def move_completed_files(fp):
    """ Moves completed MA files from analysis directory to completed directory

    Args:
        fp (str): filepath to data directory containing Full & New Databases 
    """
    print("----- Moving Full Database -----")
    completed_files = [
        x for x in os.listdir(fp + "/FullDatabase/") if x.endswith(".txt")
    ]

    for f in tqdm(completed_files):
        os.rename(fp + "/FullDatabase/" + f, fp + "/FullDatabase_Done/" + f)
        os.rename(
            fp + "/FullDatabase/" + f[:-4] + ".mol",
            fp + "/FullDatabase_Done/" + f[:-4] + ".mol",
        )

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

    print("----- Moving ReaxysRandom -----")
    completed_files = [
        x for x in os.listdir(fp + "/ReaxysRandom/") if x.endswith(".txt")
    ]
    for f in tqdm(completed_files):
        os.rename(fp + "/ReaxysRandom/" + f, fp + "/ReaxysRandom_Done/" + f)
        os.rename(
            fp + "/ReaxysRandom/" + f[:-4] + ".mol",
            fp + "/ReaxysRandom_Done/" + f[:-4] + ".mol",
        )

    print("----- Moving ReaxysProductPercentiles -----")
    completed_files = [
        x for x in os.listdir(fp + "/ReaxysProductPercentiles/")
        if x.endswith(".txt")
    ]
    for f in tqdm(completed_files):
        os.rename(fp + "/ReaxysProductPercentiles/" + f,
                  fp + "/ReaxysProductPercentiles_Done/" + f)
        os.rename(
            fp + "/ReaxysProductPercentiles/" + f[:-4] + ".mol",
            fp + "/ReaxysProductPercentiles_Done/" + f[:-4] + ".mol",
        )

    print("----- Moving ReaxysReferencePercentiles -----")
    completed_files = [
        x for x in os.listdir(fp + "/ReaxysReferencePercentiles/")
        if x.endswith(".txt")
    ]
    for f in tqdm(completed_files):
        os.rename(fp + "/ReaxysReferencePercentiles/" + f,
                  fp + "/ReaxysReferencePercentiles_Done/" + f)
        os.rename(
            fp + "/ReaxysReferencePercentiles/" + f[:-4] + ".mol",
            fp + "/ReaxysReferencePercentiles_Done/" + f[:-4] + ".mol",
        )

    print("----- Moving Author Cpds -----")
    completed_files = [
        x for x in os.listdir(fp + "/AuthorCpds/")
        if x.endswith(".txt")
    ]
    for f in tqdm(completed_files):
        try:
            os.rename(fp + "/AuthorCpds/" + f,
                    fp + "/AuthorCpds_Done/" + f)
        except OSError as e :
            print(e)
        try:
            os.rename(
            fp + "/AuthorCpds/" + f[:-4] + ".mol",
            fp + "/AuthorCpds_Done/" + f[:-4] + ".mol",
            )
        except OSError as e:
            print(e)

    print("----- Moving Assignee Cpds -----")
    completed_files = [
        x for x in os.listdir(fp + "/AssigneeCpds/")
        if x.endswith(".txt")
    ]
    for f in tqdm(completed_files):
        os.rename(fp + "/AssigneeCpds/" + f,
                  fp + "/AssigneeCpds_Done/" + f)
        os.rename(
            fp + "/AssigneeCpds/" + f[:-4] + ".mol",
            fp + "/AssigneeCpds_Done/" + f[:-4] + ".mol",
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

    # ## New database to mol files
    # print("------ New Database -----")
    # new_df = pd.read_csv("Data/AssemblyValues/newCpds_AssemblyValues.csv")

    # new_df.progress_apply(lambda x: inchi_to_mol(x["inchi"], x["month"] + "_" +
    #                                              str(x.name), "NewDatabase"),
    #                       axis=1)

    # ## Author compounds to mol files
    # print("------ Author Compounds -----")
    # cpd_df = pd.read_csv("Data/Patents/author_cpds_structures.csv")

    # cpd_df.progress_apply(
    #     lambda x: inchi_to_mol(x["InChI"], x["SureChEMBL_ID"], "AuthorCpds"),
    #     axis=1)

    # ## Author compounds to mol files
    # print("------ Assignee Compounds -----")
    # cpd_df = pd.read_csv("Data/Patents/assignee_cpds_structures.csv")

    # cpd_df.progress_apply(
    #     lambda x: inchi_to_mol(x["InChI"], x["SureChEMBL_ID"], "AssigneeCpds"),
    #     axis=1)

    ## Move completed MA files on Agave
    move_completed_files("Data/AssemblyValues")


if __name__ == "__main__":
    main()
