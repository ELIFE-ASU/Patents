import assemblycalculator as ac
import multiprocessing as mp
import pickle
import pandas as pd
import os


def calculate_assembly(inchi):
    """ Calculate the assembly value of an inchi string

    Args:
        month (string): YYYY-MM description of the month where the compound was sampled from
        inchi (string): inchi representation of the SureChemBL compound

    Returns:
        dict: values of month, inchi, and the assembly index
    """
    ai = ac.calculate_ma(inchi,
                         120,
                         "monte-carlo",
                         num_frags_hist=10000,
                         path_samples=20000)
    return {"inchi": inchi, "ai": ai}


def read_cpds(fp):
    """ Read inchi compounds from csv files

    Args:
        fp (string): Relative file path to csv file with SureChemBL cpd data

    Returns:
        list: list of all inchis contained in the csv file
    """
    data = pd.read_csv(fp)
    return data["InChI"].tolist()


def get_changing_percentileFiles(side, precision):
    """ Find the files which correspond to change percentiles

    Args:
        side (string): min/max, referring to the largest negative/positive changes, respectively
        precision (string): string representation of percentile precision

    Returns:
        list: list of file names which fit the appropriate criteria
    """
    files = []
    for f in os.listdir("Data/Cpd_Data/"):
        if f.startswith("ids_change_" + side + "Percentile_" +
                        precision) and f.endswith(".csv"):
            files.append(f)

    return files


def get_top_percentileFiles():
    """ Returns the csv files corresponding to compounds above the 99.99th percentile
    of total attachment values

    Returns:
        list: list of file names which fit the appropriate criteria
    """
    files = []
    for f in os.listdir("Data/Cpd_Data/"):
        if f.startswith("ids_above99_99percentile"):
            files.append(f)

    return files


def calculate_MAs(files):
    """ Wrapper function for MA calculation

    Args:
        files (list): list of all files which contain relevant data

    Returns:
        Writes a file containing inchis linked with assembly values
    """
    for f in files:
        cpds = read_cpds("Data/Cpd_Data/" + f)

        #Set up parallelization - a bit of overhead for setting it up, but that's fine
        pool = mp.Pool(64)

        #Calculate assembly values using MC method
        assemblies = pool.map(calculate_assembly, cpds)

        pool.close()

        pickle.dump(assemblies,
                    file=open("Data/Cpd_Data/" + f[:-4] + "_assembly.p", "wb"))


def main():
    """ Steps
    1. Read in compounds from a specific file
    2. MC assembly algorithm
    3. Save calculations using name  + "assembly" - link inchi & MA
    4. Link assembly values to csv file (eventually, probably do this in a separate script)
    """

    ### SMALLEST & LARGEST CHANGE VALUES ###
    #Options: min/max, 0.1/0.01
    for option in [("min", "0.1"), ("min", "0.01"), ("max", "0.1"),
                   ("max", "0.01")]:
        files = get_changing_percentileFiles(option[0], option[1])
        calculate_MAs(files)

    ### TOP ATTACHMENT VALUES ###
    files = get_top_percentileFiles()
    calculate_MAs(files)


if __name__ == "__main__":
    main()
