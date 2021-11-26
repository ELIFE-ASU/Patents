import assemblycalculator as ac
import multiprocessing as mp
import pickle
import numpy as np


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
    return {"month": month, "inchi": inchi, "ai": ai}


def main():
    """ Steps
    1. Read in compounds from a specific file
    2. MC assembly algorithm
    3. Save calculations using name  + "assembly" - link inchi & MA?
    4. Link assembly values to csv file (eventually, probably do this in a separate script)
    """

    #Set up parallelization - a bit of overhead for setting it up, but that's fine
    pool = mp.Pool(64)

    #Calculate assembly values using MC method
    assemblies = pool.map(calculate_assembly, data)

    pool.join()
    pool.close()

    pickle.dump(assemblies,
                file=open("results_MC.p",
                            "wb"))



if __name__ == "__main__":
    main()
