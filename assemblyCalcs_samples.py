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


def calculate_assembly(month, inchi):
    """ Calculate the assembly value of an inchi string

    Args:
        month (string): YYYY-MM description of the month where the compound was sampled from
        inchi (string): inchi representation of the SureChemBL compound

    Returns:
        dict: values of month, inchi, and the assembly index
    """
    ai = ac.calculate_ma(inchi,
                         60,
                         "monte-carlo",
                         num_frags_hist=5000,
                         path_samples=10000)
    return {"month": month, "inchi": inchi, "ai": ai}


def main():

    #Read in sampled compounds (updated for full 1000 compounds)
    #NOTE: add "_NEW" for new compounds found in each year (remove for all compounds)
    cpds = pickle.load(file=open("Data/sample_inchi_1000.p", "rb"))

    for year in np.arange(1987, 2020, 1):
        #Set up parallelization - a bit of overhead for setting it up, but that's fine
        pool = mp.Pool(64)

        #Build months in a specific year
        months = build_month_increments(year, year)

        date_cpd_sets = []
        for key, value in cpds.items():
            if key in months:
                for cpd in value:
                    date_cpd_sets.append((key, cpd))

        #Calculate assembly values for all inchis, save in a list holding dictionaries
        assemblies = [
            pool.apply(calculate_assembly, args=(s[0], s[1]))
            for s in date_cpd_sets
        ]

        pool.close()
        pool.join()

        #NOTE: include '_FULL_' when sampling all compounds
        pickle.dump(assemblies,
                    file=open("Data/assembly_values_1000_FULL_" + str(year) + ".p",
                              "wb"))


if __name__ == "__main__":
    main()
