import pickle
from tqdm import tqdm

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

def count_cpds(months):
    all_patents = []
    for month in tqdm(months):
        unique_patents = pickle.load(open(f"/../../../mnt/Archive/Shared/PatentData/SureChemBL/CpdPatentIdsDates/Unique_Patents/unique_patents_{month}.p", "rb"))
        all_patents = list(set(all_patents + unique_patents))
        
    print(f"Total patents, 1962-1979: {len(all_patents)}")

def main():
    count_cpds(build_month_increments(1962,1979))

if __name__ == "__main__":
    main()