import pickle
from tqdm import tqdm


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


def read_ids(fp):
    """ Reads all SureChemBL ids from a certain month.

    Args: 
        fp, filepath to pickled ids

    Returns: 
        list of lists, each containing a cluster of SureChemBL ids ([0] being the lcc)
    """
    return pickle.load(file=open(fp, "rb"))


def getIds(ids, allIds, month):
    """ Finds all new SureChemBL IDs added in a given month, as well as 
    which of those compounds were added to the LCC (and those which weren't)

    Args:
        ids: monthly id 2D list (each sublist is a cluster, [0] is LCC)
        allIds: list of all SureChemBL ids seen thus far
        month: month of analysis

    Returns:
        newIds (list): all new ids added
    """
    #Get new ids
    newIds = list(set(ids) - set(allIds))
    pickle.dump(
        newIds,
        file=open(
            "/Volumes/Macintosh HD 4/SureChemBL/CpdPatentIdsDates/New_Ids/newIds_updated_"
            + month + ".p", "wb"))

    # #All new LCC IDs that month
    # newlccIds = list(set(ids[0]) - set(lccIds))
    # pickle.dump(newlccIds,
    #             file=open(
    #                 "/scratch/jmalloy3/Patents/CpdPatentIdsDates/newlccIds_" +
    #                 month + ".p", "wb"))

    # #Existing compounds which were added to LCC
    # existing_newlccIds = list(set(newlccIds) - set(newIds))

    # #Finds the newIds which ARE NOT in the LCC
    # nonlcc_newIds = list(set(newIds) - set(ids[0]))

    return newIds


def calculate_LCC_stats(allIds, lccIds, newIds, newlccIds, existing_newlccIds,
                        nonlcc_newIds, month):
    """ Left over from largest connected component analysis

    Args:
        allIds (list): all IDs found up to a given month
        lccIds (list): all IDs in the LCC
        newIds (list): novel IDs
        newlccIds (list): novel IDs, only those in the largest connected component
        existing_newlccIds (list): new-to-the-LCC ids
        nonlcc_newIds (list): new IDs which are not part of the LCC
        month (str): month to analyze

    Returns:
        list: descriptive stats of the LCC
    """

    return [
        month,  #Month of analysis
        len(newIds),  #Number of new ids added that month
        len(newlccIds),  #Number of NEW ids added to the LCC that month
        len(nonlcc_newIds),  #Number of NEW ids that were not added to the LCC 
        len(existing_newlccIds),
        1 -
        (float(len(nonlcc_newIds)) /
         len(newIds)),  #Percentage of new ids which immediately entered the LCC
        float(len(existing_newlccIds)) /
        (len(allIds) - len(lccIds)
        ),  #Percentage of existing ids (outside LCC) which were added to LCC
        len(allIds),  #Total number of compounds
        len(lccIds)  #Total number of compounds in the LCC
    ]


def main():
    #Set up first month
    fp = "/Volumes/Macintosh HD 4/SureChemBL/CpdPatentIdsDates/Unique_Cpds/"
    ids = read_ids(fp + "unique_cpds_1962-01.p")
    allIds = list(set(ids))
    # lccIds = ids[0] ## Leftover from largest connected component analysis

    for month in tqdm(build_month_increments(1963, 2022)):
        ### FIND NEW COMPOUNDS ###
        ids = read_ids(fp + "unique_cpds_" + month + ".p")

        newIds = getIds(ids, allIds, month)

        allIds.extend(newIds)

    # ### LCC OUTPUT ###
    # df = pd.DataFrame(stats,
    #                   columns=[
    #                       "month", "newIds", "newIds_newLCC", "newIds_nonLCC",
    #                       "oldIds_newLCC", "newIdsinLCC_percentage",
    #                       "oldIdsNewinLCC_percentage", "totalIds", "totalLCCIds"
    #                   ])
    # df.to_csv("/scratch/jmalloy3/Patents/NetworkStats/lcc_stats_TEST.csv")


if __name__ == "__main__":
    main()
