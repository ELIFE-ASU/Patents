import numpy as np
import pickle
import pandas as pd


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


def getIds(ids, allIds, lccIds, month):
    """ Finds all new SureChemBL IDs added in a given month, as well as 
    which of those compounds were added to the LCC (and those which weren't)

    Args:
        ids: monthly id 2D list (each sublist is a cluster, [0] is LCC)
        allIds: list of all SureChemBL ids seen thus far
        lccIds: all ids present in LCC
        month: month of analysis

    Returns:
        newIds (list): all new ids added
        newlccIds (list): all new ids which are present in LCC
        existing_newlccIds (list): all ids which were previously found, but were added to LCC this month
        nonlcc_newIds (list): all new ids that were NOT added to the LCC
    """
    #Get new ids
    newIds = list(set([id for sublist in ids for id in sublist]) - set(allIds))
    pickle.dump(newIds,
                file=open(
                    "/scratch/jmalloy3/Patents/CpdPatentIdsDates/newIds_" +
                    month + ".p", "wb"))

    #All new LCC IDs that month
    newlccIds = list(set(ids[0]) - set(lccIds))
    pickle.dump(newlccIds,
                file=open(
                    "/scratch/jmalloy3/Patents/CpdPatentIdsDates/newlccIds_" +
                    month + ".p", "wb"))

    #Existing compounds which were added to LCC
    existing_newlccIds = list(set(newlccIds) - set(newIds))

    #Finds the newIds which ARE NOT in the LCC
    nonlcc_newIds = list(set(newIds) - set(ids[0]))

    return newIds, newlccIds, existing_newlccIds, nonlcc_newIds


def calculate_LCC_stats(allIds, lccIds, newIds, newlccIds, existing_newlccIds,
                    nonlcc_newIds, month):

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
    ids = read_ids("/scratch/jmalloy3/Patents/NetworkStats/lcc_ids_1980-01.p")
    allIds = list(set([id for sublist in ids for id in sublist]))
    lccIds = ids[0]

    stats = []

    for month in build_month_increments(1980,2019):
        ### LCC STATS ###
        ids = read_ids("/scratch/jmalloy3/Patents/NetworkStats/lcc_ids_" +
                       month + ".p")
        # ids_2 = read_ids("Data/NetworkStats/lcc_ids_1980-03.p")

        newIds, newlccIds, existing_newlccIds, nonlcc_newIds = getIds(
            ids, allIds, lccIds, month)
        stats.append(
            calculate_LCC_stats(allIds, lccIds, newIds, newlccIds,
                            existing_newlccIds, nonlcc_newIds, month))

        allIds.extend(newIds)
        lccIds.extend(newlccIds)        

    ### LCC OUTPUT ###
    df = pd.DataFrame(stats,
                      columns=[
                          "month", "newIds", "newIds_newLCC", "newIds_nonLCC",
                          "oldIds_newLCC", "newIdsinLCC_percentage",
                          "oldIdsNewinLCC_percentage", "totalIds", "totalLCCIds"
                      ])
    df.to_csv("/scratch/jmalloy3/Patents/NetworkStats/lcc_stats_TEST.csv")


if __name__ == "__main__":
    main()
