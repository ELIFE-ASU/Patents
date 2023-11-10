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

def check_iso(candidate_fragments, all_frags, month):
    """ Finds the unique fragments within a set of candidates

    Args:
        candidate_fragments (list): list of igraph graphs from a single text file
        all_frags (list): list of all unique igraph graphs from previous files
        month (str): month that is currently being analyzed

    Returns:
        list: updated list of unique igraph graphs
    """
    #Check all candidate fragments...
    novel_frags = []
    novel_frag_count = 0
    for possible_f in candidate_fragments:
        is_unique = True
        #...against all unique fragments
        for i in range(len(all_frags)):
            if possible_f.isomorphic_vf2(all_frags[i],
                                         color1=possible_f.vs["color"],
                                         color2=all_frags[i].vs["color"],
                                         edge_color1=possible_f.es["color"],
                                         edge_color2=all_frags[i].es["color"]):
                #if candidate is isomorphic, it is not unique, therefore break off check
                is_unique = False

            if not is_unique:
                break

        #If it makes it through the full list of fragments, it is unique!
        if is_unique:
            all_frags.append(possible_f)
            novel_frags.append(possible_f)
            novel_frag_count += 1

    #Write novel fragments to a pickle file
    with open(f"~/../../../mnt/Archive/Shared/PatentData/SureChemBL/AssemblyValues/Fragments/NewFrags/newFrags_{month}_updated.p", "wb") as f:
        pickle.dump(novel_frags, f)

    #Return updated full list
    return all_frags, novel_frag_count


def main():
    ### Find all novel fragments within each month (1976-2022, inclusive), meant to run on Cradle

    months = build_month_increments(1976, 2022)
    all_frags = []
    novel_frag_counts = []
    
    for month in tqdm(months):
        #Load full fragments for each month
        full_frags = pickle.load(f"~/../../../mnt/Archive/Shared/PatentData/SureChemBL/AssemblyValues/Fragments/fullFrags_{month}.p", "rb")

        #If Jan 1976, all fragments is full fragment list
        if month=="1976-01":
            all_frags = full_frags
            novel_frag_counts.append(len(full_frags))
            with open(f"~/../../../mnt/Archive/Shared/PatentData/SureChemBL/AssemblyValues/Fragments/NewFrags/newFrags_{month}_updated.p", "wb") as f:
                pickle.dump(full_frags, f)
            break

        all_frags, monthly_novel_frag_count = check_iso(full_frags, all_frags, month)
        novel_frag_counts.append(monthly_novel_frag_count)

    print(novel_frag_counts)
    with open("~/../../../mnt/Archive/Shared/PatentData/SureChemBL/AssemblyValues/Fragments/NewFrags/newFrags_count_updated.p", "wb") as f:
        pickle.dump(novel_frag_counts, f)

        


if __name__ == "__main__":
    main()
