import pandas as pd
import os
from tqdm import tqdm
import pickle

""" Reads in SureChemBL mapping data
    Input: filepath
    Output: dataframe of SureChemBL ID, no. of occurances (at time of data dump), date of patent (YYYY-MM-DD format)
"""
def get_data(fp):
    df = pd.read_csv(fp, delimiter='\t', usecols=[0,3,5], names=["ID", "Occurances", "Date"])
    return dict(zip(df.ID, df.Occurances))

""" Find the occurances for each data update (including the first)
    Input: None (standalone function - assumes data is in Data/SureChemblMAP directory)
    Output: pickle files of each update, including a final occurance dictionary
"""
def occurances():
    d = get_data("Data/SureChemblMAP/SureChEMBL_map_20141231.txt")
    pickle.dump(d, open("Data/occurance_dict_20141231.p", "wb")) #Save update

    #Update mapping as necessary with subsequent data updates
    files = os.listdir("Data/SureChemblMAP/")
    #Sort files by ascending order, remove 1st data file (already read in)
    files.sort()
    files = files[1:]
    for f in files:
        print("Analyzing:", f)
        sub_d = get_data("Data/SureChemblMAP/" + f)
        d = {**d, **sub_d}
        #Save each updated dictionary (cumulative of all updates)
        pickle.dump(d, open("Data/occurrence_dict_" + f[15:-4] + ".p", "wb"))

    pickle.dump(d, open("Data/occurrence_dict_FINAL.p", "wb"))

def main():
    #Goal = create an occurance graph of compounds
    occurances()


    #Data viz of occurance data (probably in jupyter)

if __name__ == "__main__":
    main()
