import urllib
import urllib.request as request
import urllib.error as error
import time
import json
import pickle
from tqdm import tqdm
from random import sample


def get_patentIDs(fp, n):
    """ Return a sample of patent ids from all patents in SureChemBL dataset

    Args:
        fp (str): filepath to patent_ID_index dictionary
        n (int): number of patents to sample

    Returns:
        (list): list of n patent IDs
    """
    patent_id_dict = pickle.load(file=open(fp, "rb"))

    return sample(list(patent_id_dict), n)


def save_patentJSON(url, patent):
    """ Save JSON for a given patent in Data/Patents/   

    Args:
        url (urllib object): result of patentAPI search
        patent (str): patent ID
    """
    patresp = json.loads(url.read().decode('latin-1'))
    #print(patresp)

    # Serializing json
    json_object = json.dumps(patresp, indent=4)

    # Writing to sample.json
    with open("Data/Patents/patent_{}.json".format(patent), "w") as outfile:
        outfile.write(json_object)


def get_patentData(patent):
    """ Download JSON with patent information

    Args:
        patent (str): patent ID
    """
    patentapi = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/patent/{}/JSON/".format(
        patent)
    try:
        url = request.urlopen(patentapi)

        save_patentJSON(url, patent)

    except error.HTTPError as err:
        print("tried {} will sleep on it".format(patentapi))
        time.sleep(5)
        try:
            url = request.urlopen(patentapi)
            save_patentJSON(url, patent)

        except error.HTTPError as err:
            print(err)

    except Exception as e:
        print(e)


def main():
    fp = "Data/CpdPatentIdsDates/patent_ID_index_dict.p"
    n = 10000

    patents = get_patentIDs(fp, n)

    for patent in tqdm(patents):
        get_patentData(patent)


if __name__ == "__main__":
    main()
