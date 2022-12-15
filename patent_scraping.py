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


def get_author_patents(name):
    """ Downloads a JSON of all patent records associated with a specific author

    NOTE: TAKES ~24 HOURS FOR FULL DATASET

    Args:
        name (str): name of patent author to search
    """
    #HREF is split into two parts, in order to add name in the middle
    href_1 = "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=json&query={%22download%22:%22*%22,%22collection%22:%22patent%22,%22where%22:{%22ands%22:["
    href_2 = "]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000}"

    #Names needs to be split & added partwise
    names = name.split()
    for subname in names:
        #Also surrounded by these specific characters (idk why...ask PubChem lol)
        href_1 += "{%22*%22:%22" + subname + "%22}"
        if subname != names[-1]:
            href_1 += ","

    #Merge full href, get resulting record (if something fails along the way, just skip it)
    href = href_1 + href_2
    try:
        get_url = request.urlopen(href)

        #Save record as JSON in Data/Patents
        with open("Data/Patents/" + name.replace(" ", "_") + ".json", "w") as f:
            json_data = json.loads(get_url.read())
            json.dump(json_data, f, indent=4)
    except:
        pass


def get_assignee_patents(name):
    """ Downloads a JSON of all patent records associated with a specific assignee

    NOTE: takes ~10 hours with list of assignee names

    Args:
        name (str): name of patent author to search
    """
    href_1 = "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=json&query={%22download%22:%22*%22,%22collection%22:%22patent%22,%22where%22:{%22ands%22:[{%22*%22:%22"
    href_2 = "%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000}"

    #Getting exact name - probably should do this for authors too...
    name = " ".join(name.split()[:-1])
    name = name.replace(" ", "%20")

    href_1 += name
    href = href_1 + href_2

    #Merge full href, get resulting record (if something fails along the way, just skip it)
    try:
        get_url = request.urlopen(href)

        #Save record as JSON in Data/Patents/Patent_Assignee_Records, with underscores in the name
        with open(
                "Data/Patents/Patent_Assignee_Records/" +
                name.replace("%20", "_") + ".json", "w") as f:
            json_data = json.loads(get_url.read())
            json.dump(json_data, f, indent=4)
    except:
        pass


def main():
    # #### Download n random patent IDs ####
    # fp = "Data/CpdPatentIdsDates/patent_ID_index_dict.p"
    # n = 10000

    # patents = get_patentIDs(fp, n)

    # for patent in tqdm(patents):
    #     get_patentData(patent)

    # #### Download patents associated with authors from above n patents ####
    # names = pickle.load(file=open("Data/Patents/authors.p", "rb"))

    # for name in tqdm(names):
    #     get_author_patents(name)

    #### Download patents assiciated with assignees from above n patents ####
    names = pickle.load(file=open("Data/Patents/assignees.p", "rb"))

    for name in tqdm(names):
        get_assignee_patents(name)


if __name__ == "__main__":
    main()
