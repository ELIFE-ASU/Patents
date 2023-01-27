import igraph as ig
import ast
import re
import itertools
from tqdm import tqdm
import os


def get_list(line):
    """ Evaluation a line to find the description and list containing data
        for eventual addition to an igraph graph

    Args:
        line (str): specific line to evaluate

    Returns:
        list: 2-element list of [0]: description (vertices, edges, vertexcolors, or edgecolors) 
            and [1]: list of values matching that description
    """
    #Split each line on the end of the first word
    line = line.split("s ")
    #values = json.loads(line[1].replace(" ", ", "))
    try:
        #Lists of ints (vertices & edges) can be treated with literal_eval
        values = ast.literal_eval(line[1].replace(" ", ", "))
    except ValueError:
        #Lists of strings/chars need to be treated with regex (for some reason)
        values = re.sub('[\[\]]', '', line[1])
        values = values.split()

    return [line[0] + "s", values]


def build_fragment(start, end, lines, vscolor_map, escolor_map):
    """ Builds an igraph graph of a single fragment

    Args:
        start (int): line number of the vertices entry
        end (int): line number of the edgecolours entry
        lines (list): list of all text file lines

    Returns:
        igraph object: graph containing vertex, edges, and vs/es colors
    """
    g = ig.Graph()

    #Label map from MA node labels to igraph node labels
    node_map = {}

    for line in lines[start:end]:
        #Isolate description
        desc = get_list(line)

        #Iterativley(?) build igraph fragment
        #Add vertices to igraph (guarenteed to be first)
        if desc[0] == "Vertices":
            g.add_vertices(len(desc[1]))

            #Create a map for vertex numbers to igraph values (for edges)
            for i in range(len(desc[1])):
                node_map[desc[1][i]] = i

        #Add edges based on label map built in vetex
        elif desc[0] == "Edges":
            for edge in desc[1]:
                g.add_edges([(node_map[edge[0]], node_map[edge[1]])])

        #Add vertex & edge colors
        elif desc[0] == "VertexColours":
            for element in desc[1]:
                if element not in vscolor_map:
                    vscolor_map[element] = len(vscolor_map) + 1

            g.vs["color"] = [vscolor_map[x] for x in desc[1]]

        elif desc[0] == "EdgeColours":
            for element in desc[1]:
                if element not in escolor_map:
                    escolor_map[element] = len(escolor_map) + 1

            g.es["color"] = [escolor_map[x] for x in desc[1]]

    return g


def check_iso(candidate_fragments, all_frags):
    """ Finds the unique fragments within a set of candidates

    Args:
        candidate_fragments (list): list of igraph graphs from a single text file
        all_frags (list): list of all unique igraph graphs from previous files

    Returns:
        list: updated list of unique igraph graphs
    """
    #Check all candidate fragments...
    for possible_f in candidate_fragments:
        is_unique = True
        #...against all unique fragments
        for f in all_frags:
            if possible_f.isomorphic_vf2(f,
                                         color1=possible_f.vs["color"],
                                         color2=f.vs["color"],
                                         edge_color1=possible_f.es["color"],
                                         edge_color2=f.es["color"]):
                #if candidate is isomorphic, it is not unique, therefore break off check
                is_unique = False
                break

        #If it makes it through the full list of fragments, it is unique!
        if is_unique:
            all_frags.append(possible_f)

    #Return updated full list
    return all_frags


def main():
    ### Testing - build iGraph fragments from a single output file
    fp = "Data/AssemblyValues/NewDatabase_Done/2019-12_211687.txt"

    ### Testing part 2 - all 1980 values
    fp = "Data/AssemblyValues/NewDatabase_Done/"
    files = [x for x in os.listdir(fp) if x.startswith("1980-02")]

    fragments = []
    all_frags = []
    vscolor_map = {}
    escolor_map = {}

    for f in tqdm(files[0:2000]):
        with open(fp + f) as f:
            lines = f.readlines()
            lines = [l.strip() for l in lines]

        #Isolate fragments lines
        for i in range(len(lines)):
            if lines[i] == "======":
                if lines[i + 1].startswith("Vertices"):
                    start = i + 1
                    end = i + 4

                    fragments.append(
                        build_fragment(i + 1, i + 5, lines, vscolor_map,
                                       escolor_map))

        ## More testing - check isomorphism within a single output file
        all_frags = check_iso(fragments, all_frags)

    print(vscolor_map)
    print(escolor_map)
    print("Unique fragment size:", len(all_frags))
    # for frag in all_frags:
    #     print(frag.vs["color"], frag.es["color"])
    #     print("----")


if __name__ == "__main__":
    main()
