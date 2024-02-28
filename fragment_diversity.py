import igraph as ig
import ast
import re
import itertools
from tqdm import tqdm
import os
import pickle
import pandas as pd
from random import sample


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

# Depth-first search (DFS). Choose a vertex, go to neighbors, collect info, repeat for neighbors...
def explore_component_dfs(graph, vertex, visited, component_vertices, component_edges, component_vertex_colors, component_edge_colors):

    vertex_index = vertex.index
    visited.add(vertex_index)
    component_vertices.append(vertex_index) 
    
    # Grab vertex colors
    vertex_color = graph.vs[vertex_index]["color"]
    component_vertex_colors.append(vertex_color)

    for neighbor_index in graph.neighbors(vertex):
        neighbor_vertex = graph.vs[neighbor_index]
        if neighbor_index not in visited:
            component_edges.append((vertex_index, neighbor_index))
    
            # Grab edge colors
            edge_id = graph.get_eid(vertex_index, neighbor_index) # Get the edge ID given indices of endpoints
            edge_color = graph.es[edge_id]["color"]
            component_edge_colors.append(edge_color)

            # Sweet, sweet recursivity 
            explore_component_dfs(graph, neighbor_vertex, visited, component_vertices, component_edges, component_vertex_colors, component_edge_colors)


# Collecting graph properties of a component
def dfs(graph, start_vertex):
    visited = set()
    component_vertices = []
    component_edges = []
    component_vertex_colors = []
    component_edge_colors = []

    explore_component_dfs(graph, start_vertex, visited, component_vertices, component_edges, component_vertex_colors, component_edge_colors)

    return component_vertices, component_edges, component_vertex_colors, component_edge_colors


# From a graph component, create an igraph object
def create_igraph_component(component_vertices, component_edges, component_vertex_colors, component_edge_colors):
    
    # Create a new igraph object
    g = ig.Graph()

    # Add vertices to the graph
    g.add_vertices(component_vertices)

    # Because vertices reset starting at 0 when new graph is made, 
    # map original indices from parent graph to new indices, starting at 0. 
    vertex_ids = g.vs['name']
    indices = range(len(vertex_ids))
    vertex_dict = dict(zip(vertex_ids, indices))

    # Convert component edges with new vertices
    new_component_edges = [(vertex_dict[edge[0]], vertex_dict[edge[1]]) for edge in component_edges]
    
    # Add edges to the graph
    g.add_edges(new_component_edges)

    # Set vertex colors
    g.vs["color"] = component_vertex_colors

    # Set edge colors
    g.es["color"] = component_edge_colors

    return g

# Gather all graphs in the remnant disjointed graph
def get_remnant_fragments(graph):
    og_components = []
    remnant_fragments = []
    visited = set()

    for vertex in graph.vs:
        if vertex.index not in visited:
            component_vertices, component_edges, component_vertex_colors, component_edge_colors = dfs(graph, vertex)
            og_components.append((component_vertices, component_edges, component_vertex_colors, component_edge_colors))
            component_igraph = create_igraph_component(component_vertices, component_edges, component_vertex_colors, component_edge_colors)
            remnant_fragments.append(component_igraph)
            visited.update(component_vertices)
            
    return remnant_fragments


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


def check_iso(candidate_fragments, all_frags, frag_count):
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
        for i in range(len(all_frags)):
            if possible_f.isomorphic_vf2(all_frags[i],
                                         color1=possible_f.vs["color"],
                                         color2=all_frags[i].vs["color"],
                                         edge_color1=possible_f.es["color"],
                                         edge_color2=all_frags[i].es["color"]):
                #if candidate is isomorphic, it is not unique, therefore break off check
                is_unique = False

            if not is_unique:
                frag_count[i] += 1
                break

        #If it makes it through the full list of fragments, it is unique!
        if is_unique:
            all_frags.append(possible_f)
            frag_count[len(all_frags) - 1] = 1  #Start count at 1

    #Return updated full list
    return all_frags, frag_count


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


def get_files(fp, labels):
    """ Find all file paths to specific log files to analyze

    These all correspond to fully completed MAs

    Args:
        fp (str): filepath to overall AssemblyValue directory
        labels (list): all labels of completed MAs

    Returns:
        (list): list of all filepaths
    """
    possible_dataSources = ["AuthorCpds_Done/", "AssigneeCpds_Done/"]

    fps = []
    for label in labels:
        if "SCHEMBL" in label:
            for dataSource in possible_dataSources:
                if os.path.exists(fp + dataSource + label + ".txt"):
                    fps.append(fp + dataSource + label + ".txt")
                    break
        elif label.endswith("_new"):
            fps.append(fp + "NewDatabase_Done/" + label[:-4] + ".txt")
        elif label.endswith("_full"):
            fps.append(fp + "FullDatabase_Done/" + label[:-5] + ".txt")

    return fps


def main():
    ### Using all completed MAs to build fragments
    MA_df_completed = pd.read_csv(
        "AssemblyValues/ALLSAMPLEDcpds_AssemblyGo_COMPLETED.csv")

    months = build_month_increments(1976, 2022)

    all_IDs = {}
    vscolor_map = {}
    escolor_map = {}

    for month in months:
        new_frag_count = 0
        print("--- Analyzing", month, "---")

        # #Old newdatabase code
        # files = [x for x in os.listdir(fp) if x.startswith(month)]
        # files = [x for x in files if x.endswith(".txt")]
        # print(len(files))

        fragments = []
        all_frags = []
        frag_count = {}

        ## Sample month compounds
        month_df = MA_df_completed[
            MA_df_completed["earliest_date"].str.startswith(month)]
        n = 1000
        if len(month_df) > n:
            month_df = month_df.sample(n)

        files = get_files("AssemblyValues/", list(month_df["label"]))

        #Store ids (for agave parallelization)
        all_IDs["month"] = files

        for file in tqdm(files):
            if os.path.exists(file):
                with open(file) as f:
                    lines = f.readlines()
                    lines = [l.strip() for l in lines]

            else:
                pass

            #Isolate fragments lines
            for i in range(len(lines)):
                if lines[i] == "======":
                    if lines[i + 1].startswith("Vertices"):
                        #i+1 is the start of the fragment definition, i+5 is the end
                        fragments.append(
                            build_fragment(i + 1, i + 5, lines, vscolor_map,
                                           escolor_map))

                if lines[i] == 'Remnant Graph':
                    # i+1 is the start of fragment definition,
                    remnant_graph = build_fragment(i + 1, i + 5, lines, vscolor_map, escolor_map)
                    fragments.extend(get_remnant_fragments(remnant_graph))

            ## More testing - check isomorphism within a single output file
            all_frags, frag_count = check_iso(fragments, all_frags, frag_count)

        #To ensure these are the same values
        print("Num Frags:", len(all_frags))
        print("Num Frags Count", len(frag_count))

        print("\n")

        pickle.dump(all_frags,
                    file=open(
                        "AssemblyValues/Fragments/fullFrags_" + month +
                        "_withRemnants.p", "wb"))
        pickle.dump(frag_count,
                    file=open(
                        "AssemblyValues/Fragments/fullFragsCount_" +
                        month + "_withRemnants.p", "wb"))

    pickle.dump(
        all_IDs,
        file=open(
            "AssemblyValues/Fragments/all_ids_updated_withRemnants.p",
            "wb"))
    pickle.dump(vscolor_map,
                file=open(
                    "AssemblyValues/Fragments/vscolor_map_withRemnants.p",
                    "wb"))
    pickle.dump(escolor_map,
                file=open(
                    "AssemblyValues/Fragments/escolor_map_withRemnants.p",
                    "wb"))

    print(vscolor_map)
    print(escolor_map)
    print("Unique fragment size:", len(all_frags))
    # for frag in all_frags:
    #     print(frag.vs["color"], frag.es["color"])
    #     print("----")


if __name__ == "__main__":
    main()
