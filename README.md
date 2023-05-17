
# Patent Networks

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs)


Uses data from the SureChEMBL database to analyze
compound usage over time within patents.

## Data

The data for this project (will be) stored on Cradle in the `SureChemBL` directory, as well as backed up on Alice Lab Mac. This describes the file structure of all data generated from SureChemBL, as well as descriptions. 

### Cpd_Data

- Raw .txt compound files from SureChemBL, containing a SureChemBL ID, SMILES, InChI, and InChIKey representation of compounds in a tab-separated csv format. These are updated quarterly from ([SureChemBL](https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/data/map/)), and if additional data is needed, it must be downloaded and manually added to Cradle.

- `SureChemBL_allCpds.p`: pandas dataframe of all combined data from .txt files.

- `attachment_values_allIDs.p`: all preferential attachment values (in 5 year increments?) associated with each SureChemBL ID.

- `cpd_ID_index_dict.p`: dictionary linking ID to SureChemBL index (CHECK)

- `master_cpd_date_df.p`: pandas dataframe links cpd IDs to dates of appearaence (CHECK)

- `master_cpd_date_dict.p`: dictionary linking cpd IDs to dates of appearaence (CHECK)

- `master_cpd_date_index_df.p`: pandas dataframe linking cpd ID, dates, and SureChemBL index (CHECK)

### AssemblyValues

Contains all assembly & fragment data - some old data (specifically, 1000 sampled compounds per month from both the full database and new compounds) from earlier iterations of the assembly calculator are in the `Old Assembly Results` directory.

- `AssigneeCpds_Done`: XX compounds, from patents associated from a random set of YY patent assignees. The name of every .txt file is the SureChemBL Id of the compound.

- `AuthorCpds_Done`: XX compounds, from patents associated from a random set of YY patent authors. The name of every .txt file is the SureChemBL Id of the compound.

- `FullDatabase_Done`: XX compounds, randomly sampled from the full set of compounds available in the SureChemBL database at every month. Labeled by the sampled month.

- `NewDatabase_Done`: XX compounds, randomly sampled from the set of compounds which were new (not previously entered) to theSureChemBL database at every month. Labeled by the sampled month.

- `CostRandom`: 50,000 compounds, from a subset of Reaxys which Hessam Mehr & Dario Caramelli were able to extract cost data from.

- **Fragments**: Fragment data - contains all sampled IDs (including file paths) (`all_ids.p`), test fragment samples from the author MA compounds only, as well as full fragments (fullFrags*) and counts (fullFragsCounts*) of fragments from all completed data, regardless of intitial sampling. Also includes two dictionaries mapping iGraph color values to elements & bonds (`vscolor.p` and `escolor.p`)

### Attachment_noNetworks

Contains preferential attachment values & compound degrees for all 5-year increments. These values were calculated using degrees only, with no iGraph networks involved.

### Cost

Summary data, including .csv files with MAs & dates associated with all 50,000 sampled compounds from Hessam's & Dario's sampling.

### CpdPatentIdsDates

This directory holds connections between compounds & patents, as well as network groupings (e.g., largest connected component IDs) and unique compounds & patents.

- **Cpd_Date_Dict**: For every month, a dictionary linking compounds (SureChemBL IDs) to date of appearence in a patent. Overlaps are allowed.

- **Cpd_Patent_Edges**: For every month, a dictionary linking compounds (SureChemBL IDs) to patents they appear in. Overlaps allowed.

- `index_edgelist_bipartite.p`: Network edgelist of entire graph (CHECK).

- **New_Ids**: For every month, lists all novel compound IDs. These are compounds which have not been seen in the database at the month of first appearence - this does not necessarily mean invention.

- **New_LCC_Ids**: For every month, lists all new compound IDs which have been newly added to the largest connected cluster of the overall SureChemBL bipartite network. These could be new IDs, or "old" IDs which are newly connected to patents within the LCC.

- **Patent_Cpd_Edges**: For every month, a dictionary linking all patents registered to compounds (SureChemBL IDs) used within them.

- **Patent_Date_Dict**: For every month, a dictionary linking all patents registered to the exact date of registration.

- **Patent_ID_Edges**: For every month, a dictionary from the SureChembBL network linking all patent iGraph IDs to the compound iGraph IDs connected to to (CHECK).

- `patent_ID_index_dict`: CHECK

- **Unique_Cpds**: For every month, a list of all unique compounds used in patents that month (the set of all compounds used).

- **Unique_Patents**: For every month, a list of all unique patents registered that month.

### Degrees

Lists of network degrees (number of connections) of various components of the SureChemBL iGraph network.

- `full_id_degrees_XXXX_XXXX.p`: Pickle of all degrees of all nodes (represented by integer ID values) within the network.

- **Months**: Contains degrees of compounds only, patents only, and all nodes (IdDegrees, POSSIBLE DEGREES) within every month. The network grows every month, so these degrees increase over time (connections are not broken).

### Graphs

Graph files (full & month-by-month) containing iGraph descriptions of the SureChemBL network. The network is bipartite, with compounds and patents both considered different types of nodes, with undirected edges between patents & compounds if a compound is found within a patent.

*Note: Needs to be extended to include recent 2020-2022 data*

- `cpd_patent_G_FULL.p`: Full iGraph network.

- `cpd_patent_G.p`: CHECK

- **G**: For every month, the current bipartite iGraph network is stored. These contain only patents (& compounds) registered up through the given month.

- **G_Cpd**: For every month, the current compound-only iGraph network is stored. This is the unipartite projection of the bipartite network stored in **G**.

- `index_edgelist_bipartite.p`: CHECK

### Network Stats

Contains various network statistics calculated from the graphs in **Graphs**. 

- `lcc_stats.csv`: month-by-month largest connected cluster statistics, including total ID count (in LCC and overall network), number of IDs added (from outside of LCC & newly added to the network), and percentages of each.

- `networkStats_byMonth.csv`: month-by-month network stats, including number of nodes (including patent & cpd), edges, average degree (including patent & cpd), and LCC size.

- **PageRank**: for every month, page rank analysis of all nodes to rank importance.

- `patent_stats.csv`: for every month, gives the number of total patents & new patents.

- `stats_XXXX_XXXX.csv`: for every month in each 5-year period between 1980-2020, gives number of nodes, edges, avg degree, max degree, LCC size, and clustering coefficient.

### Pref_Attachment

Holds preferential attachment values for each 5-year period between 1980-2020.

### Reaxys

Holds randomly sampled compounds (compounds only, no metadata) from Reaxys. This was originally thought of as a null model for MA analysis, but there is no verifiable time data with compounds nor a complete network, so this was not used.

- `sample_inchi_calcs.p`: CHECK


## Usage/Examples

#### Occurrence data

```python
python get_occurrences.py
```

Takes pre-downloaded data from SureChEMBL 
([available here](https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/data/map/)) 
and counts the number of times a compound appears
in a specific date range. This script assumes
data is present in a `Data\SureChemblMAP\`
directory, and writes data to `Data\Occurrences`.

These data are analyzed in ` occurance_viz.ipynb`

#### Network
```python
python build_network.py
```

Builds a bipartite, undirected network using 
[igraph](https://igraph.org/python/). Compounds
and patents are the nodes, and an edge exists 
between a compound and patent if that compound
is found within a patent. 

  