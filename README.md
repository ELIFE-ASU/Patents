
# Patent Networks

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs)


Uses data from the SureChEMBL database to analyze
compound usage over time within patents.

## Code

This code performs a variety of tasks relating to network & molecular assembly analysis of the SureChemBL database. The various scripts are organized here by tasks - individual files & methods have comments within them describing their exact usages.

### Chemical Properties (MA, MW, bonds, etc...) Visualizations

- `MA_distribution.ipynb`

Visualizes MA & Molecular Weight distribution (x-axis: MA, y-axis: count) for compounds sampled from different 5-year intervals. 

- `assembly_analysis.ipynb`

Workhorse jupyter notebook containing various ideas & graphs for visualizing MA data. These visualizations mainly resulted in the chemical properties & social factors figures in the manuscript, as well as a number of SI figures. The topics found in this notebook include: New Compound MA vs Full Database MA (and various tests, including an expected MA from molecular weight, wasserstein tests between new & full compounds, and normal distribution tests for each month of sampling); MA increases across all types of sampling; MA distributions over time intervals (similar to `MA_distribution.ipynb`); Molecular Weight time series analyses; Molecular Weight & MA spearman correlations & plots, including exact & inexact MA values; Bonds & MA spearman correlations and plots, also including exact & inexact MA values; MA predictions using ML techniques, particularly exponential smoothing; MA & Degree (usage) correlations; and preliminary fragment analysis.

- `fragment_diversity.ipynb`

Tests various fragment hypothesis - number of unique fragments over time (more detailed figure is in chemical properties manuscript figure), different types of fragments over time, time-ordering of fragment analyses (forward in time, reversed, and mixed) which tests if innovative fragments are based on historical contigencies, and fragment/copy number log/log plots. Also has some potential tests (e.g., relative frequency of high-bond count fragments) at the end which could be interesting.

### Chemical Properties (MW, MW, bonds, etc...) Analysis & Data Collection

- `assemblyCalcs_cheminformatics.py`

_Outdated_ Calculates molecular weight from sampled compounds (MA is already calculated for these compounds).

Uses RDKit MolDescriptors library to calculate molecular weight

Outdated because there are more integrated methods of calculating molecular weight (along with bonds & other chemical descriptors) in other scripts.

- `assemblyCalcs_degrees.py`

Links MA values of compounds with network degree. These data are used in the social factors figure, where degree represents usage. The hypothesis is that degree (usage) has an influence on the MA of the compound. Also links different quantiles of degrees (e.g., 90th percentile of degrees) with MA.

- `assemblyCalcs_percentiles.py`

** Outdated ** Calculates MA values of compounds with high preferential attachment value. The hypothesis was that high-attachment compounds have lower MA, but this line of questioning was replaced by a degree analysis. Also uses the old Monte-Carlo based MA calculation instead of the updated AssemblyGo package.

- `assemblyCalcs_samples.py`

** Main Assembly Calculation Script ** Calculates MA values in parallel using AssemblyGo package with a timeout of 300 seconds. This script was used to calculate MA values over Full samples, New samples, cost data, author/assignee samples, and other data sources.

- `assembly_analysis.py`

Links MA calculations with molecular weight & bonds. 

- `build_mol_files.py`

Takes InChI representations of data (standard SureChemBL storage) and converts them into molecular object files (.mol, also can used as .sdf), which is necessary for AssemblyGo MA calculations. This script also takes completed MA .txt files and moves them into the correct directory structure on Agave.

- `cpd_similarity.ipynb`

Calculates the similarity between compounds invented within different 5-year time intervals. This was not used in the manuscript, but was based on the hypothesis that compounds invented near each other are more similar than those invented further away.

- `fragment_diversity.py`

Workhorse script which: finds fragments from MA .txt result files, builds fragments into iGraph objects, and determines if a given fragment is a unique fragment (given a specific time-ordering - see `fragment_diversity.ipynb` for different sequencing) using iGraph's isomorphic_vf2 test.

- `fragment_diversity_parallel.py`

_Outdated_ Failed attempt to parallelize `fragment_diversity.py` - abandoned due to sequencing, since each sampled set of MAs needed to have a specific ordering.

- `get_MA_results.py`

Takes .txt output files from the AssemblyGo algorithm and stores the label (identifying number/ID of the compound, listed as the first part of the filename), MA value, and time it took for the MA algorithm to complete as a .csv file.


### Cost Visualization & Analysis

- `cost_influence.ipynb`

Uses cost data (collected by Hessam Mehr & Dario Caramelli) & MA values to build correlations between price & MA. Also includes time-series analysis of cost over time. 

- `cost_timewise_influence.py`

Script to merge cost compounds (sampled separately from SureChemBL) with SureChemBL time-series data to build a time-series analysis (results & vizualization found in `cost_influence.ipynb`)


### Network Visualizations

- `bipartite_network_analysis.ipynb`

** Outdated ** Basic plots describing network stats (average degree, largest connected component size, number of edges, average clustering coefficient) for SureChemBL networks over time. Updated & better plots can be found in `network_analysis.ipynb` 

- `network_analysis.ipynb`

Plots describing bipartite network statistics, such as degree distributions & powerlaw fits, avg degree of both compounds and patents, sizes of compounds & patents and exponential fits, and largest connected compound sizes. Many of these plots are found in the SI, and the network stats manuscript figure is based on ideas first graphed here.

- `occurance_viz.ipynb`

_Outdated_ Basic plots representing the occurance of compounds within patents - these were an attempt to see if compounds followed a "rich-get-richer" model, which was found using preferential attachment instead.

- `pref_attachment_analysis.ipynb`

_Outdated_ Initial attempts to analyze & graph preferential attachment values, given 5-year intervals. Updated code is in `pref_attachment_noNetworks.ipynb`

- `pref_attachment_noNetworks.ipynb`

Calculates & graphs preferential attachment values for compounds over 5-year intervals without loading iGraph networks (uses dictionaries instead). The final version of these graphs are in the preferential attachment manuscript figure. Also tracks individual compound changes over time (e.g., green chemistry compounds) and how preferential attachment changes for specific percentiles of compounds (e.g., how the 99th percentile of compounds changes over time). 

### Network Analysis

- `get_bipartite_network_data.py`

Calculates network statistics across the cpd-patent SureChemBL bipartite network by month, including number of nodes, edges, cpds, patents,  average degrees of both, and the largest connected component size. These data are used to build the network statistics manuscript figure.

- `get_cpd_network_data.py`

Builds & calculates network statistics for a compound-only SureChemBL network, which are individually built for each month. Statistics include nodes, edges, average & max degree, as well as largest connected component size and average clustering coefficient.

- `network_analysis.py`

Finds both new IDs (IDs which were previously not in the network previously - these are used as the basis for the New ID sampling for MA calculations) and largest connected component IDs.

- `pageRank.ipynb`

Basic plots & descriptions of PageRank analysis - this was an attempt to track innovative compounds over time, with the hypothesis of comopunds that have a high connectivity (or a highly variable, positive connectivity) are innovative.

- `pageRank.py`

Calculates PageRank statistics for all compounds & patents within the SureChemBL networks over months.


### SureChemBL Data

- `cpd_analysis.py`

Samples compounds over months from novel compounds (`sample_compounds_unique()` method) and full database (`sample_compounds` method). Also includes various tests to find the highest degree compounds (see `assemblyCalcs_degrees.py`) and compounds from Llanos et al, 2019, as a hypothesis as to where specific compounds were listed in the database.

- `sample_sanityTesting.ipynb`

Testing to make sure SureChemBL sampling was accurate.

### Societal Factors Visualizations

- `patent_MA_analysis.ipynb`

Workhorse notebook, analyzes changes in MA across authors, assignees, and classifications - also organizes and samples from patent author/assignee/classification lists, as well as performs dropout tests. Various plots from here are found in the SI, and the delta-MA plots are the initial versions of the ones found in the societal factors manuscript figure.

### Societal Factors Analysis

- `patent_authors_assignees.ipynb`

Notebook which scrapes lists of patents for authors & assignees, as well as links patent IDs to compounds. Additionally links compounds with structures and prepares patent data for MA analysis (the results of which are found in `patent_MA_analysis.ipynb`)

- `patent_prefAttach_analysis.ipynb`

_Outdated_ Attempts to link preferential attachment values of compounds to compounds associated with patents & assignees over time. This was replaced by MA values instead.

- `patent_scraping.py`

Scrapes PubChem, using url-based JSON downloads, to find patents associated with randomly selected authors & assignees.

- `reaxys_influence.ipynb`

Initial attempts to use the Reaxys database as a way to measure societal impact on MA changes. Various metrics, such as the number of times a compound is used as a product or reactant, are correlated with MA. This was dropped in favor of the authors/assignees/classifications analysis found in the manuscript.

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

  