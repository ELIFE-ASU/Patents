
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

  