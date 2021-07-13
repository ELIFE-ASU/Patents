
# Patent Networks

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs)


Uses data from the SureChEMBL database to analyze
compound usage over time within patents.



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

  