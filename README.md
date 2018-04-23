Portrait Divergence
===================

Python code for computing **portrait divergences**, a general-purpose tool for comparing networks.

Please see the paper for more details:

**An information-theoretic, all-scales approach to comparing networks**  
James P. Bagrow and Erik M. Bollt, 2018  
[arXiv:1804.03665](https://arxiv.org/abs/1804.03665)


## Usage

The output of each calculation is a float between 0 and 1 describing how similar the two
networks are (0 = identical, 1 = maximally different).

Portrait divergences can be computed at the command line or within Python scripts:

#### Command line

1. Basic example:  
    `python portrait_divergence.py data/net1.edgelist data/net2.edgelist`

1. Directed networks stored in graphml files:  
    `python portrait_divergence.py -d --graphml digraph_time1.graphml digraph_time2.graphml`

1. Use C++ code (assuming it's installed):  
    `python portrait_divergence.py --cpp big_g.edgelist big_h.edgelist`

1. Weighted graphs (`strength` edge attribute) w 10-percentile bins on Dijkstra's paths lengths:  
    `python portrait_divergence.py --weighted=strength -b 10 --graphml g.graphml h.graphml`


See the help string for more: `python portrait_divergence.py -h`

#### Python

Here's a script to compare an 
[Erdős-Rényi](https://en.wikipedia.org/wiki/Erdős–Rényi_model) graph
and a 
[Barabási-Albert](https://en.wikipedia.org/wiki/Barabási–Albert_model) graph:

```Python
import networkx as nx
from portrait_divergence import portrait_divergence

G = nx.erdos_renyi_graph(100, 3/99)
H = nx.barabasi_albert_graph(100, 3)

Djs = portrait_divergence(G, H)
print("Djs =", Djs)
```


## Requirements

* [Python 3.x](https://www.python.org) with packages:
    + [numpy](http://numpy.scipy.org/)
    + [scipy](http://www.scipy.org/)
    + [networkx](https://networkx.github.io)

A recent install of [Anaconda Python](https://www.anaconda.com) should come with everything you need.
