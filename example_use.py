#!/usr/bin/env python
# -*- coding: utf-8 -*-

# example_use.py
# Jim Bagrow
# Last Modified: 2018-04-22

import sys, os
import itertools
import networkx as nx
import numpy as np
from portrait_divergence import portrait_divergence


# make n ER graphs and n BA graphs:
n = 10
list_ER = [ nx.erdos_renyi_graph(100, 3/99)  for _ in range(n) ]
list_BA = [ nx.barabasi_albert_graph(100, 3) for _ in range(n) ]


# compare every pair of ER graphs:
Djs_sameER = []
for ERi, ERj in itertools.combinations(list_ER, 2):
    Djs = portrait_divergence(ERi, ERj)
    Djs_sameER.append(Djs)

# compare every pair of BA graphs:
Djs_sameBA = []
for BAi, BAj in itertools.combinations(list_BA, 2):
    Djs = portrait_divergence(BAi, BAj)
    Djs_sameBA.append(Djs)

# compare every ER with every BA:
Djs_ERvBA = []
for ER, BA in itertools.product(list_ER, list_BA):
    Djs = portrait_divergence(ER, BA)
    Djs_ERvBA.append(Djs)


try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit(0)

# plot histograms:
hargs = dict(bins='auto', density=True, histtype='stepfilled')
plt.hist(Djs_sameER, label='Same ER',   alpha=0.7, **hargs)
plt.hist(Djs_sameBA, label='Same BA',   alpha=0.6, **hargs)
plt.hist(Djs_ERvBA,  label='ER vs. BA', alpha=0.7, **hargs)

plt.xlabel("Portrait divergence $D_\mathrm{JS}$")
plt.ylabel("Prob. density")
plt.legend()
plt.tight_layout()
plt.show()
