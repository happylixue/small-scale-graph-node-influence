#!/usr/env/bin python3
# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.pyplot as plt
import spectrum
from SIR import *


G = nx.Graph()
G.add_edge(1, 4)
G.add_edge(2, 4)
G.add_edge(3, 4)
G.add_edge(3, 4)
G.add_edge(4, 5)
G.add_edge(4, 6)
G.add_edge(4, 7)
G.add_edge(5, 6)
G.add_edge(6, 7)

nx.draw(G, with_labels=True)
plt.show()



def run_SIR_simulation(G, beta, gammar, seeds, iterations=500):
    t_max_list = []
    for i in range(iterations):
        sim = SIR(G, infection_rate=beta, recovery_rate=gammar, infection_seeds=seeds)
        sim.simulate()
        max_t = sim._get_max_t()
        t_max_list.append(max_t)
    return sum(t_max_list)/iterations

for n in G.nodes():
    t_max = run_SIR_simulation(G, 0.04, 0.0, n)
    print(n, t_max, spectrum.spectrum_coefficient(G)[n], spectrum.spectrum_power(G)[n])
