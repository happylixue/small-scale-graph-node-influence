#!/usr/env/bin python3
# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.pyplot as plt
import spectrum
from SIR import *
import numpy as np

SW = nx.davis_southern_women_graph()
# nx.draw(G, with_labels=True)
# plt.show()

women = [n for n in SW.nodes() if  SW.node[n]['bipartite'] == 0]
event = [n for n in SW.nodes() if  SW.node[n]['bipartite'] == 1]


# for n in SW.nodes():
#     if n in event:
#         print(n, spectrum.spectrum_coefficient(SW)[n],spectrum.spectrum_power(SW)[n])


def run_SI_simulation(G, beta, seeds, iterations=100):
    num_of_idx = G.order()
    sim_i, sim_t = np.zeros(num_of_idx), np.zeros(num_of_idx)

    for i in range(iterations):
        sim = SIR(G, infection_rate=0.04, recovery_rate=0.00, infection_seeds='E7')
        sim.simulate()
        i, t = sim.get_i_of_t()
        sim_i += i
        sim_t += t
    return np.array([sim_i/iterations, sim_t/iterations])


i_1, t_1  = run_SI_simulation(SW, 0.06, 'E7')
i_2, t_2  = run_SI_simulation(SW, 0.06, 'E8')
i_3, t_3  = run_SI_simulation(SW, 0.06, 'E9')
plt.plot(t_1, i_1, 'r', label='E7')
plt.plot(t_2, i_2, 'g', label='E8')
plt.plot(t_3, i_3, 'b', label='E9')
plt.legend(loc=0)
plt.show()


