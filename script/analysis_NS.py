#!/usr/env/bin python3
# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx
import matplotlib.font_manager as font_manager
ch_font = font_manager.FontProperties(fname='C:/Windows/Fonts/msyh.ttc')
import matplotlib.pyplot as plt
from SIR import *
import spectrum


# def run_SIR_simulation(G, beta, gammar, seeds, iterations=100):
#     t_max_list = []

#     for i in range(iterations):
#         sim = SIR(G, infection_rate=beta, recovery_rate=gammar, infection_seeds=seeds)
#         sim.simulate()
#         max_t = sim._get_max_t()
#         t_max_list.append(max_t)
#     return sum(t_max_list)/iterations


def get_top_n(dict, n):
    return sorted(dict.items(), key=lambda i: i[1], reverse=True)[:10]



def run_SI_simulation(G, beta,seeds):
    sim = SIR(G, infection_rate=beta, recovery_rate=0, infection_seeds=seeds)
    sim.simulate()
    i, t = sim.get_i_of_t()
    return i, t

NS = nx.read_gml('netscience.gml')
i_1, t_1  = run_SI_simulation(NS, 0.6, 'BARABASI, A')
i_2, t_2  = run_SI_simulation(NS, 0.6, 'NEWMAN, M')
i_3, t_3  = run_SI_simulation(NS, 0.6, 'CAGNEY, G')
plt.plot(t_1, i_1, 'r', label='BARABASI, A')
plt.plot(t_2, i_2, 'g', label='NEWMAN, M')
plt.plot(t_3, i_3, 'b', label='CAGNEY, G')
plt.legend(loc=0)
plt.show()

