#!/usr/env/bin python3
# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
ch_font = font_manager.FontProperties(fname='C:/Windows/Fonts/msyh.ttc')
from SI import *


def run_SI_simulation(G, beta, seeds, iterations=500):
    number_of_index = G.order()
    total_i, total_t = np.zeros(number_of_index), np.zeros(number_of_index)
    for i in range(iterations):
        sim = SI(G, infection_rate=beta, infection_seeds=seeds,verbose=False)
        sim.simulate()
        i, t = sim.get_i_of_t()
        total_i += i
        total_t += t
    average_i = total_i/iterations
    average_t = total_t/iterations
    return np.array([average_i, average_t])

NS = nx.read_gml('netscience.gml')
i_1, i_t_1 = run_SI_simulation(NS, 0.004, 'BARABASI, A')
i_2, i_t_2 = run_SI_simulation(NS, 0.004, 'NEWMAN, M')
i_3, i_t_3 = run_SI_simulation(NS, 0.004, 'CAGNEY, G')

# plot
plt.plot(i_t_1, i_1, 'r', label='度数前10节点')
plt.plot(i_t_2, i_2, 'b', label='介数前10节点')
plt.plot(i_t_3, i_3, 'g', label='传播功率前10节点')
plt.xlabel("感染天数", fontproperties=ch_font)
plt.ylabel("感染比例", fontproperties=ch_font)
plt.title("悲惨世界人物关系网络SI仿真", fontproperties=ch_font)
plt.legend(loc=0, prop=ch_font)
# plt.savefig('SI_LM_1.png')
plt.show()



