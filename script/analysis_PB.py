#!/usr/env/bin python3
# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx
import matplotlib.font_manager as font_manager
ch_font = font_manager.FontProperties(fname='C:/Windows/Fonts/msyh.ttc')
import matplotlib.pyplot as plt
from SIR import *
import spectrum


def run_SIR_simulation(G, beta, gammar, seeds, iterations=100):
    t_max_list = []

    for i in range(iterations):
        sim = SIR(G, infection_rate=beta, recovery_rate=gammar, infection_seeds=seeds)
        sim.simulate()
        max_t = sim._get_max_t()
        t_max_list.append(max_t)
    return sum(t_max_list)/iterations




def get_top_n(dict, n):
    return sorted(dict.items(), key=lambda i: i[1], reverse=True)[:10]


PB = nx.read_gml('polbooks.gml')
idx = []
power_list = []
t_max_list = []
for (id,n) in enumerate(PB.nodes()):
    t_max = run_SIR_simulation(PB, 0.04, 0.00, n)
    idx.append(id+1)
    t_max_list.append(t_max)
    power_list.append(spectrum.spectrum_power(PB)[n])
    print(id+1, t_max, spectrum.spectrum_power(PB)[n])

plt.plot(power_list, t_max_list, 'k.')
# plt.xlabel('Rank')
# plt.ylabel('Degree')
# plt.ylim(1,max(degree_sequence)+1)
# plt.xlim(.9,10001)
plt.show()



# # plot
# plt.plot(t1, i1, 'r', label='度数前十节点')
# plt.plot(t2, i2, 'b', label='介数前十节点')
# plt.plot(t3, i3, 'g', label='传播功率前十节点')
# plt.xlabel("感染天数")
# plt.ylabel("感染比例")
# plt.title("网络科学科研协作关系网络SIR仿真图")
# plt.legend(loc=0)
# plt.show()

