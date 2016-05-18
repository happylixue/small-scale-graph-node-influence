#!/usr/env/bin python3
# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.pyplot as plt
import spectrum

KC = nx.karate_club_graph()
# nx.draw(G, with_labels=True)
# plt.show()

for n in KC.nodes():
    print(n, spectrum.spectrum_coefficient(KC)[n],spectrum.spectrum_power(KC)[n])
