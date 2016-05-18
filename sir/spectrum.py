#!/usr/env/bin python3
# -*- coding: utf-8 -*-

import networkx as nx


def spectrum_coefficient(G):
    spectrum_coefficient = {}
    for n in G.nodes():
        neighbor_list = nx.neighbors(G, n)
        neighbor_deg_list = []
        neighbor_deg_dict = {}

        for neighbor in neighbor_list:
            neighbor_deg_list.append(G.degree(neighbor))

        for deg in neighbor_deg_list:
            if deg not in neighbor_deg_dict:
                neighbor_deg_dict[deg] = 1
            else:
                neighbor_deg_dict[deg] += 1

        n1, n2 = 0, 0
        for (deg, num) in neighbor_deg_dict.items():
            n1 += deg*num
            n2 += num
        if n2:
            spectrum_coefficient[n] = (n1/n2+len(neighbor_list))/2
    return spectrum_coefficient


def spectrum_power(G):
    spectrum_power = {}
    for n in G.nodes():
        neighbor_list = nx.neighbors(G, n)
        neighbor_deg_list = []
        neighbor_deg_dict = {}

        for neighbor in neighbor_list:
            neighbor_deg_list.append(G.degree(neighbor))

        for deg in neighbor_deg_list:
            if deg not in neighbor_deg_dict:
                neighbor_deg_dict[deg] = 1
            else:
                neighbor_deg_dict[deg] += 1
        power = 0
        for (deg, num) in neighbor_deg_dict.items():
            power += deg * num * num
        spectrum_power[n] = power
    return spectrum_power
