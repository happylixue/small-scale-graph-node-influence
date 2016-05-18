#!/usr/env/bin python3
# -*- coding: utf-8 -*-

import random
import numpy as np
import networkx as nx
from Node import Node
import matplotlib.pyplot as plt


class SIR:

    def __init__(self,
                 G,
                 infection_rate,
                 recovery_rate,
                 infection_seeds=1,
                 verbose=False):

        if nx.is_directed(G):
            raise TypeError("G is directed but needs to be undirected")

        self.G = G
        self.infection_rate = infection_rate
        self.recovery_rate = recovery_rate
        self.verbose = verbose
        self.rates = np.array([infection_rate, recovery_rate])
        self.infected = set()
        self.recovered = set()
        self.nodes = set(G.nodes())
        self.Node = {n: Node() for n in self.nodes}
        self.t_count = 0

        # add infected seed nodes
        if type(infection_seeds) == int:
            seed_nodes = random.sample(self.get_susceptibles(), infection_seeds)
        elif type(infection_seeds) == str:
            seed_nodes = [infection_seeds]

        self.infected.update(seed_nodes)
        for n in seed_nodes:
            self.Node[n].set_infected()

        # process new edges for SI links
        self.SI_links = set()
        for newly_inf in self.infected:
            new_edges = [(newly_inf, n) for n in G.neighbors(newly_inf)]
            removed_SI_links, new_SI_links = self._get_removed_and_new_SI_links_from_edge_list(new_edges)
            self.SI_links.update(new_SI_links)
            self.SI_links.difference_update(removed_SI_links)

        self.t = 0.
        self.s_of_t = [[0., self.s()]]
        self.i_of_t = [[0., self.i()]]
        self.r_of_t = [[0., self.r()]]

    def _get_event_rates(self):
        return np.array([
            self.number_of_SI_links() * self.infection_rate,
            self.number_of_infected() * self.recovery_rate,
        ], dtype=float)

    def _get_removed_and_new_SI_links_from_edge_list(self, edgelist):
        new_SI = []
        removed_SI = []

        [
            new_SI.append(e)
            if ((self.Node[e[0]].is_infected() and self.Node[e[1]].is_susceptible()) or
                (self.Node[e[1]].is_infected() and self.Node[e[0]].is_susceptible()))
            else removed_SI.extend([e, (e[1], e[0])])
            for e in edgelist
        ]
        return removed_SI, new_SI

    def _recover_event(self):
        if self.verbose:
            print("============ recover event============")
        recovered = random.sample(self.infected, 1)[0]
        self.infected.remove(recovered)
        self.recovered.add(recovered)
        self.Node[recovered].set_recovered()

        deleted_edges = []
        [deleted_edges.extend([(recovered, n), (n, recovered)]) for n in self.G.neighbors(recovered)]

        if self.verbose:
            print("deleted", deleted_edges)

        self.SI_links.difference_update(deleted_edges)

    def _infection_event(self):
        if self.verbose:
            print("============= infection event============")
            print("infected:", self.infected)
            print("recovered:", self.recovered)
            print("SI link", self.SI_links)

        infective_link = random.sample(self.SI_links, 1)[0]

        if self.verbose:
            print("infective_link:", infective_link)

        if self.Node[infective_link[0]].is_infected() and self.Node[infective_link[1]].is_susceptible():
            newly_inf = infective_link[1]
            self.infected.add(infective_link[1])
        elif self.Node[infective_link[1]].is_infected() and self.Node[infective_link[0]].is_susceptible():
            newly_inf = infective_link[0]
            self.infected.add(infective_link[0])
        else:
            raise ValueError("There was a non SI-link in the array of SI links. This shouldn't happen.")

        self.Node[newly_inf].set_infected()

        # process new edges for SI links
        new_edges = [(newly_inf, n) for n in self.G.neighbors(newly_inf)]
        removed_SI_links, new_SI_links = self._get_removed_and_new_SI_links_from_edge_list(new_edges)

        if self.verbose:
            print("now infected:", self.infected)
            print("potential new_SI_links:", new_edges)
            print("new_SI_links:", new_SI_links)
            print("removed_SI_links:", removed_SI_links)

        self.SI_links.update(new_SI_links)
        self.SI_links.difference_update(removed_SI_links)

    def _choose_tau_and_event(self):
        rates = self._get_event_rates()
        total_rate = rates.sum()
        tau = np.random.exponential(1./total_rate)
        try:
            event = np.random.choice(len(rates), p=rates/total_rate)
        except ValueError as e:
            print("rates:", rates)
            print("total_rate:", total_rate)
            raise ValueError(e)

        return tau, event

    def _event(self):
        tau, event = self._choose_tau_and_event()
        self.t += tau
        self.t_count += 1

        if event == 0:
            self._infection_event()
            self.s_of_t.append([self.t, self.s()])
            self.i_of_t.append([self.t, self.i()])
        elif event == 1:
            self._recover_event()
            self.r_of_t.append([self.t, self.r()])
            self.i_of_t.append([self.t, self.i()])

    def simulate(self):
        # when there is no infected nodes or susceptibles to infected, simulation stops
        while self.number_of_infected() > 0 and self.number_of_susceptibles() > 0:
            self._event()

    def number_of_SI_links(self):
        return len(self.SI_links)

    def get_susceptibles(self):
        return (self.nodes - self.infected) - self.recovered

    def number_of_susceptibles(self):
        return self.G.number_of_nodes() - self.number_of_infected() - self.number_of_recovered()

    def number_of_infected(self):
        return len(self.infected)

    def number_of_recovered(self):
        return len(self.recovered)

    def S(self):
        return self.number_of_susceptibles()

    def I(self):
        return self.number_of_infected()

    def R(self):
        return self.number_of_recovered()

    def s(self):
        return self.number_of_susceptibles() / float(self.G.number_of_nodes())

    def i(self):
        return self.number_of_infected() / float(self.G.number_of_nodes())

    def r(self):
        return self.number_of_recovered() / float(self.G.number_of_nodes())

    def get_outbreak_size(self):
        """return the size of the outbreak (number of susceptibles and recovered)"""

        return self.G.number_of_nodes() - self.number_of_susceptibles()

    def _get_max_t(self):
        """return the time of the last event"""

        return max([
            self.s_of_t[-1][0],
            self.i_of_t[-1][0],
            self.r_of_t[-1][0],
        ])

    def _get_x_of_t(self, arr):
        """get the time of the last event, append it to the list and pass back an nd.array"""

        t_max = self._get_max_t()
        arr = list(arr)

        if arr[-1][0] < t_max:
            arr.append([t_max, arr[-1][1]])

        arr = np.array(arr)
        return arr[:, 1], arr[:, 0]

    def get_s_of_t(self):
        return self._get_x_of_t(self.s_of_t)

    def get_i_of_t(self):
        return self._get_x_of_t(self.i_of_t)

    def get_r_of_t(self):
        return self._get_x_of_t(self.r_of_t)




