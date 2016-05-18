#!/usr/env/bin python3
# -*- coding: utf-8 -*-

from SIR import *
from Node import Node


class SI(SIR):

    def __init__(self,
                 G,
                 infection_rate,
                 infection_seeds=1,
                 verbose=False
                 ):
        SIR.__init__(
                self,
                G,
                infection_rate,
                recovery_rate=0.,
                infection_seeds=infection_seeds,
                verbose=verbose
        )
