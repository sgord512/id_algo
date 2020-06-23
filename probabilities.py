#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:43:49 2019

@author: spencergordon
"""

# This should work with Python 3.7.4 and igraph version 0.7.1
import igraph
from igraph import Graph as IGraph
from collections import defaultdict
import utilities

class Probability():
    def __init__(self,
                 var = [],
                 cond = [],
                 sumset = [],
                 do =[],
                 product=False,
                 children=[],
                 fraction=False,
                 den=[],
                 num=[],
                 domain=0,
                 is_sum=False,
                 weight=1):

        # All of the collection variables should be lists and not sets!
        # This will help when we eventually print an expression.
        # To do this, I need to ensure that I always topological sort everything along the way!
        self.var = var
        self.cond = cond
        self.sumset = sumset
        self.do = do
        self.product = product
        self.children = children
        self.fraction = fraction
        self.den = den
        self.num = num
        self.domain = domain
        self.is_sum = is_sum
        self.weight = weight

    def copy(self):
        return Probability(
            self.var[:],
            self.cond[:],
            self.sumset[:],
            self.do[:],
            self.product,
            self.children[:],
            self.fraction,
            self.den[:],
            self.num[:],
            self.domain,
            self.is_sum,
            self.weight
        )

    def parse_joint(self, v, cond, var, top_ord):
        P_new = Probability()
        P_num = self
        num_sumset = set(self.num_sumset) | (var - (v | cond))
        P_num.sumset = utilities.sort_subset_by_list(sumset, top_ord)
        if len(cond) > 0:
            P_den = P
            den_sumset = set(self.sumset) | (var - cond)
            P_den.sumset = utilities.sort_subset_by_list(den_sumset, top_ord)
            P_new.fraction = True
            P_new.num = P_num
            P_new.den = P_den
        else:
            P_new = P_num
        return P_new

    def __setattr__(self, name, value):
        assert(type(value) is not set)
        super().__setattr__(name, value)

class UnidentifiableEffect(Exception):
    def __init__(self, x, y, G, S):
        self.x = x
        self.y = y
        self.G = G
        self.S = S

    def __str__(self):
        output = "UnidentifiableEffect: {:}. Hedge is\n\tF: {:}\n\tF': {:}".format(
            causalEffectStr(self.x, self.y), self.G, self.S
        )
        return output
