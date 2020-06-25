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
import re

class Probability():
    def __init__(self,
                 var = [],
                 cond = [],
                 sumset = [],
                 is_product=False,
                 children=[],
                 is_fraction=False,
                 den=[],
                 num=[],
                 domain=0):

        # All of the collection variables should be lists and not sets!
        # This will help when we eventually print an expression.
        # To do this, I need to ensure that I always topological sort everything along the way!
        self.var = var
        self.cond = cond
        self.sumset = sumset
        self.is_product = is_product
        self.children = children
        self.is_fraction = is_fraction
        self.den = den
        self.num = num
        # This will be filled in by the ID algo and used to initialize variable contexts correctly.
        self.query = None

    def copy(self):
        return Probability(
            self.var[:],
            self.cond[:],
            self.sumset[:],
            self.is_product,
            self.children[:],
            self.is_fraction,
            self.den[:],
            self.num[:]
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
            P_new.is_fraction = True
            P_new.num = P_num
            P_new.den = P_den
        else:
            P_new = P_num
        return P_new

    def __setattr__(self, name, value):
        assert(type(value) is not set)
        super().__setattr__(name, value)

    # `variables` will be a list of strings,
    # `context` will be a mapping from strings->number represented by a dict.
    # If increase_depth is True then I increment depth for both new and old vars. If not, I just add new vars with depth 0.
    @staticmethod
    def update_bound_depth(variables, context=None):
        if context is None:
            context = defaultdict(int)
        else:
            # Clone the dict to avoid overwriting the values the caller may need.
            context = context.copy()
        for v in variables:
                context[v] += 1
        return context

    @staticmethod
    def generate_superscript_from_bound_depth(context):
        superscript = defaultdict(str)
        for (var, count) in context.items():
            if count > 1:
                superscript[var] = "^{{{:}}}".format("".join(["\\prime" for _ in range(count - 1)]))
            else:
                superscript[var] = ""
        return superscript

    @staticmethod
    def apply_superscript_to_vars(variables, superscript):
        return [var + superscript[var] for var in variables]

    @staticmethod
    def _format_vars(variables, superscript=None):
        pat = re.compile('[0-9]+$')
        formatted_vars = []
        for var in variables:
            sup = superscript[var]
            match = pat.search(var)
            if match:
                ix = match.span()[0]
                var_prefix, var_suffix = var[0:ix], var[ix:]
                var = "{:}_{{{:}}}".format(var_prefix, var_suffix)
            formatted_vars.append(var + sup)
        return formatted_vars

    def get_expression(self, mark_bound_vars=False):
        query_vars = self.query['X'] | self.query['Y']
        initial_context = Probability.update_bound_depth(query_vars)
        return self._get_expression(mark_bound_vars, initial_context, parenthesize_expression=False)

    def _get_expression(self, mark_bound_vars, context, parenthesize_expression):
        P = ""
        has_summation_vars = len(self.sumset) > 0
        summation_string = ""
        variable_string = ""
        conditioning_string = ""
        superscript = None
        if has_summation_vars:
            if mark_bound_vars:
                context = Probability.update_bound_depth(self.sumset, context)
                superscript = Probability.generate_superscript_from_bound_depth(context)
            summation_string = ",".join(Probability._format_vars(self.sumset, superscript))

            if parenthesize_expression:
                P = P + "\\left(\\sum_{{{:}}}".format(summation_string)
            else:
                P = P + "\\sum_{{{:}}}".format(summation_string)

        if self.is_fraction:
            numerator_string = self._get_expression(mark_bound_vars, context, parenthesize_expression)
            denominator_string = self._get_expression(mark_bound_vars, context, parenthesize_expression=True)
            P = P + "\\frac{{{:}}}{{{:}}}".format(numerator_expression, denominator_expression)

        if self.is_product:
            for child in self.children:
                P = P + child._get_expression(mark_bound_vars, context, parenthesize_expression=True)

        if not (self.is_fraction or self.is_product):
            P = P + "P("
            if mark_bound_vars:
                superscript = Probability.generate_superscript_from_bound_depth(context)
            var_string = ",".join(Probability._format_vars(self.var, superscript))

            P = P + var_string

            if len(self.cond) > 0:
                if mark_bound_vars:
                    superscript = Probability.generate_superscript_from_bound_depth(context)
                cond_string = ",".join(Probability._format_vars(self.cond, superscript))
                cond_string = "|" + cond_string + ")"
            else:
                cond_string = ")"
            P = P + cond_string
        if has_summation_vars and parenthesize_expression:
            P = P + "\\right)"
        return P

    def __str__(self):
        return self.get_expression(mark_bound_vars=True)
