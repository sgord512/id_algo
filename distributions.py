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

class DistributionExpression():
    # TODO: Figure out when I can have multiple summations and make sure that I always merge them!
    # TODO: Use d-separation to clean up complicated expressions (particularly those resulting from lines 6 or 7). 

    # All distribution expressions are defined over a set of variables,
    #    self.variable_set, which represents the entire visible vertex set of the graph over which the expr is defined.
    # There are two types of distribution expressions, atomic distributions and recursive distributions.
    # Atomic DistributionExpressions correspond to expressions of the form:
    #    \sum_{marginalized_set} P'(variables|conditioning_set)
    # where P'(variables|conditioning_set) may itself be of the form:
    #    \frac{\sum_{all_variables-variables-conditioning_set} P"(all_variables)}
    #         {\sum_{all_variables-conditioning_set) P"(all_variables)}

    # Recursive DistributionExpressions correspond to expressions of the form:
    #    \sum_{marginalized_set} \prod_{i=1}^k P_i(...)
    # where each of the P_i terms are distribution expressions.

    # The two types of distribution expressions can be distinguished as follows:
    # Recursive distexprs have children and no conditioning set or output set, and atomic distexprs have an (optional) conditioning set, an output set and no children.

    # For atomic distributions, all_variables = marginalized_set + output_set + conditioning_set.
    # For recursive distributions. 

    # TODO: FINISH THIS
    def __init__(self, all_variables, output_set, conditioning_set, children, marginalized_set, context, trivial):
        self.all_variables = set(all_variables)
        self.output_set = output_set and set(output_set)
        self.conditioning_set = conditioning_set and set(conditioning_set)
        self.children = children
        self.marginalized_set = set(marginalized_set) if marginalized_set is not None else set()
        self.context = context
        self.trivial = trivial
        assert(not self.children or
               (not self.conditioning_set and not self.output_set and not self.context))
        assert(self.output_set or self.children)

    def has_nontrivial_context(self):
        return self.context and not self.context.trivial

    def is_atomic(self):
        return self.children == None

    def is_fraction(self):
        return (self.is_atomic() and
                self.has_nontrivial_context() and
                self.conditioning_set is not None)

    @staticmethod
    def summation_over_variables(output_set):
        return "\\sum_{{{:}}}".format(",".join(output_set))

    @staticmethod
    def parenthesize(distexpr, depthmap, always_parenthesize=False):
        subexpr = distexpr.to_string(depthmap)
        if (distexpr.is_fraction() or
            (not distexpr.is_atomic() and len(distexpr.children) > 1) or
            always_parenthesize):
            return "\\left[" + subexpr + "\\right]"
        return subexpr

    @staticmethod
    def augmented_depthmap(depthmap, output_set):
        dm = depthmap.copy()
        for v in output_set:
            dm[v] += 1
        return dm

    @staticmethod
    def apply_depthmap(output_set, depthmap):
        return [v + ("'" * depthmap[v]) for v in output_set]

    @staticmethod
    def variable_set_expr(variables):
        # We'll assume a variable is a single character followed by any number of subscripts.

        return ",".join(map(DistributionExpression.texify_variable,variables))

    def texify_variable(var):
        ch = var[0]
        if len(var) > 1:
            return ch + '_{{{:}}}'.format(var[1:])
        else:
            return ch

    @staticmethod
    def fraction_expr(numerator_expr, denominator_expr):
        return "\\left(\\frac{{{:}}}{{{:}}}\\right)".format(numerator_expr, denominator_expr)

    @staticmethod
    def concatenate_exprs(exprs):
        return "".join(exprs)

    @staticmethod
    def Recursive(all_variables, children, marginalized_set=None):
        return DistributionExpression(all_variables,
                                      output_set=None,
                                      conditioning_set=None,
                                      children=children,
                                      marginalized_set=marginalized_set,
                                      context=None,
                                      trivial=False)

    @staticmethod
    def Atomic(all_variables, output_set, conditioning_set=None, marginalized_set=None, context=None, trivial=False):
        return DistributionExpression(all_variables=all_variables,
                                      output_set=output_set,
                                      conditioning_set=conditioning_set,
                                      children=None,
                                      marginalized_set=marginalized_set,
                                      context=context,
                                      trivial=trivial)

    def copy(self):
        return DistributionExpression(
            all_variables=self.all_variables,
            output_set=self.output_set,
            conditioning_set=self.conditioning_set,
            children=self.children,
            marginalized_set=self.marginalized_set,
            context=self.context,
            trivial=self.trivial)

    def marginalize(self, new_marginalized_set):
        marginalized_set = self.marginalized_set or set()
        if new_marginalized_set & marginalized_set:
            raise Exception("Repeated marginalization!")
        self.marginalized_set = new_marginalized_set | marginalized_set

    def simplify(self):
        if self.context:
            self.context = self.context.simplify()
        elif self.children:
            self.children = [child.simplify for child in self.children]

        if self.is_atomic():
            pass

        # TODO: Finish implementing this method.

    def numerator_expr(self):
        assert(self.is_atomic() and self.has_nontrivial_context())
        # This is the case where we have a non-trivial context.
        all_variables = self.all_variables
        output_set = self.output_set
        conditioning_set = self.conditioning_set or set()
        numerator_marginalized_set = all_variables - output_set - conditioning_set
        numerator_distexpr = self.context.copy()
        numerator_distexpr.marginalize(numerator_marginalized_set)
        return numerator_distexpr

    def denominator_expr(self):
        all_variables = self.all_variables
        conditioning_set = self.conditioning_set or set()
        denominator_marginalized_set = all_variables - conditioning_set
        denominator_distexpr = self.context.copy()
        denominator_distexpr.marginalize(denominator_marginalized_set)
        return denominator_distexpr

    def to_string(self, depthmap=None):
        if depthmap is None:
            depthmap = defaultdict(int)
        if len(self.marginalized_set) > 0:
            marginalized_set = sorted(self.marginalized_set)
            depthmap = DistributionExpression.augmented_depthmap(depthmap, marginalized_set)
            marked_marginalized_set = self.apply_depthmap(marginalized_set, depthmap)
            summation_expr = DistributionExpression.summation_over_variables(marked_marginalized_set)
        else:
            summation_expr = ""
        if self.is_atomic():
            output_set = self.output_set
            conditioning_set = self.conditioning_set or set()
            if self.has_nontrivial_context():
                # This is the case where we have a non-trivial context.
                numerator_distexpr = self.numerator_expr()
                numerator_expr = numerator_distexpr.to_string(depthmap.copy())

                if not conditioning_set:
                    # If the conditioning set is empty, then we can just return the expression computed thus far.
                    return summation_expr + numerator_expr
                else:
                    # If the conditioning set is non-empty, then we need to have a denominator as well.
                    denominator_distexpr = self.denominator_expr()
                    denominator_expr = denominator_distexpr.to_string(depthmap.copy())
                    return summation_expr + DistributionExpression.fraction_expr(numerator_expr, denominator_expr)
            else:
                # When we have either a trivial context or no context, we can just construct expressions in the obvious way.
                output_set = self.output_set
                marked_output_set = DistributionExpression.apply_depthmap(output_set, depthmap)
                output_set_expr = DistributionExpression.variable_set_expr(marked_output_set)
                if conditioning_set:
                    marked_conditioning_set = DistributionExpression.apply_depthmap(conditioning_set, depthmap)
                    conditioning_set_expr = "|" + DistributionExpression.variable_set_expr(marked_conditioning_set)
                else:
                    conditioning_set_expr = ""
                return "P({:}{:})".format(output_set_expr,conditioning_set_expr)
        else:
            if len(self.children) > 1:
                children_exprs = [DistributionExpression.parenthesize(child, depthmap.copy()) for child in self.children]
            else:
                children_exprs = [child.to_string(depthmap.copy()) for child in self.children]
            return summation_expr + DistributionExpression.concatenate_exprs(children_exprs)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.__str__()

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

