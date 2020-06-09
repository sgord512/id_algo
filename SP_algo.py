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

class Graph(IGraph):

    # These are all helper methods to make constructing graphs easier
    @staticmethod
    def _double_and_flip_list(edge_list):
        return [(v,u) if flip else (u,v) for (u,v) in edge_list for flip in range(2)]

    @staticmethod
    def _edge_dict_to_list(ed):
        return [(k, v) for k, ls in ed.items() for v in ls]

    @staticmethod
    def _values_to_lists_in_edge_dict(ed):
        return {k: (v if type(v) == list else [v]) for k, v in ed.items()}

    def _edge_list_from_dict(ed):
        return Graph._edge_dict_to_list(Graph._values_to_lists_in_edge_dict(ed))

    def get_num_nodes(self):
        return len(self.vs)

    # This is the main method used to easily construct a new graph
    @staticmethod
    def FromDicts(vertices, observed_edge_dict=None, hidden_edge_dict=None):
        observed_edges = Graph._edge_list_from_dict(observed_edge_dict) if observed_edge_dict else []
        hidden_edges = Graph._edge_list_from_dict(hidden_edge_dict) if hidden_edge_dict else []
        g = Graph()
        g.add_vertices(vertices)
        g.add_observed_edges(observed_edges)
        g.add_hidden_edges(hidden_edges)
        return g

    def __init__(self):
        super().__init__(directed=True)

    def add_edges_from(self, source, targets):
        self._add_edges_base([(source, v) for v in targets], observed=True)

    def _add_edges_base(self, edge_list, **kwds):
        for (u,v) in edge_list:
            super().add_edge(u, v, **kwds)

    def add_observed_edge(self, u, v):
        self.add_observed_edges([(u,v)])

    def add_observed_edges(self, edge_list):
        self._add_edges_base(edge_list, observed=True)

    def add_hidden_edges(self, edge_list):
        full_edge_list = Graph._double_and_flip_list(edge_list)
        self._add_edges_base(full_edge_list, observed=False)

    def __str__(self):
        return super().__str__()

    def __repr__(self):
        return super().__repr__()

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
        return ",".join(variables)

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

# Get connected components of a graph.
def getComponents(G):
    components = G.components(mode=igraph.STRONG)
    return [set(G.vs[c]["name"]) for c in components]

def getAncestors(G, y):
    # This returns the ancestors of the vertices in the set y
    n = G.get_num_nodes()
    ancestorIndices = utilities.flatten(G.obs.neighborhood(y, order=n, mode=igraph.IN))
    return set(G.obs.vs[ancestorIndices]["name"])

def constructObservedAndConfounded(G):
    G.obs = G.subgraph_edges(G.es.select(observed_eq=True), delete_vertices=False)
    G.confounded = G.subgraph_edges(G.es.select(observed_eq=False), delete_vertices=False)

def prettyPrint(formatObj, *args, depth=0, verbose):
    if not verbose:
        return
    if len(args) == 0:
        print(("\t" * depth) + str(formatObj))
    else:
        print(("\t" * depth) + formatObj.format(*args))

def causalEffectStr(x, y):
    return "P_{{{:}}}({:})".format(",".join(x), ",".join(y))

def conditionalProbStr(y, x):
    return "P({:}|{:})".format(",".join(y), ",".join(x))

# We'll use the name attribute of vertices to label vertices.
# We'll indicate that an edge is unobserved by an attribute observed = False
def ID(y, x, G, P=None, topordering=None, depth=0, verbose=False):
    try:
        if P == None:
            P = DistributionExpression.Atomic(
                all_variables = G.vs["name"],
                output_set=G.vs["name"],
                trivial=True
            )
        y = [y] if type(y) == str else y # If I just have a single vertex, I need to wrap it in a list.
        x = [x] if type(x) == str else x
        prettyPrint("Computing " + causalEffectStr(x,y) + " from " + str(P), depth=depth, verbose=verbose)
        constructObservedAndConfounded(G)
        y = set(y)
        x = set(x)
        xindices = G.vs.select(name_in=x)
        v = set(G.vs["name"])
        n = G.get_num_nodes()
        anc = getAncestors(G, y)
        if topordering == None:
            topordering = G.vs[G.obs.topological_sorting(mode=igraph.OUT)]["name"]
            prettyPrint("Topological ordering: {:}", topordering, depth=depth, verbose=verbose)

        # LINE 1
        if len(set(x)) == 0:
            finalP = DistributionExpression.Atomic(
                all_variables=v,
                output_set=y,
                marginalized_set=v-y,
                conditioning_set=None,
                context=P)
            prettyPrint("LINE 1", depth=depth, verbose=verbose)
            prettyPrint(finalP, depth=depth, verbose=verbose)
            #print()
            return finalP

        # LINE 2
        nonanc = v - anc
        #prettyPrint("Nonancestors of {:}: {:}", y, nonanc, depth=depth, verbose=verbose)
        if len(nonanc) > 0:
            for vi in nonanc:
                topordering.remove(vi)
            prettyPrint("LINE 2: non-ancestors are {:}", nonanc, depth=depth, verbose=verbose)
            newP = DistributionExpression.Atomic(
                all_variables=anc,
                output_set=anc,
                marginalized_set=nonanc,
                conditioning_set=None,
                context=P)
            output = ID(y,x & anc, G.induced_subgraph(anc), newP, topordering, depth=depth+1, verbose=verbose)
            prettyPrint(output, depth=depth, verbose=verbose)
            #print()
            return output

        # LINE 3
        xindices = G.vs.select(name_in=x)
        G.postIntervention = G.subgraph_edges(G.es.select(_to_notin=xindices), delete_vertices=False)
        constructObservedAndConfounded(G.postIntervention)
        ancestorsAvoidingX = getAncestors(G.postIntervention, y)
        w = v - x - ancestorsAvoidingX
        if len(w) > 0:
            prettyPrint("LINE 3", depth=depth, verbose=verbose)
            output = ID(y, x | w, G, P, depth=depth+1, verbose=verbose)
            prettyPrint(output, depth=depth, verbose=verbose)
            #print()
            return output

        # LINE 4
        # I want to get all the C-components of G[V\X].
        G.withoutX = G.induced_subgraph(v - x)
        components = getComponents(G.confounded.induced_subgraph(v-x))
        numComponents = len(components)
        if numComponents > 1:
            prettyPrint("LINE 4", depth=depth, verbose=verbose)
            marginalized_set = v - (y | x)
            children = []
            for si in components:
                prettyPrint("Executing ID({:}, {:}, G, P)", si, v-si, depth=depth, verbose=verbose)
                children.append(ID(si, v - si, G, P, topordering[:], depth=depth+1, verbose=verbose))
            output = DistributionExpression.Recursive(
                all_variables=v,
                children=children,
                marginalized_set=marginalized_set)
            prettyPrint(output, depth=depth, verbose=verbose)
            return output

        # LINE 5
        # By the time we get here, we know that G[V\X] has a single C-component
        s = components[0]
        # prettyPrint("G[V\X] has only a single C-component: {{{:}}}".format(",".join(s)), depth=depth, verbose=verbose)
        componentsG = getComponents(G.confounded)
        # If G has only a single C-component, then the effect is unidentifiable.
        if len(componentsG) == 1:
            prettyPrint("LINE 5", depth=depth, verbose=verbose)
            raise UnidentifiableEffect(x, y, componentsG[0], s)

        sprime = None
        for si in componentsG:
            if s <= si:
                sprime = si
                break

        # LINE 6
        # This is the case where S is a C-component of G itself.
        if len(s ^ sprime) == 0:
            marginalized_set = s - y
            children = []
            for vi in s:
                ix = topordering.index(vi)
                predecessorsi = topordering[0:ix]
                # prettyPrint(conditionalProbStr([vi], predecessorsi), depth=depth, verbose=verbose)
                children.append(DistributionExpression.Atomic(
                    all_variables=v,
                    marginalized_set=v-set([vi])-set(predecessorsi),
                    output_set=[vi],
                    conditioning_set=predecessorsi,
                    context=P)
                )
            prettyPrint("LINE 6", depth=depth, verbose=verbose)
            output = DistributionExpression.Recursive(
                all_variables=v,
                children=children,
                marginalized_set=marginalized_set)
            prettyPrint(output, depth=depth, verbose=verbose)
            #print()
            return output

        # LINE 7
        # This is the case where S is a subset of larger C-component S'
        children = []
        for vi in sprime:
            ix = topordering.index(vi)
            predecessorsi = topordering[0:ix]
            children.append(DistributionExpression.Atomic(
                all_variables=v,
                output_set=[vi],
                marginalized_set=v-set([vi])-set(predecessorsi),
                conditioning_set=predecessorsi,
                context=P)
            )
        G.sprime = G.induced_subgraph(sprime)
        newP = DistributionExpression.Recursive(
            all_variables=G.sprime.vs["name"],
            children=children)
        prettyPrint("LINE 7", depth=depth, verbose=verbose)
        output = ID(y, x & sprime, G.sprime, newP, depth=depth+1, verbose=verbose)
        prettyPrint(output, depth=depth, verbose=verbose)
        #print()
        return output
    except UnidentifiableEffect as e:
        if depth == 0:
            return e
        else:
            raise e

A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,W,X,Y,Z = "a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,w,x,y,z".split(',')
X1,X2,X3,X4,X5 = "x_1,x_2,x_3,x_4,x_5".split(',')
Y1,Y2,Y3,Y4,Y5 = "y_1,y_2,y_3,y_4,y_5".split(',')
Z1,Z2,Z3,Z4,Z5 = "z_1,z_2,z_3,z_4,z_5".split(',')
W1,W2,W3,W4,W5 = "w_1,w_2,w_3,w_4,w_5".split(',')

if __name__ == "__main__":
    SingleEdge = Graph.FromDicts([X, Y], {X: Y}, {})
    SmokingExample = Graph.FromDicts([X, Y, W], {X: Y, W: X}, {W: Y})

    singleEdgeOut = ID(Y, X, SingleEdge)
    smokingOut = ID(Y, X, SmokingExample)

    PaperFigure1 = Graph.FromDicts(
        [X, Y, W, Z],
        {W: [X, Z],
         X: Z,
         Z: Y},
        {X: Y}
    )

    paperFigure1Out = ID(Y, X, PaperFigure1)

    G = Graph.FromDicts(
        [A, W, X, Y, B],
        {A: Y,
         X: Y,
         W: [X, B],
         B: Y},
        {A: [W, Y],
         W: B,
         B: X})
    gOut = ID(Y, X, G)
    G = Graph()
    G.add_vertices([X,Y])
    G.add_observed_edge(X,Y)

    G2 = Graph()

    G2.add_vertices([X,Y,W,Z])
    G2.add_observed_edges([
        (W,X),
        (W,Z),
        (Z,Y),
        (X,Z)
    ])
    G2.add_hidden_edges([(X,Y)])
    expr2 = ID(Y,X,G2)

    G3 = Graph()
    G3.add_vertices([X,Y,W,Z])
    G3.add_observed_edges([
        (X,Y),
        (Z,Y),
        (W,X),
        (Z,X)
    ])
    expr3 = ID(Y,X,G3)

    G4 = Graph()
    G4.add_vertices([X,Y,W,Z])
    G4.add_observed_edges([
        (X,Y),
        (Z,Y),
        (W,X),
        (Z,X)
    ])
    G4.add_hidden_edges([(W,Y)])
    expr4 = ID(Y,X,G4)

    G5 = Graph()
    G5.add_vertices([X,Y,W,Z,A])
    G5.add_observed_edges([
        (X,A),
        (A,Y),
        (Z,Y),
        (W,X),
        (Z,X)
    ])
    G5.add_hidden_edges([(Z,W)])
    expr5 = ID(Y,X,G5)

    G6 = Graph()
    G6.add_vertices([X1,Y1,X2,Y2])
    G6.add_observed_edges([
        (X1,Y1),
        (X2,Y2),
        (Y1,X2)
    ])
    G6.add_observed_edges([(Y1,Y2)])
    expr6 = ID([Y1,Y2],[X1,X2], G6)

    G7 = Graph()
    G7.add_vertices([X1, X2, X3, Z1, Z2, W1, W2, Y])
    G7.add_observed_edges([
        (X1, Z1),
        (Z1, W1),
        (W1, X2),
        (X2, Z2),
        (Z2, W2),
        (W2, X3),
        (Z1, Y),
        (Z2, Y),
        (X3, Y)
    ])
    G7.add_hidden_edges([
        (Z1, X2),
        (Z1, Z2),
        (X1, W1),
        (W1, W2),
        (W1, Z2),
    ])

    expr7 = ID(Y, [X1, X2, X3], G7)

    print(expr2)
    print(expr3)
    print(expr4)
    print(expr5)
    print(expr6)
    print(expr7)
