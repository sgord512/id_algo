#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:43:49 2019

@author: spencergordon
"""

# This should work with Python 3.7.4 and igraph version 0.7.1
from addict import Dict
from graph import Graph
from probabilities import Probability, UnidentifiableEffect
import utilities

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

class Call(Dict):
    def __init__(self, y = None, x = None, P = None, G = None, line = None):
        super().__init__()
        self.y = y
        self.x = x
        self.P = P
        self.G = G
        self.line = line

class Tree():
    def __init__(self):
        self.call = Call()
        self.root = None
        self.branch = None

class IDOutput():
    def __init__(self, P, tree):
        self.P = P
        self.tree = tree

def ID(y, x, G, verbose=False):
    G._construct_observed_and_confounded()
    P = Probability()
    y = [y] if type(y) == str else y # If I just have a single vertex, I need to wrap it in a list.
    x = [x] if type(x) == str else x
    y = set(y)
    x = set(x)
    tree = Tree()
    top_ord = G.get_topological_order()
    prettyPrint("Topological order: {:}", top_ord, depth=0, verbose=verbose)
    return ID_internal(y, x, G, Probability(), top_ord, tree=tree, depth=0, verbose=verbose)

def ID_internal(y, x, G, P, top_ord, tree, depth, verbose):
    try:
        prettyPrint("Computing " + causalEffectStr(x,y) + " from " + str(P), depth=depth, verbose=verbose)

        # The following variables will be useful throughout the algorithm
        v = set(G.get_vertices())
        n = G.get_num_nodes()
        ancestors = set(G.get_ancestors(y))

        # Record the basic information about this call to the ID algorithm
        if len(P.var) == 0 and not (P.product or P.fraction):
            tree.call = Call(y, x, Probability(var=list(v)), G, line=None)
        else:
            tree.call = Call(y, x, P, G, line=None)

        # The topological order we were given may contain variables that don't exist in the subgraph we're working with.
        # We'll filter it to make sure it fits our subgraph.
        # We'll still keep the original around to sort sets we deal with throughout.
        new_top_ord = [v_i for v_i in top_ord[:] if v_i in v]

        ########## LINE 1 ##########
        if len(set(x)) == 0:
            if P.product or P.fraction:
                P_sumset = (v - y) | set(P.sumset)
                P.sumset = utilities.sort_subset_by_list(P_sumset, top_ord)
            else:
                P.var = utilities.sort_subset_by_list(y, top_ord)
            tree.call.line = 1
            tree.root = P
            return IDOutput(P, tree)

        ########## LINE 2 ##########
        nonancestors = v - ancestors
        #prettyPrint("Nonancestors of {:}: {:}", y, nonanc, depth=depth, verbose=verbose)
        if len(nonancestors) > 0:
            #prettyPrint("LINE 2: non-ancestors are {:}", nonanc, depth=depth, verbose=verbose)
            G_ancestors = G.induced_subgraph(ancestors)
            if (P.product or P.fraction):
                P_sumset = (v - ancestors) | set(P.sumset)
                P.sumset = utilities.sort_subset_by_list(P_sumset, top_ord)
            else:
                P.var = utilities.sort_subset_by_list(ancestors, top_ord)
            out = ID_internal(y, x & ancestors, G_ancestors, P, top_ord, tree=Tree(), depth=depth+1, verbose=verbose)
            tree.branch = [out.tree]
            tree.call.line = 2
            tree.call.anc = ancestors
            prettyPrint(out.P, depth=depth, verbose=verbose)
            return IDOutput(P=out.P, tree=tree)

        ########## LINE 3 ##########
        G.post_intervention = G.construct_post_intervention_subgraph(x)
        ancestors_avoiding_x = G.post_intervention.get_ancestors(y)
        w = (v - x) - ancestors_avoiding_x
        if len(w) > 0:
            prettyPrint("LINE 3", depth=depth, verbose=verbose)
            w_connected = w & G.post_intervention.get_all_connected_vertices(y)

            if len(w_connected) < len(w):
                v_new = v - (w - w_connected)
                G = G.induced_subgraph(v_new)

            out = ID_internal(y, x | w_connected, G, P, top_ord, Tree(), depth=depth+1, verbose=verbose)
            tree.branch = [out.tree]
            tree.call.line = 3
            tree.call.w = w
            tree.call.anc_post_intervention = ancestors_avoiding_x
            prettyPrint(out.P, depth=depth, verbose=verbose)
            return IDOutput(P=out.P, tree=tree)

        # I want to get all the C-components of G[V\X].
        G.withoutX = G.induced_subgraph(v - x)
        components = G.withoutX.get_c_components()
        num_components = len(components)

        ########### LINE 4 ##########
        if num_components > 1:
            tree.call.line = 4
            prettyPrint("LINE 4", depth=depth, verbose=verbose)
            recursive_calls = []
            for s_i in components:
                prettyPrint("Executing ID_internal({:}, {:}, G, P)", s_i, v-s_i, depth=depth, verbose=verbose)
                recursive_calls.append(ID_internal(s_i, v - s_i, G, P, top_ord[:], Tree(), depth=depth+1, verbose=verbose))
            tree.branch = [child.tree for child in recursive_calls]
            sumset = v - (y | x)
            P = Probability(
                sumset = utilities.sort_subset_by_list(sumset, top_ord),
                product = True,
                children = [call.P for call in recursive_calls])
            prettyPrint(P, depth=depth, verbose=verbose)
            return IDOutput(P=P, tree=tree)

        ########### LINE 5 ##########
        # By the time we get here, we know that G[V\X] has a single C-component
        s = components[0]
        # prettyPrint("G[V\X] has only a single C-component: {{{:}}}".format(",".join(s)), depth=depth, verbose=verbose)
        components_G = G.get_c_components()
        # If G has only a single C-component, then the effect is unidentifiable.
        if len(components_G) == 1:
            prettyPrint("LINE 5", depth=depth, verbose=verbose)
            raise UnidentifiableEffect(x, y, components_G[0], s)

        # Here we figure out which c-component of G is a superset of the c-component of G[V\X]
        s_prime = None
        for s_i in components_G:
            if s <= s_i:
                s_prime = s_i
                break

        ########### LINE 6 ##########
        # This is the case where S is a C-component of G itself.
        if len(s ^ s_prime) == 0:
            prettyPrint("LINE 6", depth=depth, verbose=verbose)
            tree.call.line = 6
            tree.call.s = s
            product_list = [None for _ in range(len(s))]
            P_prod = Probability()
            s = utilities.sort_subset_by_list(s, top_ord)
            for (i, v_i) in enumerate(s):
                ix = new_top_ord.index(v_i)
                cond_set = new_top_ord[0:ix]
                if P.product:
                    P_prod = P.parse_joint(set([v_i]), set(cond_set), v, top_ord)
                else:
                    P_prod = P.copy()
                    P_prod.var = [v_i]
                    P_prod.cond = cond_set

                product_list[len(s) - i - 1] = P_prod

            if len(s) > 1:
                P_new_sumset = set(s) - y
                P_new = Probability(
                    sumset = utilities.sort_subset_by_list(P_new_sumset, top_ord),
                    product = True,
                    children = product_list
                )
                tree.root = P_new
                prettyPrint(P_new, depth=depth, verbose=verbose)
                return IDOutput(P=P_new, tree=tree)

            if P_prod.product or P_prod.fraction:
                P_prod_sumset = set(P_prod.sumset) | (set(s) - y)
                P_prod.sumset = utilities.sort_subset_by_list(P_prod_sumset, top_ord)
            else:
                P_prod_var = set(P_prod.var) - (set(P_prod.sumset) | (set(s) - y))
                P_prod.var = utilities.sort_subset_by_list(P_prod_var, top_ord)
            tree.root = P_prod_var
            prettyPrint(P_prod, depth=depth, verbose=verbose)
            return IDOutput(P=P_prod, tree=tree)

        ########## LINE 7 ##########
        # This is the case where S is a subset of larger C-component S'
        tree.call.s = s
        tree.call.line = 7
        # Note that we're going to set s <- s_prime to match the code from Tikka and Karvanen.
        # I think it's confusing, but I'll try it this way first. 
        s = utilities.sort_subset_by_list(s_prime, top_ord)
        tree.call.s_prime = s

        G.s = G.induced_subgraph(s)
        product_list = [None for _ in range(len(s))]
        for (i, v_i) in enumerate(s):
            ix = new_top_ord.index(v_i)
            cond_set = new_top_ord[0:ix]
            P_prod = P.copy()
            P_prod.var = [v_i]
            P_prod.cond = cond_set
            product_list[len(s) - i - 1] = P_prod

        x_new = set(s) & x
        out = None
        if len(s) > 1:
            P_recursive = Probability(product=True, children=product_list)
            out = ID_internal(y, x_new, G.s, P_recursive, top_ord, Tree(), depth=depth+1, verbose=verbose)
        else:
            out = ID_internal(y, x_new, G.s, product_list[0], top_ord, Tree(), depth=depth+1, verbose=verbose)
        tree.branch = [out.tree]
        prettyPrint(out.P, depth=depth, verbose=verbose)
        return IDOutput(P=out.P, tree=tree)

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
