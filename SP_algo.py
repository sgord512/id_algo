#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:43:49 2019

@author: spencergordon
"""

# This should work with Python 3.7.4 and igraph version 0.7.1
from graph import Graphs
from collections import defaultdict
from distributions import DistributionExpression, UnidentifiableEffect
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
        anc = G.getAncestors(y)
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
        G.postIntervention = G.construct_post_intervention_subgraph(x)
        G.postIntervention.constructObservedAndConfounded()
        ancestorsAvoidingX = G.postIntervention.getAncestors(y)
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
        components = G.confounded.induced_subgraph(v-x).getComponents()
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
