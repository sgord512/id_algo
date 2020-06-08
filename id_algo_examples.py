from id_algo_implementation import *

catalog = dict()
G1 = Graph()

class Entry():
    def __init__(self, G, Y, X):
        self.G = G
        self.Y = Y
        self.X = X

G = Graph.FromDicts([Z, W, X1, X2, Y],
    {Z: [Y, W], X1: [W], W: [X2], X2: [Y]},
    {X1: W, X2: W}
)

#ID(Y, [X1, X2], G)

G2 = Graph.FromDicts([W, X1, X2, Y],
    {X1: W, W: X2, X2: Y},
    {X1: W, W: Y}
)

G3 = Graph.FromDicts([X1, W, X2, Y],
    {X1: W, W: X2, X2: Y},
    {X1: W})

G4 = Graph.FromDicts([Z, W, X1, X2, Y],
    {Z: Y, X1: W, W: X2, X2: Y},
    {X1: W, X2: W, Z: W}
)

G5 = Graph.FromDicts([A, B, W, X1, X2, Y],
                     {A: Y, W: X2, B: Y, X1: Y, X2: Y},
                     {A: [Y,W], W: B, B: X1})

G6 = Graph.FromDicts([A, B, X1, X2, Y],
                     {A: Y, X1: Y, X2: Y, B: Y},
                     {A: [Y, B], B: X1})
