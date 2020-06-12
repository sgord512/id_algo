from SP_algo import *

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

G7 = Graph.FromDicts([W1, X1, Y1, Z1, W2, X2, Y2, Z2, W3, X3, Y3, Z3, W4, X4, Y4, Z4],
                     {W1: [W2, Z2],
                      X1: [X2, W2],
                      Y1: [Y2, X2],
                      Z1: [Z2, Y2],
                      W2: [W3, Z3],
                      X2: [X3, W3],
                      Y2: [Y3, X3],
                      Z2: [Z3, Y3],
                      W3: [W4, Z4],
                      X3: [X4, W4],
                      Y3: [Y4, X4],
                      Z3: [Z4, Y4]
                     },
                     {W1: [X2, Y2],
                      X1: [Y2, Z2],
                      Y1: [Z2, W2],
                      Z1: [W2, X2],
                      W2: [X3, Y3],
                      X2: [Y3, Z3],
                      Y2: [Z3, W3],
                      Z2: [W3, X3],
                      W3: [X4, Y4],
                      X3: [Y4, Z4],
                      Y3: [Z4, W4],
                      Z3: [W4, X4]
                     }
)
gOut = ID(Z4, W1, G7)
