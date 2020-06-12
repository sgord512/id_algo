import igraph
from igraph import Graph as IGraph
from collections import defaultdict
from id_algo_implementation import Graph, flatten
import copy
import utilities

class CADMG(Graph):
    # This class represents a CADMG. I'm going to try to be vigilant about all methods always return a set of names of vertices as opposed to raw indices.

    def __init__(self):
        super().__init__()

    @staticmethod
    def FromDicts(random_vertices, fixed_vertices, observed_edge_dict=None, hidden_edge_dict=None):
        observed_edges = Graph._edge_list_from_dict(observed_edge_dict) if observed_edge_dict else []
        hidden_edges = Graph._edge_list_from_dict(hidden_edge_dict) if hidden_edge_dict else []
        g = CADMG()
        g.add_vertices(random_vertices)
        g.vs.select(name_in=random_vertices)["fixed"]=False
        g.add_vertices(fixed_vertices)
        g.vs.select(name_in=fixed_vertices)["fixed"]=True
        g.add_observed_edges(observed_edges)
        g.add_hidden_edges(hidden_edges)
        return g

    def get_indices_from_vertices(self, vertices):
        return self.vs.select(name_in=vertices)

    def get_random_vertices(self):
        return set(self.vs.select(fixed_eq=False)['name'])

    def get_fixed_vertices(self):
        return set(self.vs.select(fixed_eq=True)['name'])

    def is_random(self, v):
        return v in self.get_random_vertices()

    def is_fixed(self, v):
        return v in self.get_fixed_vertices()

    def generate_confounded(self):
        return self.subgraph_edges(self.es.select(observed_eq=False), delete_vertices=False)

    def generate_visible(self):
        return self.subgraph_edges(self.es.select(observed_eq=True), delete_vertices=False)

    def get_districts(self):
        confounded = self.generate_confounded()
        components = confounded.components(mode=igraph.STRONG)
        return [set(self.vs[c]["name"]) for c in components]

    def get_district(self, v):
        districts = self.get_districts()
        for district in districts:
            if v in district:
                return district
        raise Exception("Vertex doesn't have a district!")

    def get_descendants(self, v):
        n = self.get_num_nodes()
        visible = self.generate_visible()
        ancestorIndices = utilities.flatten(visible.neighborhood(v, order=n, mode=igraph.OUT))
        return set(visible.vs[ancestorIndices]["name"])

    def get_parents(self, v):
        return set(self.neighbors(v, mode=igraph.IN)["name"])

    def get_markov_blanket(self, v):
        nodes = set()
        district = self.get_district(v)
        for u in district:
            nodes = nodes | self.get_parents(u)
        nodes = nodes | district
        return nodes - set([v])

    def is_fixable(self, v):
        if self.is_fixed(v):
            return False
        else:
            descendants = self.get_descendants(v)
            district = self.get_district(v)
            return len(descendants & district) == 1

    def fix_vertex(self, v):
        G = self.deepcopy()
        #G.es.select()

    def __str__(self):
        return super().__str__()

    def __repr__(self):
        return super().__repr__()

# Initially, I'll have Kernels not be recursive. Each time I generate a new Kernel, I'm going to give it a distinct label and then I can carefully things by only allowing a single method for constructing a new expression from a current kernel.
class Kernel():
    def __init__(self, V, W, label=None):
        self.V = set(V)
        self.W = set(W)
        self.Vprime = None
        self.Wprime = None
        self.label = label or q

    def is_raw(self):
        return self.Vprime is None and self.Wprime is None

    def is_specialized(self):
        return not self.is_raw()

    def specialize(self, Vprime, Wprime):
        if self.is_specialized():
            raise Exception("Can't specialize a kernel that is already specialized!")
        elif self.V >= Vprime and self.W <= Wprime:
            raise Exception("Can't construct a kernel in this way!")
        else:
            spec_kernel = Kernel(self.V, self.W, label=self.label)
            spec_kernel.Vprime = Vprime
            spec_kernel.Wprime = Wprime
            return spec_kernel

    @staticmethod
    def MatchesCADMG(ker, G):
        return False

class CADMGKernel():
    def __init__(self, G, kernel):
        self.G = G
        self.kernel = kernel

    def fix_vertex(self, v):
        if not self.G.is_fixable(v):
            raise Exception("Chosen vertex isn't fixable.")
        

def IDAlgo(cadmgKernel, Y, A):
    Ystar = None
    ## To Be Continued



if __name__ == "__main__":
    A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,W,X,Y,Z = "a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,w,x,y,z".split(',')
    X1,X2,X3,X4,X5 = "x_1,x_2,x_3,x_4,x_5".split(',')
    Y1,Y2,Y3,Y4,Y5 = "y_1,y_2,y_3,y_4,y_5".split(',')
    Z1,Z2,Z3,Z4,Z5 = "z_1,z_2,z_3,z_4,z_5".split(',')
    W1,W2,W3,W4,W5 = "w_1,w_2,w_3,w_4,w_5".split(',')

    SingleEdge = CADMG.FromDicts([Y], X, {X: Y}, {})
