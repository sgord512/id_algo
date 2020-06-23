import igraph
from igraph import Graph as IGraph

import utilities

class Graph(IGraph):

    TYPE_CAUSAL_MODEL = 1
    TYPE_DAG = 2
    TYPE_C_SKELETON = 3

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
        # I need to do this after adding in all the edges.
        g._construct_observed_and_confounded()
        return g

    def __init__(self, graph_type=None):
        super().__init__(directed=True)
        self.obs = None
        self.confounded = None
        self.graph_type = graph_type or Graph.TYPE_CAUSAL_MODEL
        if self.graph_type == Graph.TYPE_CAUSAL_MODEL:
            self._construct_observed_and_confounded()

    def add_edges_from(self, source, targets):
        self._add_edges_base([(source, v) for v in targets], observed=True)

    def _add_edges_base(self, edge_list, **kwds):
        for (u,v) in edge_list:
            super().add_edge(u, v, **kwds)

    def add_observed_edge(self, u, v):
        self.add_observed_edges([(u,v)])

    def add_observed_edges(self, edge_list):
        self._add_edges_base(edge_list, observed=True)

    def build_vertex_id_to_name_map(self):
        names = dict()
        for (i,v) in enumerate(self.vs):
            names[i] = v["name"]
        return names

    def _convert_edge_tuples_to_named_pairs(self, edges):
        names = self.build_vertex_id_to_name_map()
        return list(map(lambda tuple: (names[tuple[0]], names[tuple[1]]), edges))

    def _get_observed_edges(self):
        if len(self.es) == 0: return []
        obs_edges = list(map(lambda e: e.tuple, self.es.select(observed_eq=True)))
        return self._convert_edge_tuples_to_named_pairs(obs_edges)

    def _get_hidden_edges(self):
        if len(self.es) == 0: return []
        hidden_edges = list(map(lambda e: e.tuple, self.es.select(observed_eq=False)))
        return self._convert_edge_tuples_to_named_pairs(hidden_edges)

    def _construct_observed_subgraph(self):
        assert(self.graph_type==Graph.TYPE_CAUSAL_MODEL
               or self.graph_type==Graph.TYPE_DAG)
        names = self.build_vertex_id_to_name_map()
        G = Graph(graph_type=Graph.TYPE_DAG)
        G.add_vertices(list(names.values()))
        G.add_observed_edges(self._get_observed_edges())
        return G

    def _construct_confounded_subgraph(self):
        assert(self.graph_type==Graph.TYPE_CAUSAL_MODEL
               or self.graph_type==Graph.TYPE_C_SKELETON)
        names = self.build_vertex_id_to_name_map()
        G = Graph(graph_type=Graph.TYPE_C_SKELETON)
        G.add_vertices(list(names.values()))
        hidden_edges = self._get_hidden_edges()
        G.add_hidden_edges(utilities.unique_tuples(hidden_edges))
        return G

    def get_post_intervention_edges(self, v):
        names = self.build_vertex_id_to_name_map()
        hidden_edge_list = []
        observed_edge_list = []
        for edge in self.es:
            source, target = edge.tuple
            if (names[target] in v or
                (names[source] in v and edge["observed"] == False)):
                continue
            if edge["observed"]:
                observed_edge_list.append((source, target))
            else:
                hidden_edge_list.append((source, target))
        return observed_edge_list, utilities.unique_tuples(hidden_edge_list)

    def construct_post_intervention_subgraph(self, v):
        names = self.build_vertex_id_to_name_map()
        G = Graph()
        observed_edge_list, hidden_edge_list = self.get_post_intervention_edges(v)
        G.add_vertices(list(names.values()))
        G.add_observed_edges(observed_edge_list)
        G.add_hidden_edges(hidden_edge_list)
        G._construct_observed_and_confounded()
        return G

    def add_hidden_edges(self, edge_list):
        full_edge_list = Graph._double_and_flip_list(edge_list)
        self._add_edges_base(full_edge_list, observed=False)

    def induced_subgraph(self, v):
        names = self.build_vertex_id_to_name_map()
        hidden_edge_list = []
        observed_edge_list = []
        for edge in self.es:
            source, target = edge.tuple
            if (names[target] not in v or
                (names[source] not in v)):
                continue
            if edge["observed"]:
                observed_edge_list.append((names[source], names[target]))
            else:
                hidden_edge_list.append((names[source], names[target]))
        G = Graph()
        G.add_vertices(list(v))
        G.add_observed_edges(observed_edge_list)
        G.add_hidden_edges(utilities.unique_tuples(hidden_edge_list))
        G._construct_observed_and_confounded()
        return G

    def __str__(self):
        return super().__str__()

    def __repr__(self):
        return super().__repr__()

    def get_c_components(self):
        components = self.confounded.components(mode=igraph.STRONG)
        return [set(self.confounded.vs[c]["name"]) for c in components]

    # Get connected components of a graph.
    def get_components(self):
        components = self.components(mode=igraph.STRONG)
        return [set(self.vs[c]["name"]) for c in components]

    def get_ancestors(self, y):
        assert(self.obs is not None)
        # This returns the ancestors of the vertices in the set y
        n = self.get_num_nodes()
        ancestorIndices = utilities.flatten(self.obs.neighborhood(y, order=n, mode=igraph.IN))
        return set(self.obs.vs[ancestorIndices]["name"])

    def _construct_observed_and_confounded(self):
        assert(self.graph_type==Graph.TYPE_CAUSAL_MODEL)
        self.obs = self._construct_observed_subgraph()
        self.confounded = self._construct_confounded_subgraph()

    def get_all_connected_vertices(self, y):
        n = self.get_num_nodes()
        reachableIndices = utilities.flatten(self.neighborhood(y, order=n, mode=igraph.ALL))
        return set(self.vs[reachableIndices]["name"])

    def get_topological_order(self):
        return self.vs[self.obs.topological_sorting(mode=igraph.OUT)]["name"]

    def get_vertices(self):
        return self.vs["name"]
