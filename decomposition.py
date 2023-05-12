import networkx as nx
import sys, os
from optparse import OptionParser

def onion_decomposition(graph):
    # Creates a copy of the graph (to be able to remove vertices and edges)
    the_graph = nx.Graph(graph).to_undirected()
    #the_graph.remove_edges_from( the_graph.selfloop_edges() )
    # Dictionaries to register the k-core/onion decompositions. 
    coreness_map = {}
    layer_map = {}
    # Performs the onion decomposition.
    current_core = 1
    current_layer = 1
    while the_graph.number_of_nodes() > 0:
        # Sets properly the current core.
        degree_sequence = [d for n,d in the_graph.degree()]
        min_degree =  min(degree_sequence)
        if min_degree >= (current_core+1):
            current_core = min_degree
        # Identifies vertices in the current layer.
        this_layer_ = []
        for v in the_graph.nodes():
            if the_graph.degree()[v] <= current_core:
                this_layer_.append(v)
        # Identifies the core/layer of the vertices in the current layer.
        for v in this_layer_:
            coreness_map[v] = current_core
            layer_map[v] = current_layer
            the_graph.remove_node(v)
        # Updates the layer count.
        current_layer = current_layer + 1
    # Returns the dictionaries containing the k-shell and onion layer of all vertices.
    return (layer_map, coreness_map)

def team_structure(graph):
    # Dictionary to register the k-core/onion decompositions.
    team_size_map = {}
    # Loop over all nodes and extract team size from node label
    for v in graph.nodes():
        x = v.split("_")
        team_size_map[v] = x[1]
    # Returns the dictionaries containing the team size of all vertices.
    return team_size_map

def merge_nodes(G,nodes, new_node):
    # Add the merged node
    G.add_node(new_node)
    # For all edges related to one of the nodes to merge,
    # make an edge going to or coming from the new node.
    for n1,n2 in G.edges():
        if n1 in nodes:
            G.add_edge(new_node,n2)
        elif n2 in nodes:
            G.add_edge(n1,new_node)
    # remove the merged nodes
    for n in nodes:
        G.remove_node(n)

def clique_structure(graph):
    # Dictionary to register the k-core/onion decompositions.
    clique_size_map = {}
    # Initialize with teams of 1
    for v in graph.nodes():
        clique_size_map[v] = 1
    # Create list of cliques in arbitrary order
    cliques_it = nx.find_cliques(G)
    cliques = []
    for c in cliques_it:
        if len(c)>2:
            cliques.append(c)
    # Sort cliques from largest to smallest
    cliques.sort(key = len, reverse=True)
    # Merge cliques if nodes haven't been assigned to clique yet
    for clique in cliques:
        new_clique = True
        for node in clique:
            if graph.has_node(node) == False:
                new_clique = False
        if new_clique:
            clique_size_map[" ".join(str(x) for x in clique)] = len(clique)
            merge_nodes(graph,clique," ".join(str(x) for x in clique))

    # Returns the dictionaries containing the team size of all vertices.
    return graph, clique_size_map

def update_sparse_matrix(dictionary, values, key, index):
    # Creates entry if it doesn't exist or else update the corresponding value
    if(not key in dictionary):
        index_i = index
        dictionary[key] = index_i
        values.append(1)
        index += 1
    else:
        index_i = dictionary[key]
        values[index_i]+=1
    # Returns new index (i.e. the nb of nonzero entry in sparse matrix)
    return index

def create_sparse_matrix_for_graph(onion, kcore, team, graph):
    # Store the index of every type, and the number of elements of that type
    cone_sparse_node_matrix = {}
    cone_sparse_edge_matrix = {}
    node_values = []
    node_index = 0
    edge_values = []
    edge_index = 0
    # Caracterize nodes
    for n in graph.nodes():
        s1 = team[n]
        d1 = graph.degree(n)
        l1 = onion[n]
        c1 = kcore[n]
        node_index = update_sparse_matrix(cone_sparse_node_matrix, node_values, (s1,d1,l1,c1), node_index)
    # Caracterize edges
    for e in graph.edges():
        d1 = graph.degree(e[0])
        d2 = graph.degree(e[1])
        l1 = onion[e[0]]
        l2 = onion[e[1]]
        edge_index = update_sparse_matrix(cone_sparse_edge_matrix,edge_values,(d1,l1,d2,l2),edge_index)
    return(cone_sparse_node_matrix, node_values, cone_sparse_edge_matrix, edge_values)


if __name__ == '__main__':
    parser = OptionParser(usage="usage: %prog [options]")
    parser.add_option("-n", "--network", dest="n", default=None, type=str, help="Specify path to edgelist")
    (options, args) = parser.parse_args()
    G = nx.read_edgelist(options.n, delimiter=',', nodetype=str)
    G = nx.Graph(G).to_undirected()
    G.remove_edges_from(nx.selfloop_edges(G))

    # Analyze clique structure
    (G, team_size_map) = clique_structure(G)
    # Simplify network after merging of cliques
    G = nx.Graph(G).to_undirected()
    G.remove_edges_from(nx.selfloop_edges(G))
    # Analyze onion structure
    (layer_map, coreness_map) = onion_decomposition(G)
    # Classify nodes and edges based on CONE types
    (node_matrix, node_values, edge_matrix, edge_values) = create_sparse_matrix_for_graph(layer_map, coreness_map, team_size_map, G)

    # Output type counts for nodes and edges
    fout1 = "./matrices/Pnkl.txt"
    fo1 = open(fout1, "w")
    for k, v in node_matrix.items():
        fo1.write(str(k[0]) + '  '+ str(k[1]) + '  '+ str(k[2]) + '  '+ str(k[3]) + '  '+ str(node_values[node_matrix[k]]) + '\n')
    fo1.close()
    fout2 = "./matrices/Lmat.txt"
    fo2 = open(fout2, "w")
    for k, v in edge_matrix.items():
        fo2.write(str(k[0]) + '  '+ str(k[1]) + '  '+ str(k[2]) + '  '+ str(k[3]) + '  '+ str(edge_values[edge_matrix[k]]) + '\n')
    fo2.close()

