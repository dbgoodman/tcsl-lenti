import leidenalg
import igraph as ig
import pandas as pd
import numpy as np
import math

def leiden_clust(nearest_neighbors, edge_weights):

    edges = []
    n_neighbors = nearest_neighbors.shape[1]

    for i, nn_row in enumerate(nearest_neighbors):
        edge_weights_i = edge_weights[i]
        
        for j in range(1, n_neighbors):
            edges.append((int(nn_row[0]), int(nn_row[j]), edge_weights_i[j]))
    
    knn_graph = ig.Graph.TupleList(edges, directed=True, weights=True)
    
    part = leidenalg.find_partition(knn_graph, 
            leidenalg.ModularityVertexPartition, weights='weight')
    
    cluster_membership = [x for _,x in 
            sorted(zip(knn_graph.vs['name'], part.membership))]
            
    return(cluster_membership)
