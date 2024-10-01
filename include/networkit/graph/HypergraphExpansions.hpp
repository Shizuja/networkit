#ifndef NETWORKIT_GRAPH_HYPERGRAPH_EXPANSIONS_HPP_
#define NETWORKIT_GRAPH_HYPERGRAPH_EXPANSIONS_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/HypergraphTools.hpp>
#include <map>
#include <set>
#include <vector>
#include <iostream>

namespace NetworKit {
namespace HypergraphExpansions {

/*
    Converts a hypergraph into its clique expansion (a simple graph)

    @param hypergraph Input hypergraph
    @return cliqueExpansion Clique expansion of hypergraph
    
*/
Graph cliqueExpansion(Hypergraph &hypergraph);

/*
    Converts a hypergraph into its line expansion (a simple graph)

    @param hypergraph Input hypergraph
    @return lineExpansion Line expansion of hypergraph
    @return nodeMap Map with nodes from lineExpansion as keys and {node,hyperedge} from hypergraph as values
    
*/
std::pair<Graph, std::map<node, std::pair<node, edgeid>>> lineExpansion(Hypergraph &hypergraph);

/*
    Reconstructs the original hypergraph from its line expansion

    @param lineExpansionGraph Input line expansion graph
    @param nodeMap Map with nodes from lineExpansionGraph as keys and {node,hyperedge} from hypergraph as values
    @return hypergraph Reconstructed hypergraph from line expansion
*/
Hypergraph reconstructHypergraphFromLineExpansion(Graph &lineExpansionGraph, std::map<node, std::pair<node, edgeid>> &nodeMap);

/*
    Returns the intersection of two hyperedges as a set

    @param hypergraph Input hypergraph
    @param eid1 first hyperedge
    @param eid2 second hyperedge
    @return intersection
*/
std::set<node> getIntersection(Hypergraph &hypergraph, edgeid eid1, edgeid eid2);

/*
    Returns a map with weighted centrality scores for each node based on the line expansion of a hypergraph and the intersection size of its hyperedges

    @param G line expansion graph of a hypergraph
    @param nodeMap node mapping for G
    @param normalized Set this parameter to true if scores should be normalized in the interval [0,1].
    @return betweenness scores
*/
std::vector<nodeweight> lineExpansionWeightedBetweenness(Graph &G, std::map<node, std::pair<node, edgeid>> &nodeMap, bool normalized = false);

/*
    Returns a map with centrality scores for each node based on the line expansion of a hypergraph

    @param G line expansion graph of a hypergraph
    @param nodeMap node mapping for G
    @param normalized Set this parameter to true if scores should be normalized in the interval [0,1].
    @param additive If true all centrality scores representing one node are summed up, if false an average will be calculated
    @return betweenness scores
*/
std::vector<nodeweight> lineExpansionBetweenness(Graph &G, std::map<node, std::pair<node, edgeid>> &nodeMap, bool normalized = false, bool additive = true);


/*
    Returns the number of nodes in a hypergraph from the nodeMap of its line Expansion

    @param nodeMap node mapping of the line expansion of a hypergraph
    @return number of nodes
*/
size_t numberOfNodesFromNodeMap(std::map<node, std::pair<node, edgeid>> &nodeMap);

/*
    Returns a vector with the amounts of hyperedges each nodes is part of

    @param nodeMap node mapping of the line expansion of a hypergraph
    @return vector with amount of hyperedges each node is part of
*/
std::vector<size_t> memberOfHyperedges(std::map<node, std::pair<node, edgeid>> &nodeMap);

/*
    Returns map of node members of a hyperedge

    @param hypergraph Input hypergraph
    @param edgeid Input hyperedge
    @return map of node members of eid
    
*/
std::set<node> getEdgeMembers(Hypergraph &hypergraph, edgeid eid);

/*
    Returns map of node members of each hyperedge

    @param hypergraph Input hypergraph
    @return map of node members of each hyperedge
    
*/
std::map<edgeid, std::set<node>> getAllEdgeMembers(Hypergraph &hypergraph);

/*
    Converts a hypergraph into its line graph (a simple graph)

    @param hypergraph Input hypergraph
    @param weighted should the edges of the line graph be weighted based on the sizes and intersection sizes of the hyperedges
    @return lineGraph line graph of hypergraph
    
*/
Graph lineGraph(Hypergraph &hypergraph, bool weighted = true);

/*
    Converts a hypergraph into its line graph (a simple graph) and computes betweenness centrality scores for its nodes

    @param hypergraph Input hypergraph
    @param normalized Set this parameter to true if scores should be normalized in the interval [0,1].
    @return betweenness scores
    
*/
std::vector<nodeweight> lineGraphBetweenness(Hypergraph &hypergraph, bool normalized = false);

} // namespace HypergraphExpansions

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_HYPERGRAPH_EXPANSIONS_HPP_
