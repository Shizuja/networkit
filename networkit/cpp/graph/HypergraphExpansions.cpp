/*
 * HypergraphExpansions.cpp
 *
 *  Created on: 11.09.2024
 *      Author: Yanneck Dimitrov
 */
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/HypergraphExpansions.hpp>
#include <networkit/centrality/Betweenness.hpp>
#include <map>
#include <set>
#include <vector>
#include <iostream>

namespace NetworKit {
    
Graph HypergraphExpansions::cliqueExpansion(Hypergraph &hypergraph) {
    /*
    Converts a hypergraph into its clique expansion (a simple graph)

    Parameters
    ----------
    hypergraph - NetworKit::Hypergraph
        input hypergraph

    Return
    ----------
    cliqueExpansion - NetworKit::Graph
        clique expansion of hypergraph
    */

    //clique expansion is a simple networkit graph where each hyperedge is a clique
    NetworKit::Graph cliqueExpansion(hypergraph.numberOfNodes());

    //iterate over all nodes
    hypergraph.forNodes([&](NetworKit::node node) {
        //create a clique by connecting each pair of nodes in the hyperedge
        std::set<NetworKit::node> neighbors = hypergraph.getNeighbors(node);
        for(NetworKit::node neighbor : neighbors) {
            if(!cliqueExpansion.hasEdge(node, neighbor)){
                cliqueExpansion.addEdge(node, neighbor);
            };
        };
    });

    return cliqueExpansion;
}

std::pair<NetworKit::Graph, std::map<NetworKit::node, std::pair<NetworKit::node, NetworKit::edgeid>>> HypergraphExpansions::lineExpansion(NetworKit::Hypergraph &hypergraph) {
    /*
    Converts a hypergraph into its line expansion (a simple graph)

    Parameters
    ----------
    hypergraph - NetworKit::Hypergraph
        input hypergraph

    Return
    ----------
    lineExpansion - NetworKit::Graph
        line expansion of hypergraph
    */

    //line expansion is a simple networkit graph of which we still need to dertermin the amount of its nodes by mapping all combinations of a node with an edge which they are part of 
    std::map<NetworKit::node, std::pair<NetworKit::node, NetworKit::edgeid>> nodeMap;

    //ID for the newest node to be added
    NetworKit::node newNodeID = 0;

    //iterate over all hyperedges and their nodes to determine all pairs of a node with its hyperedges
    hypergraph.forEdges([&](NetworKit::edgeid eid) {
        //std::cout << "current Edge: " << eid << std::endl;
        hypergraph.forNodes([&](NetworKit::node node) {
            //std::cout << "current node: " << node << std::endl;
            if(hypergraph.hasNode(node, eid)) {
                //add pair to map if a node is part of a hyperedge
                nodeMap[newNodeID] = {node, eid};
                //std::cout << newNodeID << " " << node << " " << eid << std::endl;
                newNodeID++;
            }
        });
    });

    //now create the lineExpansion graph since we now know the amount of nodes in it
    NetworKit::Graph lineExpansion(nodeMap.size());

    //add edges to lineExpansion based on entries from nodeMap
    for(auto& iterator1 : nodeMap) {
        for(auto& iterator2 : nodeMap) {
            //access the pair of data from the original hypergraph in each nodeMap entry
            auto pair1 = iterator1.second;
            auto pair2 = iterator2.second;
            //std::cout << "current Pair 1 (Node, Edgeid): " << pair1.first << pair1.second << std::endl;
            //std::cout << "current Pair 2 (Node, Edgeid): " << pair2.first << pair2.second << std::endl;
            //only check pair that are not the same
            if(pair1.first != pair2.first || pair1.second != pair2.second){
                //std::cout << "Not the same!" << std::endl;
                //add edge if a pair of entries in nodeMap refer to the same node or hyperedge in the original hypergraph
                if(pair1.first == pair2.first || pair1.second == pair2.second){
                    //std::cout << "Connected!" << std::endl;
                    if(!lineExpansion.hasEdge(iterator1.first, iterator2.first)) {
                        lineExpansion.addEdge(iterator1.first, iterator2.first);
                    }
                }
            }
        }
    }
    return {lineExpansion, nodeMap};
}

NetworKit::Hypergraph HypergraphExpansions::reconstructHypergraphFromLineExpansion(NetworKit::Graph &lineExpansionGraph, std::map<NetworKit::node, std::pair<NetworKit::node, NetworKit::edgeid>> &nodeMap) {
    /*
    Reconstructs the original hypergraph from its line expansion

    Parameters
    ----------
    lineExpansionGraph - NetworKit::Graph
        line expansion graph of a hypergraph

    nodeMap - std::map<NetworKit::node, std::pair<NetworKit::node, NetworKit::edgeid>>
        node id's from the line expansion graph mapped to a pair of node and edgeid from the original hypergraph

    Return
    ----------
    hypergraph - NetworKit::Hypergraph
        reconstructed hypergraph
    */

    //set for the nodes of the hypergraph
    std::set<NetworKit::node> nodes;
    //map for the hyperedges of the hypergraph
    std::map<NetworKit::edgeid, std::set<NetworKit::node>> hyperedges;

    //iterate over nodes from the line expansion graph
    lineExpansionGraph.forNodes([&](NetworKit::node node) {
        //get original node and hyperedge from nodeMap
        std::pair<NetworKit::node, NetworKit::edgeid> entry = nodeMap.at(node);
        //std::cout << "Current Node: " << entry.first << " Current Edge: " << entry.second << std::endl;
        //add original node to corresponding hyperedge set
        hyperedges[entry.second].insert(entry.first);
        //std::cout << "Node " << entry.first << " inserted into hyperedge " << entry.second << std::endl;
        //if node does not exist in hypergraph already, add it
        if (nodes.find(entry.first) == nodes.end()) {
            nodes.insert(entry.first);
        }
    });

    //create hypergraph and add edges from dictionary
    NetworKit::Hypergraph hypergraph(nodes.size());
    //convert sets into vectors for hyperedges
    for (int i=0; i<hyperedges.size(); i++) {
        std::set set = hyperedges[i];
        std::vector<NetworKit::node> vec(set.begin(), set.end());
        hypergraph.addEdge(vec);
    }
    return hypergraph;
}

size_t HypergraphExpansions::getIntersectionSize(Hypergraph &hypergraph, edgeid eid1, edgeid eid2, node node) {
    if(hypergraph.hasNode(node, eid1) && hypergraph.hasNode(node, eid2)) {
        size_t intersection_size = 1;
        for (NetworKit::node neighbor : hypergraph.getNeighbors(node)) {
            if(hypergraph.hasNode(neighbor, eid1) && hypergraph.hasNode(neighbor, eid2)) {
                intersection_size++;
            }
        }
        return intersection_size;
    }
    return 0;
}

std::vector<nodeweight> HypergraphExpansions::lineExpansionWeightedBetweenness(Graph &G, std::map<node, std::pair<node, edgeid>> &nodeMap, bool normalized) {

    Betweenness centrality(G, normalized);
    centrality.run();
    Hypergraph H = HypergraphExpansions::reconstructHypergraphFromLineExpansion(G, nodeMap);
    std::vector<nodeweight> betweennessScores(H.numberOfNodes());
    std::map<node, std::set<edgeid>> intersection;

    for (std::pair<const node, std::pair<node, edgeid>> entry : nodeMap) {
        intersection[entry.second.first].insert(entry.second.second);
        betweennessScores[entry.second.first] += centrality.scores().at(entry.first);
    }
    H.forNodes([&](node node) {
        if(intersection[node].size() == 1) {
            //nothing needed
        } else if(intersection[node].size() == 2) {
            edgeid eid1 = *intersection[node].begin();
            edgeid eid2 = *std::next(intersection[node].begin(), 1);
            betweennessScores[node] /= HypergraphExpansions::getIntersectionSize(H, eid1, eid2, node);
        } else {
            throw std::invalid_argument("Hypergraph contains nodes which are part of more than 2 hyperedges. Currently this function does not support such hypergraphs.");
        }
    });
    return betweennessScores;
}

std::vector<nodeweight> HypergraphExpansions::lineExpansionBetweenness(Graph &G, std::map<node, std::pair<node, edgeid>> &nodeMap, bool normalized, bool additive) {

    Betweenness centrality(G, normalized);
    centrality.run();
    std::vector<nodeweight> betweennessScores(HypergraphExpansions::numberOfNodesFromNodeMap(nodeMap));
    std::vector<size_t> memberOfHyperedges(HypergraphExpansions::memberOfHyperedges(nodeMap));

    for(std::pair<const node, std::pair<node, edgeid>> entry : nodeMap) {
        betweennessScores[entry.second.first] += centrality.scores().at(entry.first);
    } 
    if(!additive) {
        for (size_t i = 0; i < betweennessScores.size(); i++) {
            betweennessScores[i] /= memberOfHyperedges[i];
        }
    }    
    return betweennessScores;
}

size_t HypergraphExpansions::numberOfNodesFromNodeMap(std::map<node, std::pair<node, edgeid>> &nodeMap) {
    std::set<node> seen;
    size_t numberOfNodes = 0;
    for(std::pair<node, std::pair<node, edgeid>> entry : nodeMap) {
        if(seen.find(entry.second.first) == seen.end()) {
            seen.insert(entry.second.first);
            numberOfNodes++;
        }
    }
    return numberOfNodes;
}

std::vector<size_t> HypergraphExpansions::memberOfHyperedges(std::map<node, std::pair<node, edgeid>> &nodeMap) {
    std::vector<size_t> values(HypergraphExpansions::numberOfNodesFromNodeMap(nodeMap));
    for(std::pair<node, std::pair<node, edgeid>> entry : nodeMap) {
        values[entry.second.first]++;
    }
    return values;
}

} //namespace NetworKit
