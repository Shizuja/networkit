/*
 * HypergraphExpansions.cpp
 *
 *  Created on: 11.09.2024
 *      Author: Yanneck Dimitrov
 */
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/HypergraphExpansions.hpp>
#include <networkit/centrality/Betweenness.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <map>
#include <set>
#include <vector>
#include <limits>
#include <iostream>

namespace NetworKit {
    
Graph HypergraphExpansions::cliqueExpansion(Hypergraph &hypergraph) {

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

    //line expansion is a simple networkit graph of which we still need to dertermin the amount of its nodes by mapping all combinations of a node with an edge which they are part of 
    std::map<NetworKit::node, std::pair<NetworKit::node, NetworKit::edgeid>> nodeMap;

    //ID for the newest node to be added
    NetworKit::node newNodeID = 0;

    //iterate over all hyperedges and their nodes to determine all pairs of a node with its hyperedges
    hypergraph.forEdges([&](NetworKit::edgeid eid) {
        //std::cout << "Step: " << eid << " / " << hypergraph.numberOfEdges() << std::endl;
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
        //std::cout << "Step: " << iterator1.first << " / " << nodeMap.size() << std::endl;
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

std::set<node> HypergraphExpansions::getIntersection(Hypergraph &hypergraph, edgeid eid1, edgeid eid2) {
    std::set<node> set_eid1 = HypergraphExpansions::getEdgeMembers(hypergraph, eid1);
    std::set<node> set_eid2 = HypergraphExpansions::getEdgeMembers(hypergraph, eid2);
    std::set<node> intersection;
    for (node node : set_eid1) {
        if (set_eid2.find(node) != set_eid2.end()) {
            intersection.insert(node);
        }
    }
    return intersection;
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
            betweennessScores[node] /= HypergraphExpansions::getIntersection(H, eid1, eid2).size();
        } else {
            throw std::invalid_argument("Hypergraph contains nodes which are part of more than 2 hyperedges. Currently this function does not support such hypergraphs.");
        }
    });
    return betweennessScores;
}

std::vector<nodeweight> HypergraphExpansions::lineExpansionBetweenness(Graph &G, std::map<node, std::pair<node, edgeid>> &nodeMap, bool normalized, bool additive) {
    //std::cout << "Start of betweenness calculation" << std::endl;
    Betweenness centrality(G, normalized);
    centrality.run();
    std::vector<nodeweight> betweennessScores(HypergraphExpansions::numberOfNodesFromNodeMap(nodeMap), 0);
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
    for (size_t i = 0; i < values.size(); i++) {
        if (values[i] < 1)
        {
            values[i] = 1;
            std::cout << "Following node is not part of any hyperedge: " << i << std::endl;
        }
    }
    
    return values;
}

std::set<node> HypergraphExpansions::getEdgeMembers(Hypergraph &hypergraph, edgeid eid) {
    std::set<node> members;
    hypergraph.forNodes([&](node node) {
            if (hypergraph.hasNode(node, eid)) {
                members.insert(node);
            }
    });
    return members;
}

std::map<edgeid, std::set<node>> HypergraphExpansions::getAllEdgeMembers(Hypergraph &hypergraph) {
    std::map<edgeid, std::set<node>> edge_members;
    hypergraph.forEdges([&](edgeid eid) {
        std::set<node> members;
        hypergraph.forNodes([&](node node) {
            if (hypergraph.hasNode(node, eid)) {
                members.insert(node);
            }
        });
        edge_members[eid] = members;
    });
    return edge_members;
}

Graph HypergraphExpansions::lineGraph(Hypergraph &hypergraph, bool weighted) {
    Graph lineGraph(hypergraph.numberOfEdges(), weighted);
    std::set<edgeid> done;

    if(!weighted) {
        hypergraph.forEdges([&](edgeid eid1){
            hypergraph.forEdges([&](edgeid eid2){
                //do not iterate over (1,0) and (0,1) but just one of them (for (eid1,eid2))
                if(done.find(eid2) == done.end()) {
                    if(!HypergraphExpansions::getIntersection(hypergraph, eid1, eid2).empty()) {
                        if(!lineGraph.hasEdge(eid1,eid2)) {
                            lineGraph.addEdge(eid1, eid2);
                        }
                    }
                }
            });
            done.insert(eid1);
        });
    } else {
        hypergraph.forEdges([&](edgeid eid1){
            //std::cout << "Step: " << eid1 << " / " << hypergraph.numberOfEdges() << std::endl;
            hypergraph.forEdges([&](edgeid eid2){
                //do not iterate over (1,0) and (0,1) but just one of them (for (eid1,eid2))
                if(done.find(eid2) == done.end()) {
                    std::set<NetworKit::node> nodes_in_eid1 = HypergraphExpansions::getEdgeMembers(hypergraph, eid1);
                    std::set<NetworKit::node> nodes_in_eid2 = HypergraphExpansions::getEdgeMembers(hypergraph, eid2);
                    double intersection_size = double(HypergraphExpansions::getIntersection(hypergraph, eid1, eid2).size());
                    double union_size = double(nodes_in_eid1.size() + nodes_in_eid2.size() - intersection_size);
                    if(intersection_size > 0) {
                        if(!lineGraph.hasEdge(eid1,eid2)) {
                            edgeweight weight = (1.0/3.0)*(union_size + (union_size/intersection_size));
                            //std::cout << weight << " = 1/3 * (" << union_size << " + (" << union_size << "/" << intersection_size << ")) - 1" << std::endl;
                            lineGraph.addEdge(eid1, eid2, weight);
                        }
                    }
                }
            });
            done.insert(eid1);
        });
    }
    return lineGraph;
}

std::vector<nodeweight> HypergraphExpansions::lineGraphBetweenness_alt(Hypergraph &hypergraph, bool normalized) {
    
    Graph lineGraph = HypergraphExpansions::lineGraph(hypergraph, true);
    std::vector<nodeweight> centrality_scores(hypergraph.numberOfNodes());
    size_t number_of_shortest_paths = 0;

    //Get a map of all edges a node is part of
    std::map<node, std::set<edgeid>> edges;
    hypergraph.forNodes([&](node node){
        hypergraph.forEdges([&](edgeid eid){
            if(hypergraph.hasNode(node, eid)) {
                edges[node].insert(eid);
            }
        });
    });

    //calculate centrality
    hypergraph.forNodes([&](node node1){
        hypergraph.forNodes([&](node node2){
            if(node1 != node2) {
                //variables for getting the shortest paths on the line graph
                double shortest_length = std::numeric_limits<double>::max();
                std::set<std::vector<node>> shortest_paths;
                edgeweight distance;
                //all edges node1 is part of
                for(auto edge1 : edges[node1]) {
                    Dijkstra dijkstra(lineGraph, edge1, true, false);
                    dijkstra.run();
                    //all edges node2 is part of
                    for(auto edge2 : edges[node2]) {
                        distance = dijkstra.distance(edge2);
                        //search shortest path
                        if (distance < shortest_length) {
                                shortest_length = distance;
                                shortest_paths = dijkstra.getPaths(edge2);
                        } else if (distance == shortest_length) {
                            for(auto path : dijkstra.getPaths(edge2)) {
                                shortest_paths.insert(path);
                            }
                        }
                    }
                }
                //Now that we have the shortest paths between node 1 and 2, do the betweenness scores
                for(std::vector<node> path : shortest_paths) {
                    for (size_t i = 0; i < path.size()-1; i++){
                        std::set<node> intersection = HypergraphExpansions::getIntersection(hypergraph, path.at(i), path.at(i+1));
                        //add 1 / intersection_size to centrality scores for each node in intersection and divide by number of paths if scores should be normalized
                        for(node node : intersection) {
                            centrality_scores.at(node) += 1.0 / double(intersection.size() * shortest_paths.size());
                            number_of_shortest_paths++;
                        }
                    }
                }
            }
        });
    });

    if(normalized) {
        for (size_t i = 0; i < centrality_scores.size(); i++) {
            centrality_scores.at(i) /= number_of_shortest_paths;
        }
    }
    return centrality_scores;
}

std::vector<nodeweight> HypergraphExpansions::lineGraphBetweenness(Hypergraph &hypergraph, bool normalized) {
    
    Graph lineGraph = HypergraphExpansions::lineGraph(hypergraph, true);
    std::vector<nodeweight> centrality_scores(hypergraph.numberOfNodes());
    size_t number_of_shortest_paths = 0;

    //Get a map of all edges a node is part of
    std::map<node, std::set<edgeid>> edges;
    hypergraph.forNodes([&](node node){
        hypergraph.forEdges([&](edgeid eid){
            if(hypergraph.hasNode(node, eid)) {
                edges[node].insert(eid);
            }
        });
    });

    //calculate shortest paths on lineGraph and store them in a map
    std::map<std::pair<node,node>,std::pair<std::set<std::vector<node>>,edgeweight>> line_graph_shortest_paths;
    lineGraph.forNodes([&](node edge1){
        //std::cout << "Step " << edge1 << " / " << lineGraph.numberOfNodes() << std::endl;
        Dijkstra dijkstra(lineGraph, edge1, true, false);
        dijkstra.run();
        lineGraph.forNodes([&](node edge2){
            if(edge1 < edge2) {
                line_graph_shortest_paths[{edge1, edge2}] = {dijkstra.getPaths(edge2), dijkstra.distance(edge2)};
                line_graph_shortest_paths[{edge2, edge1}] = line_graph_shortest_paths[{edge1, edge2}];
            }
        });
    });


    //calculate centrality
    hypergraph.forNodes([&](node node1){
        hypergraph.forNodes([&](node node2){
            if(node1 != node2) {
                //variables for getting the shortest paths on the line graph
                double shortest_length = std::numeric_limits<double>::max();
                std::set<std::vector<node>> shortest_paths;
                edgeweight distance;
                //all edges node1 is part of
                for(auto edge1 : edges[node1]) {
                    //all edges node2 is part of
                    for(auto edge2 : edges[node2]) {
                        distance = line_graph_shortest_paths[{edge1, edge2}].second;
                        //search shortest path
                        if (distance < shortest_length) {
                                shortest_length = distance;
                                shortest_paths = line_graph_shortest_paths[{edge1, edge2}].first;
                        } else if (distance == shortest_length) {
                            for(auto path : line_graph_shortest_paths[{edge1, edge2}].first) {
                                shortest_paths.insert(path);
                            }
                        }
                    }
                }
                //Now that we have the shortest paths between node 1 and 2, do the betweenness scores
                for(std::vector<node> path : shortest_paths) {
                    for (size_t i = 0; i < path.size()-1; i++){
                        std::set<node> intersection = HypergraphExpansions::getIntersection(hypergraph, path.at(i), path.at(i+1));
                        //add 1 / intersection_size to centrality scores for each node in intersection and divide by number of paths if scores should be normalized
                        for(node node : intersection) {
                            centrality_scores.at(node) += 1.0 / double(intersection.size() * shortest_paths.size());
                            number_of_shortest_paths++;
                        }
                    }
                }
            }
        });
    });

    if(normalized) {
        for (size_t i = 0; i < centrality_scores.size(); i++) {
            centrality_scores.at(i) /= number_of_shortest_paths;
        }
    }
    return centrality_scores;
}

} //namespace NetworKit
