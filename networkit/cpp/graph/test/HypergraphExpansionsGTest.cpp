/*
 * HypergraphExpansionsGTest.cpp
 *
 *  Created on: 11.09.2024
 *      Author: Yanneck Dimitrov
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/HypergraphExpansions.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>
#include <map>
#include <set>
#include <vector>
#include <iostream>

namespace NetworKit {

class HypergraphExpansionsGTest : public testing::TestWithParam<std::tuple<bool>> {
protected:
    bool isHypergraph() const { return !isWeighted(); }
    bool isWeightedHypergraph() const { return isWeighted(); }

    bool isWeighted() const;
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, HypergraphExpansionsGTest,
                         testing::Values(std::make_tuple(false), std::make_tuple(true)));

bool HypergraphExpansionsGTest::isWeighted() const {
    return std::get<0>(GetParam());
}

/** GRAPH EXPANSIONS **/

TEST_P(HypergraphExpansionsGTest, testCliqueExpansion) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({0, 3});
    
    Graph cGraph = HypergraphExpansions::cliqueExpansion(hGraph);

    EXPECT_EQ(cGraph.numberOfNodes(), 4);
    EXPECT_EQ(cGraph.numberOfEdges(), 5);
    ASSERT_TRUE(cGraph.hasEdge(0,1));
    ASSERT_TRUE(cGraph.hasEdge(1,2));
    ASSERT_TRUE(cGraph.hasEdge(1,3));
    ASSERT_TRUE(cGraph.hasEdge(2,3));
    ASSERT_TRUE(cGraph.hasEdge(0,3));
}

TEST_P(HypergraphExpansionsGTest, testLineExpansion) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({0, 3});
    
    auto output = HypergraphExpansions::lineExpansion(hGraph);
    Graph lGraph = output.first;
    auto nodeMap = output.second;

    EXPECT_EQ(lGraph.numberOfNodes(), 7);
    EXPECT_EQ(lGraph.numberOfEdges(), 8);

    for (auto entry : nodeMap) {
        node o_node = entry.second.first;
        edgeid o_edge = entry.second.second;
        ASSERT_TRUE(hGraph.hasNode(o_node, o_edge));
    }
}

TEST_P(HypergraphExpansionsGTest, testRecoverHypergraphFromLineExpansion) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({0, 3});
    
    auto output = HypergraphExpansions::lineExpansion(hGraph);

    EXPECT_EQ(output.first.numberOfNodes(), 7);
    EXPECT_EQ(output.first.numberOfEdges(), 8);

    Hypergraph rGraph = HypergraphExpansions::reconstructHypergraphFromLineExpansion(output.first, output.second);

    hGraph.forEdges([&](edgeid eid) {
        ASSERT_TRUE(rGraph.hasEdge(eid));
    });

    ASSERT_TRUE(rGraph.hasNode(0,0));
    ASSERT_TRUE(rGraph.hasNode(1,0));
    ASSERT_TRUE(rGraph.hasNode(1,1));
    ASSERT_TRUE(rGraph.hasNode(2,1));
    ASSERT_TRUE(rGraph.hasNode(3,1));
    ASSERT_TRUE(rGraph.hasNode(0,2));
    ASSERT_TRUE(rGraph.hasNode(3,2));
    ASSERT_FALSE(rGraph.hasNode(2,0));
}

TEST_P(HypergraphExpansionsGTest, testGetIntersection) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    std::set<node> intersection0_1{1};
    std::set<node> intersection1_2{2, 3};
    std::set<node> intersection2_0{};

    EXPECT_EQ(HypergraphExpansions::getIntersection(hGraph, 0, 1), intersection0_1);
    EXPECT_EQ(HypergraphExpansions::getIntersection(hGraph, 1, 2), intersection1_2);
    EXPECT_EQ(HypergraphExpansions::getIntersection(hGraph, 2, 0), intersection2_0);
}

TEST_P(HypergraphExpansionsGTest, testLineExpansionWeightedBetweenness) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    std::pair<Graph, std::map<node, std::pair<node, edgeid>>> lineExpansion = HypergraphExpansions::lineExpansion(hGraph);
    auto scores = HypergraphExpansions::lineExpansionWeightedBetweenness(lineExpansion.first, lineExpansion.second, false);

    EXPECT_EQ(scores.at(0), 0);
    EXPECT_EQ(scores.at(1), 26);
    EXPECT_EQ(scores.at(2), 4);
    EXPECT_EQ(scores.at(3), 4);

    scores = HypergraphExpansions::lineExpansionWeightedBetweenness(lineExpansion.first, lineExpansion.second, true);

    EXPECT_EQ(scores.at(0), 0);
    EXPECT_EQ(scores.at(1), 13.0/15.0);
    EXPECT_EQ(scores.at(2), 2.0/15.0);
    EXPECT_EQ(scores.at(3), 2.0/15.0);
}

TEST_P(HypergraphExpansionsGTest, testLineExpansionBetweenness) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    std::pair<Graph, std::map<node, std::pair<node, edgeid>>> lineExpansion = HypergraphExpansions::lineExpansion(hGraph);
    auto scores = HypergraphExpansions::lineExpansionBetweenness(lineExpansion.first, lineExpansion.second, false, true);

    EXPECT_EQ(scores.at(0), 0);
    EXPECT_EQ(scores.at(1), 26);
    EXPECT_EQ(scores.at(2), 8);
    EXPECT_EQ(scores.at(3), 8);

    scores = HypergraphExpansions::lineExpansionBetweenness(lineExpansion.first, lineExpansion.second, true, true);

    EXPECT_EQ(scores.at(0), 0);
    EXPECT_EQ(scores.at(1), 13.0/15.0);
    EXPECT_EQ(scores.at(2), 4.0/15.0);
    EXPECT_EQ(scores.at(3), 4.0/15.0);

    scores = HypergraphExpansions::lineExpansionBetweenness(lineExpansion.first, lineExpansion.second, false, false);

    EXPECT_EQ(scores.at(0), 0);
    EXPECT_EQ(scores.at(1), 13);
    EXPECT_EQ(scores.at(2), 4);
    EXPECT_EQ(scores.at(3), 4);

    scores = HypergraphExpansions::lineExpansionBetweenness(lineExpansion.first, lineExpansion.second, true, false);

    EXPECT_EQ(scores.at(0), 0);
    EXPECT_EQ(scores.at(1), 6.5/15.0);
    EXPECT_EQ(scores.at(2), 2.0/15.0);
    EXPECT_EQ(scores.at(3), 2.0/15.0);
}

TEST_P(HypergraphExpansionsGTest, testNumberOfNodesFromNodeMap) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    std::pair<Graph, std::map<node, std::pair<node, edgeid>>> lineExpansion = HypergraphExpansions::lineExpansion(hGraph);

    EXPECT_EQ(HypergraphExpansions::numberOfNodesFromNodeMap(lineExpansion.second), 4);
}

TEST_P(HypergraphExpansionsGTest, testMemberOfHyperedges) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    std::pair<Graph, std::map<node, std::pair<node, edgeid>>> lineExpansion = HypergraphExpansions::lineExpansion(hGraph);
    std::vector<NetworKit::node> values = HypergraphExpansions::memberOfHyperedges(lineExpansion.second);

    EXPECT_EQ(values[0],1);
    EXPECT_EQ(values[1],2);
    EXPECT_EQ(values[2],2);
    EXPECT_EQ(values[2],2);
}

TEST_P(HypergraphExpansionsGTest, testGetEdgeMembers) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    std::set<node> edge0{0, 1};
    std::set<node> edge1{1, 2, 3};
    std::set<node> edge2{2, 3};

    EXPECT_EQ(HypergraphExpansions::getEdgeMembers(hGraph, 0),edge0);
    EXPECT_EQ(HypergraphExpansions::getEdgeMembers(hGraph, 1),edge1);
    EXPECT_EQ(HypergraphExpansions::getEdgeMembers(hGraph, 2),edge2);
}

TEST_P(HypergraphExpansionsGTest, testGetAllEdgeMembers) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    std::set<node> edge0{0, 1};
    std::set<node> edge1{1, 2, 3};
    std::set<node> edge2{2, 3};

    std::map<edgeid, std::set<node>> edge_members = HypergraphExpansions::getAllEdgeMembers(hGraph);

    EXPECT_EQ(edge_members[0],edge0);
    EXPECT_EQ(edge_members[1],edge1);
    EXPECT_EQ(edge_members[2],edge2);
}

TEST_P(HypergraphExpansionsGTest, testLineGraph) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    Graph weightedLineGraph = HypergraphExpansions::lineGraph(hGraph, true);
    EXPECT_TRUE(weightedLineGraph.hasEdge(1,2));
    EXPECT_EQ(weightedLineGraph.weight(1, 2), 0.5);
    EXPECT_FALSE(weightedLineGraph.hasEdge(2,0));
    EXPECT_EQ(weightedLineGraph.weight(2, 0), 0.0);


    Graph unweightedLineGraph = HypergraphExpansions::lineGraph(hGraph, false);
    EXPECT_TRUE(unweightedLineGraph.hasEdge(0,1));
    EXPECT_TRUE(unweightedLineGraph.hasEdge(1,2));
    EXPECT_FALSE(unweightedLineGraph.hasEdge(2,0));
}

TEST_P(HypergraphExpansionsGTest, testWeight) {
    Graph graph(3, true);
    graph.addEdge(0, 1, 2.5);
    graph.addEdge(0, 2, 0.5);

    EXPECT_TRUE(graph.isWeighted());
    EXPECT_EQ(graph.weight(0, 1), 2.5);
    EXPECT_EQ(graph.weight(0, 2), 0.5);
}

} //namespace NetworKit
