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

TEST_P(HypergraphExpansionsGTest, testGetIntersectionSize) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    EXPECT_EQ(HypergraphExpansions::getIntersectionSize(hGraph, 0, 1, 1), 1);
    EXPECT_EQ(HypergraphExpansions::getIntersectionSize(hGraph, 0, 1, 0), 0);
    EXPECT_EQ(HypergraphExpansions::getIntersectionSize(hGraph, 0, 2, 1), 0);
    EXPECT_EQ(HypergraphExpansions::getIntersectionSize(hGraph, 1, 2, 3), 2);
    EXPECT_EQ(HypergraphExpansions::getIntersectionSize(hGraph, 1, 2, 2), 2);
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

} //namespace NetworKit
