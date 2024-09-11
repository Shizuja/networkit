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

    EXPECT_EQ(hGraph.numberOfNodes(), 4);
    EXPECT_EQ(hGraph.numberOfEdges(), 3);
    
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

    ASSERT_EQ(hGraph.numberOfNodes(), 4);
    ASSERT_EQ(hGraph.numberOfEdges(), 3);
    /*
    auto output = HypergraphExpansions::lineExpansion(hGraph);
    Graph lGraph = output.first;
    auto nodeMap = output.second;

    ASSERT_EQ(lGraph.numberOfNodes(), 7);
    ASSERT_EQ(lGraph.numberOfEdges(), 8);

    for (auto entry : nodeMap) {
        node o_node = entry.second.first;
        edgeid o_edge = entry.second.second;
        ASSERT_TRUE(hGraph.hasNode(o_node, o_edge));
    }*/
}

TEST_P(HypergraphExpansionsGTest, testRecoverHypergraphFromLineExpansion) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({0, 3});

    ASSERT_EQ(hGraph.numberOfNodes(), 4);
    ASSERT_EQ(hGraph.numberOfEdges(), 3);
    /*
    auto output = HypergraphExpansions::lineExpansion(hGraph);

    ASSERT_EQ(output.first.numberOfNodes(), 7);
    ASSERT_EQ(output.first.numberOfEdges(), 8);

    Hypergraph rGraph = HypergraphExpansions::reconstructHypergraphFromLineExpansion(output.first, output.second);

    hGraph.forEdges([&](edgeid eid) {
        ASSERT_TRUE(rGraph.hasEdge(eid));
    });

    hGraph.forNodes([&](node node) {
        ASSERT_TRUE(rGraph.getNeighbors(node) == hGraph.getNeighbors(node));
    });*/
}

} //namespace NetworKit
