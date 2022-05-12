#pragma once

#include <optional>

#include "definitions.hpp"

// Fast kruskal implementation. This is used by the framework and as a baseline contender.
// You are not allowed to call this method in your implementation.
algen::WEdgeList fast_kruskal(algen::WEdgeList edges, std::size_t n);
algen::WEdgeList fast_kruskal(algen::WEdgeList edges);

// Given the edge list edges of an entire graph, a spanning tree spanning_tree of that graph,
// and the number of vertices n in the graph, this method checks which graph edges are light
// wrt to spanning_tree. It returns a bool vector bv of size m:=edges.size()
// s.t. bv[i] iff edges[i] is light wrt spanning_tree (i = 0, .., m-1).
// You are not allowed to call this method in your implementation.
// If you work on "Aufgabe 4: Linear-Time MST Algorithmus mit Blackbox" use
// algen::getEdgeClassifier() as your black box instead.
std::vector<bool> are_edges_light(const algen::WEdgeList& edges, const algen::WEdgeList& spanning_tree, std::size_t n, bool perform_checks = false);

// Given the edge list edges of an entire graph, a spanning tree spanning_tree of that graph,
// and the number of vertices n in the graph, this method checks whether spanning_tree is an
// MST of the graph. If yes, an empty optional is returned. If no, the returned optional
// contains a graph edge that is light wrt to spanning_tree.
// You are not allowed to call this method in your implementation.
std::optional<algen::WEdge> verify_spanning_tree(const algen::WEdgeList& edges, const algen::WEdgeList& spanning_tree, std::size_t n, bool perform_checks = false);
