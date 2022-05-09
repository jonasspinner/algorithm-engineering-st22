#pragma once

#include <optional>

#include "definitions.hpp"

algen::WEdgeList fast_kruskal(algen::WEdgeList edges, std::size_t n);
algen::WEdgeList fast_kruskal(algen::WEdgeList edges);

std::vector<bool> are_edges_light(const algen::WEdgeList& edges, const algen::WEdgeList& spanning_tree, std::size_t n, bool perform_checks = false);
std::optional<algen::WEdge> verify_spanning_tree(const algen::WEdgeList& edges, const algen::WEdgeList& spanning_tree, std::size_t n, bool perform_checks = false);
