#pragma once

#include <optional>
#include <stack>
#include <cassert>
#include <algorithm>
#include "includes/utils.hpp"

// Naive MST verifier. Runs one DFS on the given ST T per graph edge e to find
// the heaviest edge in the st-path from e.head to e.tail to determine
// if e is T-light.
class NaiveDFSBasedVerifier {

    using EdgeIdx = std::uint64_t;

public:
    NaiveDFSBasedVerifier() = default;

    std::optional<algen::WEdge> operator()(const algen::WEdgeList& graph_edges, const algen::WEdgeList& st_edges,
            const algen::VertexId num_vertices) {

        auto sorted_st_edges = st_edges;
        std::vector<EdgeIdx> first_st_out_edge;
        make_inefficient_adjacency_structure(sorted_st_edges, first_st_out_edge, num_vertices);

        for (const auto& e : graph_edges) {
            if (is_st_light(e, sorted_st_edges, first_st_out_edge))
                return e;
        }
        return {};
    }

    static bool is_st_light(const algen::WEdge& e, const algen::WEdgeList& st_edges, const std::vector<EdgeIdx>& first_st_out_edge) {
        using namespace algen;
        const VertexId num_vertices = first_st_out_edge.size() - 1;

        std::vector<bool> marked(num_vertices);
        std::vector<Weight> heaviest_on_path(num_vertices, 0);
        std::vector<EdgeIdx> next_st_edge_offset(num_vertices, 0);
        std::stack<VertexId, std::vector<VertexId>> active;

        // Find the weight of the heaviest edge on the st-path from e.head to e.tail.
        active.push(e.head);
        marked[e.head] = true;
        while (!active.empty()) {
            const auto v = active.top();

            const auto next_st_edge = first_st_out_edge[v] + next_st_edge_offset[v];
            if (next_st_edge != first_st_out_edge[v + 1]) {
                const VertexId head = st_edges[next_st_edge].head;
                ++next_st_edge_offset[v];
                if (!marked[head]) {
                    heaviest_on_path[head] = std::max(heaviest_on_path[v], st_edges[next_st_edge].weight);
                    if (head == e.tail)
                        break;

                    active.push(head);
                    marked[head] = true;
                }
            } else {
                active.pop();
            }
        }
        assert(!active.empty());

        // e is st-light if e is lighter than the heaviest edge on the st-path from e.tail to e.head.
        return heaviest_on_path[e.tail] > e.weight;
    }

};
