#pragma once

#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>
#include <stack>

#include "definitions.hpp"

namespace algen {
template <typename EdgeType>
struct TailHeadOrder {
  bool operator()(const EdgeType& lhs, const EdgeType& rhs) const {
    return std::tie(lhs.tail, lhs.head) < std::tie(rhs.tail, rhs.head);
  }
};

template <typename EdgeType>
struct TailHeadWeightOrder {
  bool operator()(const EdgeType& lhs, const EdgeType& rhs) const {
    return std::tie(lhs.tail, lhs.head, lhs.weight) <
           std::tie(rhs.tail, rhs.head, rhs.weight);
  }
};

template <typename EdgeType>
struct WeightOrder {
  bool operator()(const EdgeType& lhs, const EdgeType& rhs) const {
    return lhs.weight < rhs.weight;
  }
};

template <typename Container>
void print_container(const Container& container) {
  for (const auto& elem : container) {
    std::cout << elem << std::endl;
  }
}

inline Weight sum_weights(const WEdgeList& edges) {
  Weight w = 0;
  for (const auto& [tail, head, weight] : edges) {
    w += weight;
  }
  return w;
}

inline Weight sum_weights(const WeightArray& edges) {
  Weight w = 0;
  for (const auto& weight : edges) {
    w += weight;
  }
  return w;
}

inline void add_back_edges(WEdgeList& edges) {
  WEdgeList back_edges;
  for (const auto& [tail, head, w] : edges) {
    back_edges.emplace_back(head, tail, w);
  }
  edges.insert(edges.end(), back_edges.begin(), back_edges.end());
  std::sort(edges.begin(), edges.end(), TailHeadOrder<WEdge>{});
  auto TailHeadEqual = [](const WEdge& lhs, const WEdge& rhs) {
    return std::tie(lhs.tail, lhs.head) == std::tie(rhs.tail, rhs.head);
  };
  auto it = std::unique(edges.begin(), edges.end(), TailHeadEqual);
  edges.erase(it, edges.end());
}

// Reorders edges and fills first_out_edge s.t. for a vertex v all edges e=(v,.)
// are in the range edges[first_out_edge[v]..first_out_edge[v+1]].
// More efficient representations regarding space and construction time are possible.
// This representation is only used for non-optimized framework methods.
void make_inefficient_adjacency_structure(algen::WEdgeList& edges,
                                          std::vector<EdgeIdx>& first_out_edge,
                                          const algen::VertexId num_vertices) {
    using namespace algen;

    std::sort(edges.begin(), edges.end(), TailHeadOrder<WEdge>{});
    first_out_edge.clear();
    VertexId cur_vertex = 0;
    first_out_edge.push_back(0);
    for (VertexId i = 0; i < edges.size(); ++i) {
        if (edges[i].tail != cur_vertex) {
            first_out_edge.insert(first_out_edge.end(), edges[i].tail - cur_vertex, i);
            cur_vertex = edges[i].tail;
            assert(first_out_edge.size() == cur_vertex + 1);
        }
    }
    first_out_edge.insert(first_out_edge.end(), num_vertices - first_out_edge.size() + 1, edges.size());
    assert(first_out_edge.size() == num_vertices + 1 && first_out_edge[num_vertices] == edges.size());
}

// Edge Lists have to comply with the following requirements:
// 1.   The edge list is not empty.
// 2.   No duplicate edges occur in the edge list, i.e. for each pair of vertex
//      ids u,v at most one edge e={u,v,w} and one edge e'={v,u,w} exist.
// 3.   Every edge is present in both directions, i.e. for each edge e={u,v,w}
//      in the list, e'={v,u,w} is also contained in the list.
// 4.   The vertex ids occurring in the edge list are smaller than the given
//      number of vertices.
    std::pair<bool, std::string> edge_list_format_check(
            const algen::WEdgeList& edges, const algen::VertexId num_vertices) {
        using namespace algen;
        std::stringstream msg;

        // Check that list is not empty.
        if (edges.empty()) {
            msg << "Edge list is empty!";
            return {false, msg.str()};
        }

        // Make sure the vertex ids are smaller than num_vertices
        for (const auto& edge : edges) {
            VertexId larger = std::max(edge.tail, edge.head);
            if (larger >= num_vertices) {
                msg << "Vertex id " << larger << " in edge (" << edge.tail << " - "
                    << edge.head << ") is invalid. (n = " << num_vertices << ")."
                    << std::endl;
                return {false, msg.str()};
            }
        }

        // Make adjacency structure
        auto sorted_edges = edges;
        std::vector<EdgeIdx> first_out_edge;
        make_inefficient_adjacency_structure(sorted_edges, first_out_edge, num_vertices);

        // Check that there are no duplicate edges
        for (VertexId i = 1; i < sorted_edges.size(); ++i) {
            if (sorted_edges[i - 1].head == sorted_edges[i].head &&
                    sorted_edges[i - 1].tail == sorted_edges[i].tail) {
                msg << "Edge (" << sorted_edges[i].tail << " - " << sorted_edges[i].head
                    << ") exists more than once!";
                return {false, msg.str()};
            }
        }

        // todo this can surely be done faster (n x n bit matrix is a lot of memory
        //  but maybe hash instead?)

        // For each edge check if the head also has an outgoing edge to the tail.
        std::vector<bool> checked_already(sorted_edges.size(), false);
        for (EdgeIdx i = 0; i < sorted_edges.size(); ++i) {
            const auto edge = sorted_edges[i];
            if (edge.head == edge.tail || checked_already[i]) continue;
            bool found_tail = false;
            for (EdgeIdx j = first_out_edge[edge.head];
                 j < first_out_edge[edge.head + 1]; ++j) {
                if (sorted_edges[j].head == edge.tail) {
                    // Check that reverse edge has same weight
                    if (sorted_edges[j].weight != edge.weight) {
                        msg << "Edge (" << edge.tail << " - " << edge.head << ") has weight "
                            << edge.weight << " but edge (" << sorted_edges[j].tail << " - "
                            << sorted_edges[j].head << ") has weight " << sorted_edges[j].weight << "!";
                        return {false, msg.str()};
                    }

                    found_tail = true;
                    checked_already[j] = true;
                    break;
                }
            }
            if (!found_tail) {
                msg << "Could not find reverse edge for (" << edge.tail << " - "
                    << edge.head << ")!";
                return {false, msg.str()};
            }
        }

        return {true, "Edge list has correct format."};
    }

    // Returns {true, "Is spanning tree."} if the given tree edges are a spanning forest of the given graph edges.
    // Returns {false, msg} otherwise where msg describes the reason (not spanning or not a forest).
    std::pair<bool, std::string> is_spanning_forest(const WEdgeList& graph_edges,
                                                  const WEdgeList& tree_edges,
                                                  const VertexId num_vertices) {

        auto sorted_graph_edges = graph_edges;
        std::vector<EdgeIdx> first_graph_out_edge;
        make_inefficient_adjacency_structure(sorted_graph_edges, first_graph_out_edge, num_vertices);

        auto sorted_tree_edges = tree_edges;
        std::vector<EdgeIdx> first_tree_out_edge;
        make_inefficient_adjacency_structure(sorted_tree_edges, first_tree_out_edge, num_vertices);


        std::vector<bool> vertex_seen(num_vertices, false);
        std::vector<EdgeIdx> next_out_edge(num_vertices, 0);
        std::vector<VertexId> parent(num_vertices, VERTEXID_UNDEFINED);
        std::stack<VertexId, std::vector<VertexId>> active;

        // Count connected components using a DFS. Also sets a flag indicating whether the given graph has any cycles.
        auto count_ccs = [&] (const WEdgeList& edges, const std::vector<EdgeIdx>& first_out_edge, bool& has_cycle) {

            assert(active.empty());
            VertexId cc_counter = 0;
            has_cycle = false;

            for (VertexId s = 0; s < num_vertices; ++s) {
                if (vertex_seen[s]) continue;
                ++cc_counter;
                vertex_seen[s] = true;
                active.push(s);
                parent[s] = s;
                while (!active.empty()) {
                    const auto v = active.top();
                    const EdgeIdx out_edge_idx = first_out_edge[v] + next_out_edge[v];
                    if (out_edge_idx < first_out_edge[v + 1]) {
                        ++next_out_edge[v];
                        const auto head = edges[out_edge_idx].head;
                        if (!vertex_seen[head]) {
                            vertex_seen[head] = true;
                            active.push(head);
                            parent[head] = v;
                        } else if (head != parent[v]) {
                            has_cycle = true;
                        }
                    } else {
                        active.pop();
                    }
                }
            }

            return cc_counter;
        };

        bool tree_has_cycle;
        const auto num_tree_css = count_ccs(sorted_tree_edges, first_tree_out_edge, tree_has_cycle);
        if (tree_has_cycle) {
            return {false, "Given tree edges contain cycle."};
        }

        std::fill(vertex_seen.begin(), vertex_seen.end(), false);
        std::fill(next_out_edge.begin(), next_out_edge.end(), 0);
        std::fill(parent.begin(), parent.end(),VERTEXID_UNDEFINED);

        bool graph_has_cycle;
        const auto num_graph_ccs = count_ccs(sorted_graph_edges, first_graph_out_edge, graph_has_cycle);
        if (num_tree_css != num_graph_ccs) {
            return {false, "Given tree is not spanning."};
        }

        return {true, "Is spanning tree."};
    }


}  // namespace algen
