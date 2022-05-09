#pragma once
#include "includes/definitions.hpp"
#include "naive_dfs_based_verifier.hpp"
#include "includes/binary_includes.hpp"
#include <optional>

namespace mst_verification {
    namespace params {
        template <class F>
        struct Contender {
            std::string_view name;
            F factory;
        };
        template <class F>
        Contender(std::string_view, F) -> Contender<F>;

// List of contenders.
// Each contender must have a unique name, and provide
// a factory for constructing an instance of the algorithm.
// The factory must be callable without parameters.
// The algorithm instance must be callable with two edge lists and a number
// of vertices as arguments, the first edge list being the entire graph and
// the second one being an (M)ST to be verified. The second edge list can
// be assumed to be a spanning tree (or forest) of the given graph.
// It should return a std::optional<WEdge> s.t. if the given ST T is an MST,
// the optional is empty and if T is not an MST the optional contains an
// edge e s.t. e is a graph edge that is T-light.
// In summary, the algorithm must define
//      std::optional<WEdge> operator()(const WEdgeList& graph_edge_list,
//                                        const WEdgeList& st_edge_list,
//                                        const int num_vertices) {...}.
// For the format requirements that any edge list has to fulfill, see
// edge_list_format_check() in includes/utils.hpp.

// Register your contenders in this tuple:
        constexpr std::tuple contenders{
                    // Some examples:

                    // A naive MST verifier that executes one DFS in the given ST per graph edge.
                    Contender{"naive_dfs_verify", [] { return NaiveDFSBasedVerifier(); }},

                    // An MST verifier that is implemented in the binary library.
                    Contender{"fast_verify_from_binary", [] {
                        return [](const algen::WEdgeList& edges, const algen::WEdgeList& st_edges, const int num_vertices) {
                            return verify_spanning_tree(edges, st_edges, num_vertices);
                        };
                    }},

                    // An example contender that always outputs yes (not very accurate).
                    Contender{"always_outputs_yes", [] {
                        return [](const algen::WEdgeList&, const algen::WEdgeList&, const int) {
                            return std::optional<algen::WEdge>();
                        };
                    }},

                    // An example contender that always outputs no and outputs an ST edge as a counter example
                    // (which is always an invalid counter example).
                    Contender{"always_outputs_no_with_bad_counter_example", [] {
                        return [](const algen::WEdgeList&, const algen::WEdgeList& st_edges, const int) {
                            return std::optional<algen::WEdge>(st_edges.front());
                        };
                    }}
                };
        constexpr auto num_contenders = std::tuple_size_v<decltype(contenders)>;

        constexpr std::size_t iterations = 1;

        constexpr struct {
            std::size_t log_n = 16; // number of vertices = 2^(log_n)
            std::size_t log_m = 16; // number undirected edges = 2^(log_m)
            algen::Weight max_weight = 255; // maximum weight of generated edges
        } graph_generator_params;

        constexpr std::size_t num_changed_edges = 100; // ST instance for verification algorithms is constructed by changing num_changed_edges in a valid MST

    }  // end namespace params
}
