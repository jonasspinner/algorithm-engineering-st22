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
//                    Contender{"naive_dfs_verify", [] { return NaiveDFSBasedVerifier(); }},

                    // An MST verifier that is implemented in the binary library.
                    Contender{"fast_verify_from_binary", [] {
                        return [](const algen::WEdgeList& edges, const algen::WEdgeList& st_edges, const int num_vertices) {
                            return verify_spanning_tree(edges, st_edges, num_vertices);
                        };
                    }},

                    // An example contender that always outputs yes (not very accurate).
//                    Contender{"always_outputs_yes", [] {
//                        return [](const algen::WEdgeList&, const algen::WEdgeList&, const int) {
//                            return std::optional<algen::WEdge>();
//                        };
//                    }},

                    // An example contender that always outputs no and outputs an ST edge as a counter example
                    // (which is always an invalid counter example).
//                    Contender{"always_outputs_no_with_bad_counter_example", [] {
//                        return [](const algen::WEdgeList&, const algen::WEdgeList& st_edges, const int) {
//                            return std::optional<algen::WEdge>(st_edges.front());
//                        };
//                    }}
                };
        constexpr auto num_contenders = std::tuple_size_v<decltype(contenders)>;

        constexpr std::size_t iterations = 10;

        struct Experiment {
            std::size_t log_n;
            std::size_t edge_factor;
            std::size_t num_changed_edges;
            algen::Weight max_weight;
            bool generateNewGraph;
            friend std::ostream& operator<<(std::ostream& out, const Experiment& exp) {
                return out << "Experiment Config = ("
                           << exp.log_n << ", "
                           << exp.edge_factor << ", "
                           << exp.num_changed_edges << ", "
                           << exp.max_weight << ")";
            }
        };

        struct ExperimentSuite {
            std::size_t log_n_begin = 14;
            std::size_t log_n_end = 19;
            std::size_t edge_factor_begin = 1;
            std::size_t edge_factor_end = 256;
            bool try_yes_instances = true;
            std::size_t num_changed_edges_begin = 1;
            std::size_t num_changed_edges_end = 1000;
            algen::Weight max_weight = 255;
            std::size_t step_size_log_n = 1;
            std::size_t step_size_edge_factor = 2;
            std::size_t step_size_num_changed_edges = 10;

            std::size_t cur_log_n = log_n_begin; // Generate a graph with 2^cur_log_n vertices
            std::size_t cur_edge_factor = edge_factor_begin; // Generate a graph with cur_edge_factor times as many edges as vertices
            std::size_t cur_num_changed_edges = try_yes_instances? 0 : num_changed_edges_begin; // Construct an ST instance by changing cur_num_changed_edges in a valid MST

            bool has_next() const { return cur_log_n <= log_n_end; }
            Experiment get_next() {
                const bool need_new_graph = (try_yes_instances && cur_num_changed_edges == 0) || (!try_yes_instances && cur_num_changed_edges == num_changed_edges_begin);
                Experiment exp{cur_log_n, cur_edge_factor, cur_num_changed_edges, max_weight, need_new_graph};
                cur_num_changed_edges = (try_yes_instances && cur_num_changed_edges == 0)? num_changed_edges_begin : cur_num_changed_edges * step_size_num_changed_edges;
                if (cur_num_changed_edges > num_changed_edges_end) {
                    cur_num_changed_edges = try_yes_instances? 0 : num_changed_edges_begin;
                    cur_edge_factor *= step_size_edge_factor;
                    if (cur_edge_factor > edge_factor_end) {
                        cur_edge_factor = edge_factor_begin;
                        cur_log_n += step_size_log_n;
                    }
                }
                return exp;
            }
        };


    }  // end namespace params
}
