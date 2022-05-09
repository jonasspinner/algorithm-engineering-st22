#pragma once

#include "naive_jarnik_prim.hpp"
#include "naive_kruskal.hpp"

namespace mst_construction {
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
// The algorithm instance must be callable with an edge list as its argument,
// i.e. it must define
//      WEdgeList operator()(const WEdgeList& edge_list, const VertexId num_vertices)
//      {...}.
// The call operator serves as the main interface for computing the MST.
// It must accept an edge list (in arbitrary order) that represents the input
// graph. It should return an edge list (in any order) that represents the MST
// of the input graph. For the format requirements that any edge list has to
// fulfill, see edge_list_format_check() in includes/utils.hpp.

// Register your contenders in this tuple:
constexpr std::tuple contenders{
    // Some examples:

    // Slow kruskal, deactivate for larger graphs (log_m > 16)
    Contender{"naive_kruskal", [] { return NaiveKruskal(); }},

    // Faster kruskal from the binary library
    Contender{"fast_kruskal", [] {
        return [](const algen::WEdgeList& el, const algen::VertexId n) {
            return fast_kruskal(el, n);
        };
    }},

    // Jarnik-Prim with inefficient addressable PQ
    Contender{"naive_jarnik_prim", [] { return NaiveJarnikPrim(); }},

    // An example returning a badly formatted edge list
    Contender{"outputs_badly_formatted_edge_list", [] {
        return [](const algen::WEdgeList& el, const algen::VertexId) {
            algen::WEdgeList el_copy = el;
            // Duplicate edge so edge list will fail sanity check
            el_copy.push_back(el.back());
            return el_copy;
        };
    }},

    // An example returning a spanning tree that is not an MST
    Contender{"outputs_corrupted_mst", [] {
        return [](const algen::WEdgeList& el, const algen::VertexId n) {
            using namespace algen;
            auto call_fast_kruskal = [](const WEdgeList& edges, const VertexId n) { return fast_kruskal(edges, n);};
            benchmark::CorruptedMSTGenerator<decltype(call_fast_kruskal)> instance_gen;
            instance_gen.preprocess(el, n, call_fast_kruskal);
            return instance_gen.generate_corrupted_mst(1);
        };
    }}};
constexpr auto num_contenders = std::tuple_size_v<decltype(contenders)>;

constexpr std::size_t iterations = 1;

constexpr struct {
    std::size_t log_n = 16; // number of vertices = 2^(log_n)
    std::size_t log_m = 16; // number undirected edges = 2^(log_m)
    algen::Weight max_weight = 255; // maximum weight of generated edges
} graph_generator_params;

}  // end namespace params
}  // end namespace mst_construction
