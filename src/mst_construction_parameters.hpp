#pragma once

#include "benchmark/verification_instance_generation.hpp"
#include "includes/binary_includes.hpp"
#include "naive_jarnik_prim.hpp"
#include "naive_kruskal.hpp"
#include "expected_linear_time_mst.hpp"

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
//      WEdgeList operator()(const WEdgeList& edge_list, const VertexId
//      num_vertices)
//      {...}.
// The call operator serves as the main interface for computing the MST.
// It must accept an edge list (in arbitrary order) that represents the input
// graph. It should return an edge list (in any order) that represents the MST
// of the input graph. For the format requirements that any edge list has to
// fulfill, see edge_list_format_check() in includes/utils.hpp.

// Register your contenders in this tuple:
constexpr std::tuple contenders{
    Contender{"expected_linear_time_mst",
              [] { return js::ExpectedLinearTimeMST(); }},

    // Some examples:

    // Slow kruskal, deactivate for larger graphs (log_m > 16)
    // Contender{"naive_kruskal", [] { return NaiveKruskal(); }},

    // Faster kruskal from the binary library
    Contender{"fast_kruskal",
              [] {
                return [](const algen::WEdgeList& el, const algen::VertexId n) {
                  return fast_kruskal(el, n);
                };
              }},

    // Jarnik-Prim with inefficient addressable PQ
//    Contender{"naive_jarnik_prim", [] { return NaiveJarnikPrim(); }},

    // An example returning a badly formatted edge list
    // Contender{"outputs_badly_formatted_edge_list",
    //          [] {
    //            return [](const algen::WEdgeList& el, const algen::VertexId) {
    //              algen::WEdgeList el_copy = el;
    //              // Duplicate edge so edge list will fail sanity check
    //              el_copy.push_back(el.back());
    //              return el_copy;
    //            };
    //          }},

    // An example returning a spanning tree that is not an MST
//     Contender{"outputs_corrupted_mst", [] {
//                return [](const algen::WEdgeList& el, const algen::VertexId n)
//                {
//                  using namespace algen;
//                  auto correct_mst = fast_kruskal(el, n);
//                  benchmark::CorruptedMSTGenerator instance_gen(1337);
//                  instance_gen.preprocess(el, correct_mst, n);
//                  return instance_gen.generate_corrupted_mst(1);
//                };
//              }}
};
constexpr auto num_contenders = std::tuple_size_v<decltype(contenders)>;

constexpr std::size_t iterations = 4;
struct Experiment {
  std::size_t log_n;
  std::size_t edge_factor;
  algen::Weight max_weight;
  friend std::ostream& operator<<(std::ostream& out, const Experiment& exp) {
    return out << "Experiment Config = (" << exp.log_n << ", "
               << exp.edge_factor << ", " << exp.max_weight << ")";
  }
};

struct ExperimentSuite {
  std::size_t log_n_begin = 14;
  std::size_t log_n_end = 19;
  std::size_t edge_factor_begin = 1;
  std::size_t edge_factor_end = 256;
  algen::Weight max_weight = 255;
  std::size_t step_size_n = 1;
  std::size_t step_size_edge_factor = 2;
  std::size_t cur_log_n = log_n_begin;
  std::size_t cur_edge_factor = edge_factor_begin;
  bool has_next() const { return cur_log_n <= log_n_end; }
  Experiment get_next() {
    Experiment exp{cur_log_n, cur_edge_factor, max_weight};
    cur_edge_factor *= step_size_edge_factor;
    if (cur_edge_factor > edge_factor_end) {
      cur_edge_factor = edge_factor_begin;
      cur_log_n += step_size_n;
    }
    return exp;
  }
};

}  // end namespace params
}  // end namespace mst_construction
