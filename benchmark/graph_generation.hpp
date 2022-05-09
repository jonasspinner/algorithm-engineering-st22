#pragma once

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <random>
#include <unordered_set>

#include "includes/definitions.hpp"
#include "includes/utils.hpp"

namespace benchmark {

inline void print_log_degree_stats(const algen::WEdgeList& edges, std::size_t log_n) {
  if (edges.empty()) {
    return;
  }
  std::vector<std::size_t> num_buckets(log_n + 2, 0);
  std::size_t cur_degree = 0;
  algen::VertexId prev_tail = edges.front().tail;
  for (const auto& edge : edges) {
    if (edge.tail == prev_tail) {
      ++cur_degree;
    } else {
      ++num_buckets[std::log2(cur_degree)];
      cur_degree = 1;
      prev_tail = edge.tail;
    }
  }
  const auto num_non_isolated_edges =
      std::accumulate(num_buckets.begin(), num_buckets.end(), 0ull);
  std::cout << "#isolated vertices: "
            << (1ull << log_n) - num_non_isolated_edges << std::endl;
  for (std::size_t i = 0; i < num_buckets.size(); ++i) {
    std::size_t start = 1ull << i;
    std::size_t end = 1ull << (i + 1);
    std::cout << "[" << std::setw(5) << start << ", " << std::setw(5) << end
              << "): " << num_buckets[i] << std::endl;
  }
}

inline void print_degree_stats(const algen::WEdgeList& edges, std::size_t log_n) {
  if (edges.empty()) {
    return;
  }
  std::vector<std::size_t> num_buckets(50, 0);
  std::size_t cur_degree = 0;
  algen::VertexId prev_tail = edges.front().tail;
  for (const auto& edge : edges) {
    if (edge.tail == prev_tail) {
      ++cur_degree;
    } else {
      ++num_buckets[std::min(cur_degree, (num_buckets.size() - 1))];
      cur_degree = 1;
      prev_tail = edge.tail;
    }
  }
  const auto num_non_isolated_edges =
      std::accumulate(num_buckets.begin(), num_buckets.end(), 0ull);
  std::cout << "#isolated vertices: "
            << (1ull << log_n) - num_non_isolated_edges << std::endl;
  for (std::size_t i = 0; i < num_buckets.size(); ++i) {
    std::cout << "[" << std::setw(5) << i << "]: " << num_buckets[i]
              << std::endl;
  }
}

class GNM_Generator {
 public:
  void configure(std::size_t log_n, std::size_t log_m, algen::Weight max_weight = 255,
                 std::size_t seed = 1, bool show_degree_stats = false) {
    log_n_ = log_n;
    log_m_ = log_m;
    max_weight_ = max_weight;
    show_degree_stats_ = show_degree_stats;
    seed_ = seed;
  }

  // implementation follows approach from Batagelj, Vladimir, and Ulrik Brandes.
  // "Efficient generation of large random networks." Physical Review E 71.3
  // (2005): 036113.
  algen::WEdgeList generate() const {
    using namespace algen;
    if (log_n_ == 0 || log_m_ == 0) {
      return WEdgeList{};
    }
    const std::size_t n = 1ull << log_n_;
    const std::size_t m = 1ull << log_m_;
    const std::size_t edge_universe_size =
        (n * (n - 1) / 2) - 1;  // number of possible (undirected) edges
    std::mt19937 gen(seed_);
    std::uniform_int_distribution<std::size_t> uniform_edge_dist(
        0, edge_universe_size);
    std::uniform_int_distribution<Weight> uniform_weight_dist(1, max_weight_);
    WEdgeList edge_list;
    edge_list.reserve(2 * m);
    std::unordered_set<std::size_t> selected_edges;
    for (std::size_t i = 0; i < m; ++i) {
      std::size_t edge_id = -1;
      do {
        edge_id = uniform_edge_dist(gen);
      } while (selected_edges.count(edge_id) == 1);
      selected_edges.insert(edge_id);
    }
    std::cout << "finish sample edges" << std::endl;
    for (const auto& id : selected_edges) {
      const auto edge = get_tail_head_from_id(id);
      const Weight w = uniform_weight_dist(gen);
      edge_list.emplace_back(edge.tail, edge.head, w);
      edge_list.emplace_back(edge.head, edge.tail, w);
    }
    postprocessing_and_checks(edge_list, m);
    std::cout << "finish get names" << std::endl;
    if (show_degree_stats_) {
      print_log_degree_stats(edge_list, log_n_);
    }
    return edge_list;
  }

 private:
  std::size_t log_n_ = 8;
  std::size_t log_m_ = 10;
  algen::Weight max_weight_ = 255;
  std::size_t seed_ = 1;
  bool show_degree_stats_ = false;

  algen::UnweightedEdge get_tail_head_from_id(std::size_t i) const {
    using namespace algen;
    VertexId tail = 1 + std::floor(-0.5 + std::sqrt(0.25 + 2 * i));
    VertexId head = i - (tail * (tail - 1)) / 2;
    return UnweightedEdge{tail, head};
  }
  void postprocessing_and_checks(algen::WEdgeList& edges, const std::size_t m) const {
    using namespace algen;
    std::sort(edges.begin(), edges.end(), TailHeadOrder<WEdge>{});
    auto TailHeadEqual = [](const WEdge& lhs, const WEdge& rhs) {
      return std::tie(lhs.tail, lhs.head) == std::tie(rhs.tail, rhs.head);
    };
    // following lines are just a sanity check
    auto it = std::unique(edges.begin(), edges.end(), TailHeadEqual);
    const auto num_duplicates = std::distance(it, edges.end());
    edges.erase(it, edges.end());
    assert(num_duplicates == 0);
    assert(edges.size() == 2 * m);  // back edges
  }
};

}  // namespace benchmark
