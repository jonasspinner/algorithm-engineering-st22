#pragma once

#include <queue>
#include <random>
#include "includes/definitions.hpp"

// This implementation of the JP algorithm is deliberately _very_ bad. It is
// only here to be used as a naive contender, showing how to use the contender
// interface.
class NaiveJarnikPrim {
  // An addressable priority queue that reinserts elements upon decrease_key()
  struct NaiveJPPQ {
    NaiveJPPQ(const algen::VertexId num_vertices)
        : contains(num_vertices, false),
          popped_before(num_vertices, false),
          pq() {}

    bool empty() {
      while (!pq.empty() && popped_before[pq.top().second]) {
        pq.pop();
      }
      return pq.empty();
    }

    void insert(const algen::Weight& weight, const algen::VertexId& vertex) {
      assert(!contains[vertex] && !popped_before[vertex]);
      pq.push({weight, vertex});
      contains[vertex] = true;
    }

    void decrease_key(const algen::Weight& new_weight,
                      const algen::VertexId& vertex) {
      assert(!popped_before[vertex] && contains[vertex]);
      pq.push({new_weight, vertex});
    }

    algen::WeightVertex delete_and_return_min() {
      while (popped_before[pq.top().second]) {
        pq.pop();
      }
      assert(!pq.empty());
      const auto res = pq.top();
      pq.pop();
      contains[res.second] = false;
      popped_before[res.second] = true;
      return res;
    }

   private:
    std::vector<bool> contains;
    std::vector<bool> popped_before;

    struct Compare {
      // Compare has to be in reverse order for std::priority_queue
      bool operator()(const algen::WeightVertex& wv1,
                      const algen::WeightVertex& wv2) {
        return wv1.first >= wv2.first;
      }
    };

    std::priority_queue<algen::WeightVertex,
                        std::vector<algen::WeightVertex>, Compare> pq;
  };

 public:
  algen::WEdgeList operator()(const algen::WEdgeList& edge_list,
                              const algen::VertexId num_vertices) {
    using namespace algen;

    static constexpr Weight INF_WEIGHT = std::numeric_limits<Weight>::max();
    static constexpr VertexId INVALID_ID = std::numeric_limits<VertexId>::max();

    preprocess_edge_list(edge_list, num_vertices);

    VertexId msf_size = 0;
    std::vector<bool> part_of_msf(num_vertices, false);
    std::vector<Weight> min_weight(num_vertices, INF_WEIGHT);
    std::vector<VertexId> min_parent(num_vertices, INVALID_ID);
    NaiveJPPQ open_vertices(num_vertices);

    auto gen = std::mt19937(std::random_device()());
    auto distr = std::uniform_int_distribution<VertexId>(0, num_vertices - 1);

    while (msf_size < num_vertices) {
      VertexId v;
      if (open_vertices.empty()) {
        // Draw random start vertex for next tree in forest
        do {
          v = distr(gen);
        } while (part_of_msf[v]);
      } else {
        v = open_vertices.delete_and_return_min().second;
      }

      part_of_msf[v] = true;
      ++msf_size;
      for (VertexId i = first_out_edge[v]; i < first_out_edge[v + 1]; ++i) {
        const auto& head = sorted_edges[i].head;
        if (part_of_msf[head]) continue;

        const auto& weight = sorted_edges[i].weight;

        if (min_weight[head] == INF_WEIGHT) {
          open_vertices.insert(weight, head);
          min_weight[head] = weight;
          min_parent[head] = v;
        } else if (weight < min_weight[head]) {
          open_vertices.decrease_key(weight, head);
          min_weight[head] = weight;
          min_parent[head] = v;
        }
      }
    }


    assert(std::for_each(part_of_msf.begin(), part_of_msf.end(),
                         [](const bool is_part) { return is_part; }));

    WEdgeList msf;
    for (VertexId head = 0; head < num_vertices; ++head) {
      const auto& tail = min_parent[head];
      assert(head != tail);
      if (tail == INVALID_ID) continue;
      const auto& weight = min_weight[head];
      msf.emplace_back(tail, head, weight);
      msf.emplace_back(head, tail, weight);
    }

    return msf;
  }

 private:
  // Order the edge list primarily by edge tail and secondarily by edge head and
  // store the index of the first outgoing edge of each vertex.
  void preprocess_edge_list(const algen::WEdgeList& edge_list,
                            const algen::VertexId num_vertices) {
    using namespace algen;
    sorted_edges = edge_list;
    std::sort(sorted_edges.begin(), sorted_edges.end(),
              TailHeadOrder<WEdge>());

    first_out_edge.clear();
    VertexId cur_vertex = 0;
    first_out_edge.push_back(0);
    for (VertexId i = 0; i < sorted_edges.size(); ++i) {
      if (sorted_edges[i].tail != cur_vertex) {
        first_out_edge.insert(first_out_edge.end(),
                              sorted_edges[i].tail - cur_vertex, i);
        cur_vertex = sorted_edges[i].tail;
        assert(first_out_edge.size() == cur_vertex + 1);
      }
    }
    first_out_edge.insert(first_out_edge.end(),
                          num_vertices - first_out_edge.size() + 1,
                          sorted_edges.size());
    assert(first_out_edge.size() == num_vertices + 1 &&
           first_out_edge[num_vertices] == edge_list.size());
  }

  algen::WEdgeList sorted_edges;
  std::vector<algen::VertexId> first_out_edge;
};
