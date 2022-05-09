#pragma once

#include <algorithm>
#include <cassert>
#include "datastructures/naive_union_find.hpp"
#include "includes/definitions.hpp"

// This implementation of the Kruskal algorithm is deliberately _very_ bad. It
// is only here to be used as a naive contender, showing how to use the
// contender interface.
struct NaiveKruskal {
  algen::WEdgeList operator()(const algen::WEdgeList& edge_list, const algen::VertexId num_vertices) {
    using namespace algen;
    WEdgeList el_copy = edge_list;
    WEdgeList mst;

    UnionFind<VertexId> subtrees(num_vertices);

    std::sort(el_copy.begin(), el_copy.end(), [&](const auto e1, const auto e2) {
      return e1.weight < e2.weight;
    });

    for (auto& e : el_copy) {
      const auto tail_subtree = subtrees.find(e.tail);
      const auto head_subtree = subtrees.find(e.head);
      if (tail_subtree != head_subtree) {
        mst.push_back(e);
          subtrees.do_union(tail_subtree, head_subtree);
      }
    }

    // Add reverse edges
    const auto mst_one_way_size = mst.size();
    mst.reserve(2 * mst_one_way_size);
    for (VertexId i = 0; i < mst_one_way_size; ++i) {
      mst.emplace_back(mst[i].head, mst[i].tail, mst[i].weight);
    }

    return mst;
  }
};
