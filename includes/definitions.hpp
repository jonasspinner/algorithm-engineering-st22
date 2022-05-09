#pragma once

#include <cstdint>
#include <limits>
#include <ostream>
#include <vector>

namespace algen {
using VertexId = uint64_t;
using Weight = int32_t;
using EdgeIdx = uint64_t;

static constexpr VertexId VERTEXID_UNDEFINED =
    std::numeric_limits<VertexId>::max();
static constexpr VertexId VERTEXID_MAX =
    std::numeric_limits<VertexId>::max() - 1;
static constexpr Weight WEIGHT_UNDEFINED = std::numeric_limits<Weight>::max();
static constexpr Weight WEIGHT_MAX = std::numeric_limits<Weight>::max() - 1;

using VertexArray = std::vector<VertexId>;
using WeightVertex = std::pair<Weight, VertexId>;
using WeightArray = std::vector<Weight>;

struct UnweightedEdge {
  VertexId tail;
  VertexId head;
};

struct Edge {
  Edge() = default;
  Edge(VertexId tail, VertexId head) : tail{tail}, head{head} {}
  VertexId tail;
  VertexId head;
  friend std::ostream& operator<<(std::ostream& out, const Edge edge) {
    return out << "(" << edge.tail << ", " << edge.head << ")";
  }
};

struct WEdge {
  WEdge() = default;
  WEdge(VertexId tail, VertexId head, Weight weight)
      : tail{tail}, head{head}, weight{weight} {}
  VertexId tail;
  VertexId head;
  Weight weight;
  friend std::ostream& operator<<(std::ostream& out, const WEdge edge) {
    return out << "(" << edge.tail << ", " << edge.head << ", " << edge.weight
               << ")";
  }
};

using WEdgeList = std::vector<WEdge>;
}  // namespace algen
