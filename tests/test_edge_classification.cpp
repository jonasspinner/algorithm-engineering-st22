// just some quick and dirty test to "ensure" correct functionality

#include "benchmark/graph_generation.hpp"

algen::WEdgeList filter_out_edges(const std::vector<bool>& is_edge_light,
                                  const algen::WEdgeList& edges) {
  algen::WEdgeList light_edges;
  for (std::size_t i = 0; i < is_edge_light.size(); ++i) {
    if (is_edge_light[i]) light_edges.push_back(edges[i]);
  }
  return light_edges;
}

bool test_edge_classification(std::size_t log_n, std::size_t log_m) {
  benchmark::GNM_Generator generator;
  const std::size_t max_edge_weight = 255;
  generator.configure(log_n, log_m, max_edge_weight);
  auto gen_edges = generator.generate();

  const auto mst_org = fast_kruskal(gen_edges, 1ull << log_n);
  auto is_edge_light =
      algen::getEdgeClassifier().execute(gen_edges, mst_org, 1ull << log_n);
  auto num_heavy_edges =
      std::count(is_edge_light.begin(), is_edge_light.end(), false);
  std::cout << "log_n: " << log_n << " log_m: " << log_m
            << " #mst edges: " << mst_org.size()
            << " #heavy edges: " << num_heavy_edges << std::endl;
  const auto light_edges = filter_out_edges(is_edge_light, gen_edges);
  if (gen_edges.size() - num_heavy_edges != light_edges.size()) {
    std::cout << "\t\twrong filter result" << std::endl;
    return false;
  }

  const auto mst_based_on_light = fast_kruskal(light_edges, 1ull << log_n);
  if (algen::sum_weights(mst_based_on_light) != algen::sum_weights(mst_org)) {
    std::cout << "\t\twrong mst weight" << std::endl;
    return false;
  }
  return true;
}

int main() {
  test_edge_classification(15, 10);
  test_edge_classification(15, 15);
  test_edge_classification(15, 20);
  test_edge_classification(17, 10);
  test_edge_classification(17, 15);
  test_edge_classification(17, 20);
}
