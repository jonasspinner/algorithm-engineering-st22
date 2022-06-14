
#include <exception>
#include <iostream>
#include <sstream>
#include <vector>

#include "benchmark_helpers.hpp"
#include "for_each.hpp"
#include "graph_generation.hpp"
#include "includes/binary_includes.hpp"
#include "includes/utils.hpp"
#include "instrumentation.hpp"
#include "mst_construction_parameters.hpp"
#include "verification_instance_generation.hpp"

namespace {
struct Options {
  std::filesystem::path output_file = "mst_construction.csv";
  std::size_t iterations = mst_construction::params::iterations;
  bool no_verification = false;
};
}  // end anonymous namespace

template <class I, class C, class V>
void execute_benchmark(I& instrumentation, const algen::WEdgeList& edge_list,
                       const int num_vertices, C& contender,
                       std::size_t iterations, V& verify,
                       benchmark::CsvOutput& out,
                       std::string_view input_graph_name) {
  for (std::size_t i = 0; i < iterations; ++i) {
    out.line.str("");
    out.add(contender.name, input_graph_name, i, num_vertices,
            edge_list.size() / 2);

    instrumentation.start();
    auto contender_instance = contender.factory();
    instrumentation.stop();
    for (auto m : instrumentation) out.add(m.value);

    instrumentation.start();
    auto mst_edges = contender_instance(edge_list, num_vertices);
    instrumentation.stop();
    for (auto m : instrumentation) out.add(m.value);

    auto [correct, msg] = verify(mst_edges);
    if (correct) {
      out.finish("yes");
    } else {
      std::stringstream incorrect_msg;
      incorrect_msg << "no: " << msg;
      out.finish(incorrect_msg.str());
    }
  }
}
static bool verify_input_list(const algen::WEdgeList& edge_list,
                              const algen::VertexId num_vertices) {
  const auto [correct, msg] = edge_list_format_check(edge_list, num_vertices);
  if (!correct) {
    std::cerr << "Input edge list does not fulfill requirements: " << msg
              << std::endl;
  }
  return correct;
}

std::pair<bool, std::string> verify_result(const algen::WEdgeList& graph_edges,
                                           const algen::WEdgeList& result,
                                           const algen::VertexId num_vertices) {
  const auto [format_correct, format_msg] =
      edge_list_format_check(result, num_vertices);
  if (!format_correct) {
    return {false, format_msg};
  }

  const auto [is_spanning_tree, st_msg] =
      is_spanning_forest(graph_edges, result, num_vertices);
  if (!is_spanning_tree) {
    return {false, st_msg};
  }

  const auto counter_example =
      verify_spanning_tree(graph_edges, result, num_vertices);
  if (counter_example) {
    std::stringstream msg;
    msg << "Edge (" << counter_example->tail << " - " << counter_example->head
        << ") is light wrt the ST returned.";
    return {false, msg.str()};
  }
  return {true, "Returned ST is correct MST."};
}

int main(int argc, char** argv) try {
  using namespace algen;
  const auto options = benchmark::parse<Options>(argc, argv);

  // Generate random input graph with size specified in parameters
  benchmark::GNM_Generator generator;
  mst_construction::params::ExperimentSuite experiments;
  benchmark::CsvOutput output(options.output_file);
  benchmark::TimeInstrumentation instrumentation;
  output.print_mst_construction_header(instrumentation);
  while (experiments.has_next()) {
    const auto experiment = experiments.get_next();
    std::cout << experiment << std::endl;
    const std::size_t log_n = experiment.log_n;
    const std::size_t log_m =
        experiment.log_n + std::log2(experiment.edge_factor);
    const algen::Weight max_edge_weight = experiment.max_weight;

    // generate graph
    generator.configure(log_n, log_m, max_edge_weight);
    auto gen_edges = generator.generate();
    static constexpr bool print_generated_edges = false;
    if (print_generated_edges) {
      print_container(gen_edges);
    }

    VertexId num_vertices = 1ull << log_n;
    std::cout << "\tstart verification input" << std::endl;
    assert(verify_input_list(gen_edges, num_vertices));
    std::cout << "\tstop verification input" << std::endl;

    // Prepare result verification
    const auto verify = [&](WEdgeList& result) -> std::pair<bool, std::string> {
      std::cout << "\tstart verification result" << std::endl;
      auto res = verify_result(gen_edges, result, num_vertices);
      std::cout << "\tstop verification result" << std::endl;
      return res;
    };

    // Execute all contenders on the generated input graph
    benchmark::for_each(
        mst_construction::params::contenders, [&](auto& contender) {
          execute_benchmark(instrumentation, gen_edges, num_vertices, contender,
                            options.iterations, verify, output,
                            generator.name());
        });
  }
  return 0;
} catch (std::exception& ex) {
  std::cerr << "Unhandled exception: " << ex.what() << '\n';
  return -2;
}
