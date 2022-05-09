#include <filesystem>
#include "mst_verification_parameters.hpp"
#include "includes/definitions.hpp"
#include "includes/utils.hpp"
#include "benchmark_helpers.hpp"
#include "graph_generation.hpp"
#include "instrumentation.hpp"
#include "for_each.hpp"
#include "naive_jarnik_prim.hpp"
#include "verification_instance_generation.hpp"
#include "includes/binary_includes.hpp"

namespace {
    struct Options {
        std::filesystem::path output_file = "mst_verification.csv";
        std::size_t iterations = mst_verification::params::iterations;
        bool no_verification = false;
    };
}  // end anonymous namespace

template<class I, class C, class V>
void execute_benchmark(I& instrumentation, const algen::WEdgeList & graph_edge_list, const algen::WEdgeList& st_edge_list,
                       const int num_vertices, C& contender, std::size_t iterations, V& verify_response,
                       benchmark::CsvOutput& out, std::string_view input_graph_name) {
    for (std::size_t i = 0; i < iterations; ++i) {

        out.line.str("");
        out.add(contender.name, input_graph_name, i);

        instrumentation.start();
        auto contender_instance = contender.factory();
        instrumentation.stop();
        for (auto m : instrumentation)
            out.add(m.value);

        instrumentation.start();
        std::optional<algen::WEdge> response = contender_instance(graph_edge_list, st_edge_list, num_vertices);
        instrumentation.stop();
        for (auto m : instrumentation)
            out.add(m.value);

        auto [correct, msg] = verify_response(response);
        if (correct) {
            out.finish("yes");
        } else {
            std::stringstream incorrect_msg;
            incorrect_msg << "no: " << msg;
            out.finish(incorrect_msg.str());
        }

    }
}

static bool verify_input_list(const algen::WEdgeList& edge_list, const algen::VertexId num_vertices) {
    const auto [correct, msg] = edge_list_format_check(edge_list, num_vertices);
    if (!correct) {
        std::cerr << "Input edge list does not fulfill requirements: " << msg
                  << std::endl;
    }
    return correct;
}

static algen::WEdgeList construct_st_instance(const algen::WEdgeList& graph_edges,
                                              const algen::VertexId num_vertices,
                                              const std::size_t num_changed_edges) {
    // Construct an ST instance as input for the verification contenders by computing an
    // MST and corrupting it in num_changed_edges random places.
    using namespace algen;
    auto call_fast_kruskal = [](const WEdgeList& edges, const VertexId n) { return fast_kruskal(edges, n);};
    benchmark::CorruptedMSTGenerator<decltype(call_fast_kruskal)> instance_gen;
    static constexpr std::size_t seed = 123; // seed for RNG
    instance_gen.preprocess(graph_edges, num_vertices, call_fast_kruskal, seed, true);
    return instance_gen.generate_corrupted_mst(num_changed_edges, true);
}

std::pair<bool, std::string> verify_result(const algen::WEdgeList& graph_edges,
                                           const algen::WEdgeList & tree_edges,
                                           const algen::VertexId num_vertices,
                                           const std::optional<algen::WEdge> result) {
    using namespace algen;
    if (!result) {
        // Algorithm says MST is correct so check that it actually is
        const auto is_mst_incorrect = verify_spanning_tree(graph_edges, tree_edges, num_vertices);
        if (is_mst_incorrect) {
            return {false, "Algorithm says MST is correct but it is not."};
        }
        return {true, "Correctly identified correct MST."};
    }

    // Algorithm says edge e is light wrt to the spanning tree so check that it actually is using a DFS on the tree
    const WEdge e = result.value();

    auto sorted_st_edges = tree_edges;
    std::vector<EdgeIdx> first_st_out_edge;
    make_inefficient_adjacency_structure(sorted_st_edges, first_st_out_edge, num_vertices);

    if (!NaiveDFSBasedVerifier::is_st_light(e, sorted_st_edges, first_st_out_edge)) {
        std::stringstream msg;
        msg << "Algorithm says edge (" << e.tail << " - " << e.head << ") is ST light but it actually is not.";
        return {false, msg.str()};
    }
    return {true, "Correctly identified incorrect MST."};
}

int main(int argc, char** argv) try {
    using namespace algen;
    const auto options = benchmark::parse<Options>(argc, argv);

    // Generate random input graph with size specified in parameters
    benchmark::GNM_Generator generator;
    constexpr std::size_t log_n = mst_verification::params::graph_generator_params.log_n;
    constexpr std::size_t log_m = mst_verification::params::graph_generator_params.log_m;
    constexpr Weight max_edge_weight = mst_verification::params::graph_generator_params.max_weight;
    generator.configure(log_n, log_m, max_edge_weight);
    auto gen_edges = generator.generate();
    std::stringstream graph_name;
    graph_name << "generated_graph_" << log_n << "_" << log_m << "_"
               << max_edge_weight;

    static constexpr bool print_generated_edges = false;
    if (print_generated_edges) {
        print_container(gen_edges);
    }

    constexpr VertexId num_vertices = 1ull << log_n;
    assert(verify_input_list(gen_edges, num_vertices));

    // Construct a spanning tree instance for the MST verification contenders on the generated input graph
    const auto corrupted_mst_edge_list = construct_st_instance(gen_edges, num_vertices, mst_verification::params::num_changed_edges);

    // Prepare time measurement, output and result verification
    benchmark::TimeInstrumentation instrumentation;
    benchmark::CsvOutput output(options.output_file);
    output.print_header(instrumentation);
    const auto verify = [&](const std::optional<WEdge> response) -> std::pair<bool, std::string> {
        return verify_result(gen_edges, corrupted_mst_edge_list, num_vertices, response);
    };

    // Execute all contenders on the generated spanning tree
    benchmark::for_each(
            mst_verification::params::contenders, [&](auto& contender) {
                execute_benchmark(instrumentation, gen_edges, corrupted_mst_edge_list, num_vertices,
                                             contender, options.iterations, verify,
                                             output, graph_name.str());
            });

    return 0;
} catch (std::exception& ex) {
    std::cerr << "Unhandled exception: " << ex.what() << '\n';
    return -2;
}