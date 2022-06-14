#include "command_line.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>

namespace benchmark {

template <class Options>
void print_help(const char* binary) {
  Options defaults;
  std::cout << "Usage: " << binary
            << " [option]... [output-file]\n"
               "Options:\n"
               "  -h, --help                    Prints this help message and exits.\n"
               "  -i, --iterations              Number of iterations ran per contender and input graph (default: 1).\n"
               "  -o, --output                  Path to an output file in CSV format (default: mst_construction.csv / mst_verification.csv).\n"
               "\n";
}

template <class Options>
Options parse(int argc, char** argv) try {
  using command_line::opt;
  Options options = command_line::parse_options<Options>(
      argc, argv,
      opt{"-h", "--help", [&] {
            print_help<Options>(argv[0]);
            std::exit(0);
          }},
          opt{"-i", "--iterations", &Options::iterations},
          opt{"", "--no-verification", &Options::no_verification},
          opt{"-o", "--output", &Options::output_file});

  if (options.output_file.extension() != ".csv") {
      options.output_file += ".csv";
  }
  return options;
} catch (std::exception& ex) {
  std::cerr << "Error while parsing command line: " << ex.what() << '\n';
  std::exit(-1);
} catch (...) {
  std::cerr << "Unknown exception occurred while parsing command line.\n";
  std::exit(-1);
}

struct CsvOutput {
      std::ofstream file;
      std::stringstream line;

      explicit CsvOutput(const std::filesystem::path& path) {
          file.exceptions(std::ios_base::badbit | std::ios_base::failbit);
          file.open(path);
          line.precision(3);
      }

      template<class I>
      void print_mst_construction_header(const I& instr) {
          const std::vector<std::string> config_column_names = {"contender", "input_graph", "iteration", "num_vertices", "num_edges"};
          print_header(config_column_names, instr);
      }

    template<class I>
    void print_mst_verification_header(const I& instr) {
        const std::vector<std::string> config_column_names = {"contender", "input_graph", "iteration", "num_vertices", "num_edges", "num_changed_edges"};
        print_header(config_column_names, instr);
    }

      template <class I>
      void print_header(const std::vector<std::string>& config_column_names, const I& instr) {
          line.str("");
          for (const auto& column_name : config_column_names)
              line << column_name << ", ";
          for (auto m : instr)
              line << "factory_" << m.key << ", ";
          for (auto m : instr)
              line << "execution_" << m.key << ", ";
          line << "result_correct\n";
          file << line.str();
      }



      template <class... Ts>
      void add(Ts... fields) {
          for_each(std::tie(fields...), [&](auto val) {
              line << val << ", ";
          });
      }

      void finish(std::string_view field) {
          line << field << '\n';
          file << line.str();
      }
    };

} // namespace benchmark
