#pragma once

#include <chrono>
#include <cstddef>
#include <string_view>

namespace benchmark {

struct Measurement {
  std::string_view key;
  std::size_t value;
};

class TimeInstrumentation {
public:
  void start() { start_ = std::chrono::steady_clock::now(); }

  void stop() {
    using namespace std::chrono;
    const auto now = steady_clock::now();
    result_[0].value = duration_cast<nanoseconds>(now - start_).count();
  }

  auto begin() const { return std::begin(result_); }
  auto end() const { return std::end(result_); }

private:
  std::chrono::steady_clock::time_point start_;
  Measurement result_[1] = {{"time_ns", 0}};
};

} // namespace benchmark
