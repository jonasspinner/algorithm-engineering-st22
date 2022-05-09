#pragma once

#include <charconv>
#include <cstddef>
#include <exception>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace command_line {
namespace detail {

template <class T>
using remove_cv_ref_t = std::remove_cv_t<std::remove_reference_t<T>>;

template <class T>
struct function_traits : function_traits<decltype(&std::remove_reference_t<T>::operator())> {};

template <class R, class... As>
struct function_traits<R(As...)> {
  using args_t = std::tuple<As...>;
};

template <class R, class... As>
struct function_traits<R (*)(As...)> : function_traits<R(As...)> {};

template <class C, class R, class... As>
struct function_traits<R (C::*)(As...)> : function_traits<R(As...)> {};

template <class C, class R, class... As>
struct function_traits<R (C::*)(As...) const> : function_traits<R(As...)> {};

template <class W, class T>
struct count_in_tuple : std::integral_constant<std::size_t, 0> {};
template <class W, class H, class... Ts>
struct count_in_tuple<W, std::tuple<H, Ts...>>
    : std::integral_constant<std::size_t, (std::is_same_v<W, remove_cv_ref_t<H>> ? 1 : 0) +
                                              count_in_tuple<W, std::tuple<Ts...>>{}> {};

template <class T, class O>
decltype(auto) parse_arg([[maybe_unused]] O* opts, char** args, std::size_t idx, std::size_t& offset) {
  if constexpr (std::is_same_v<T, O>) {
    ++offset;
    return static_cast<O&>(*opts);
  } else {
    const std::string_view arg = args[idx - offset];
    if constexpr (std::is_arithmetic_v<T>) {
      // Integer and floating point values
      std::conditional_t<std::is_same_v<T, bool>, int, T> value{};
      const auto r = std::from_chars(arg.data(), arg.data() + arg.size(), value);
      if (r.ec != std::errc() || r.ptr != arg.data() + arg.size())
        throw std::system_error{std::make_error_code(r.ec), std::string(arg)};
      return value;
    } else if constexpr (std::is_same_v<T, const char*> || std::is_same_v<T, char*>) {
      // Strings
      return arg;
    } else if constexpr (std::is_same_v<T, std::vector<std::string_view>> ||
                         std::is_same_v<T, std::vector<std::string>>) {
      // Vector of comma-separated strings
      T vec;
      for (std::size_t i = 0; i < arg.size();) {
        auto j = arg.find(',', i);
        if (j == std::string_view::npos)
          j = arg.size();
        vec.emplace_back(arg.data() + i, j - i);
        i = j + 1;
      }
      return vec;
    } else if constexpr (std::is_constructible_v<T, std::string_view>) {
      // Strings, and things constructible from strings
      return T{arg};
    } else {
      static_assert(!std::is_same_v<T, T>, "Unsupported argument type");
    }
  }
}

template <class O, class T>
struct arg_parser;
template <class O, class... As>
struct arg_parser<O, std::tuple<As...>> {
  template <std::size_t... Is>
  static auto parse_args(O* opts, char** args, std::index_sequence<Is...>) {
    std::size_t offset = 0;
    return std::tuple<decltype(parse_arg<remove_cv_ref_t<As>>(opts, args, Is, offset))...>{
        parse_arg<remove_cv_ref_t<As>>(opts, args, Is, offset)...};
  }
};

template <class T>
struct ptr_handler {
  using type = T;
};

template <class T>
struct ptr_handler<T*> {
  struct type {
    type() = delete;
    explicit type(T* p) : ptr(p) {}
    T* ptr;
    void operator()(T&& v) const { *ptr = std::move(v); }
  };
};

template <class C, class M>
struct ptr_handler<M C::*> {
  struct type {
    type() = delete;
    explicit type(M C::*p) : ptr(p) {}
    M C::*ptr;
    void operator()(C& c, M&& v) const { c.*ptr = std::move(v); }
  };
};

template <class C, class... As>
struct ptr_handler<std::vector<As...> C::*> {
  struct type {
    type() = delete;
    explicit type(std::vector<As...> C::*p) : ptr(p) {}
    std::vector<As...> C::*ptr;
    void operator()(C& c, std::vector<As...>&& v) const {
      (c.*ptr).insert((c.*ptr).end(), std::move(v).begin(), std::move(v).end());
    }
  };
};

template <class C, class R, class... As>
struct ptr_handler<R (C::*)(As...)> {
  struct type {
    type() = delete;
    explicit type(R (C::*p)(As...)) : ptr(p) {}
    R (C::*ptr)(As...);
    void operator()(C& c, As&&... as) const { (c.*ptr)(std::move(as)...); }
  };
};

template <class T>
using ptr_handler_t = typename ptr_handler<T>::type;

template <class T, class F, std::size_t I = 0>
constexpr bool arg_type_is = std::is_same_v<remove_cv_ref_t<std::tuple_element_t<I, typename F::args_t>>, T>;

} // namespace detail

/**
 * Defines a command line option, using short form, long form and option
 * handler. Handler may be a function, a pointer, or a member pointer.
 */
template <class F>
struct opt {
  using handler_t = detail::ptr_handler_t<F>;
  using args_t = typename detail::function_traits<handler_t>::args_t;
  static constexpr std::size_t num_args = std::tuple_size_v<args_t>;

  opt(std::string_view so, std::string_view lo, F f, bool repeat = false)
      : short_opt(so), long_opt(lo), handler(std::move(f)), can_repeat(repeat) {}

  std::string_view short_opt;
  std::string_view long_opt;
  handler_t handler;
  bool can_repeat;
};

/**
 * Parses command line options.
 */
template <class O = void, class... Fs>
O parse_options(int argc, char** argv, opt<Fs>&&... fs) {
  std::conditional_t<std::is_void_v<O>, bool, O> opts;

  int current = 1;
  [[maybe_unused]] const auto try_parse = [&](auto&& f) -> bool {
    using F = detail::remove_cv_ref_t<decltype(f)>;
    const std::string_view arg = argv[current];

    if (arg[0] != '-') {
      // Positional arguments
      if (!f.short_opt.empty() || !f.long_opt.empty())
        return false;
      --current;
    } else if (arg != f.short_opt && arg != f.long_opt) {
      return false;
    }

    if (!f.can_repeat) {
      f.short_opt = " ";
      f.long_opt = " ";
    }

    // Handle boolean flags
    if constexpr (F::num_args == 0) {
      // Call handler without arguments
      std::invoke(f.handler);
      return true;
    } else if constexpr (F::num_args == 1) {
      if constexpr (detail::arg_type_is<bool, F>) {
        // Call handler with single boolean argument. Usually used for pointers.
        std::invoke(f.handler, true);
        return true;
      }
    } else if constexpr (F::num_args == 2 && !std::is_void_v<O>) {
      if constexpr (detail::arg_type_is<O, F, 0> && detail::arg_type_is<bool, F, 1>) {
        // Call handler with Options + boolean arguments. Usually used for
        // member pointers.
        std::invoke(f.handler, opts, true);
        return true;
      }
    }

    static constexpr int req_args = F::num_args - detail::count_in_tuple<O, typename F::args_t>();

    if (argc - current - 1 < req_args)
      throw std::runtime_error{"Missing argument to " + std::string(arg)};

    try {
      std::apply(f.handler, detail::arg_parser<O, typename F::args_t>::parse_args(
                                &opts, argv + current + 1, std::make_index_sequence<F::num_args>{}));
    } catch (std::runtime_error& ex) {
      std::throw_with_nested(
          std::runtime_error{"Failed to parse argument to " + std::string(arg) + ": " + std::string(ex.what())});
    }
    current += req_args;
    return true;
  };

  for (; current < argc; ++current)
    if (!(try_parse(fs) || ...))
      throw std::runtime_error{"Unknown option: " + std::string(argv[current])};

  if constexpr (!std::is_void_v<O>)
    return opts;
}

} // namespace command_line
