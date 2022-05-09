#pragma once

namespace benchmark {

    template <class F, class T, T... vals>
    constexpr void for_each(std::integer_sequence<T, vals...>, F&& func) {
        (..., func(std::integral_constant<T, vals>{}));
    }

    template <class F, class... Ts, std::size_t... Is>
    constexpr void for_each(const std::tuple<Ts...>& values, std::index_sequence<Is...>, F&& func) {
        (..., func(std::get<Is>(values)));
    }

    template <class F, class... Ts>
    constexpr void for_each(const std::tuple<Ts...>& values, F&& func) {
        for_each(values, std::index_sequence_for<Ts...>{}, std::forward<F>(func));
    }

} // end namespace benchmark
