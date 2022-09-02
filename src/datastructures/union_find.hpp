#pragma once

#include <vector>
#include <numeric>
#include <cstdint>
#include <cassert>

namespace js {
    template<typename T>
    class UnionFind {
        static_assert(std::is_integral_v<T>);
    public:
        UnionFind() = default;

        explicit UnionFind(T num_elements) : m_parents_ranks(num_elements, num_elements), m_num_representatives(num_elements) {}

        void reset(T num_elements) {
            m_parents_ranks.clear();
            m_parents_ranks.resize(num_elements, num_elements);
            m_num_representatives = num_elements;
        }

        constexpr T find(T element) noexcept {
            assert(element < num_elements());
            if (m_parents_ranks[element] >= num_elements()) {
                return element;
            }
            return m_parents_ranks[element] = find(m_parents_ranks[element]);
        }

        constexpr bool do_union(T a, T b) noexcept {
            assert(a < num_elements());
            assert(b < num_elements());

            a = find(a);
            b = find(b);

            assert(m_parents_ranks[a] >= num_elements());
            assert(m_parents_ranks[b] >= num_elements());

            if (a == b) { return false; }

            auto rank_a = m_parents_ranks[a];
            auto rank_b = m_parents_ranks[b];
            if (rank_a > rank_b) std::swap(a, b);
            m_parents_ranks[a] = b;

            if (rank_a == rank_b) {
                m_parents_ranks[b]++;
            }
            return true;
        }

        bool is_singleton(T element) const {
            return m_parents_ranks[element] == num_elements();
        }

        bool is_representative(T element) const {
            return m_parents_ranks[element] >= num_elements();
        }

        [[nodiscard]] constexpr bool is_non_singleton_representative(T element) const noexcept {
            return m_parents_ranks[element] > num_elements();
        }

        [[nodiscard]] size_t num_elements() const {
            return m_parents_ranks.size();
        }

    private:
        std::vector<T> m_parents_ranks;
        T m_num_representatives{0};
    };
}