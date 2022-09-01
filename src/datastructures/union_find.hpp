#pragma once

#include <vector>
#include <numeric>
#include <cstdint>
#include <cassert>

namespace js {

    template<typename T>
    class UnionFindOld {
        static_assert(std::is_integral_v<T>);
    public:
        explicit UnionFindOld(const T &num_elements) : m_parents(num_elements), m_ranks(num_elements, 0) {
            std::iota(m_parents.begin(), m_parents.end(), T());
        }

        T find(T element) {
            assert(element < m_parents.size());

            if (m_parents[element] != element)
                m_parents[element] = find(m_parents[element]);
            return m_parents[element];
        }

        bool do_union(T a, T b) {
            assert(a < m_parents.size());
            assert(b < m_parents.size());

            a = find(a);
            b = find(b);

            if (a == b) { return false; }
            if (m_ranks[a] > m_ranks[b]) std::swap(a, b);
            m_parents[a] = b;

            if (m_ranks[a] == m_ranks[b]) {
                m_ranks[b]++;
            }
            return true;
        }

        bool is_singleton(T element) const {
            return m_parents[element] == element && m_ranks[element] == 0;
        }

        [[nodiscard]] size_t num_elements() const {
            assert(m_parents.size() == m_ranks.size());
            return m_parents.size();
        }

    private:
        std::vector<T> m_parents;
        std::vector<uint8_t> m_ranks;
    };

    template<typename T>
    class UnionFind {
        static_assert(std::is_integral_v<T>);
    public:
        explicit UnionFind(const T &num_elements) : m_parents_ranks(num_elements, num_elements) {}

        T find(T element) {
            assert(element < num_elements());
            if (m_parents_ranks[element] >= num_elements()) {
                return element;
            }
            m_parents_ranks[element] = find(m_parents_ranks[element]);
            return m_parents_ranks[element];
        }

        bool do_union(T a, T b) {
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

        bool is_non_singleton_representative(T element) const {
            return m_parents_ranks[element]  > num_elements();
        }

        /**
         * Call f for all elements that have at least one child.
         * @tparam F
         * @param f
         */
        template<class F>
        void for_cluster_representatives(F f) const {
            if (m_cluster_representatives.empty()) {
                // This function is called for the first time.
                m_cluster_representatives.reserve(num_elements() / 2);
                for (T element = 0; element < num_elements(); ++element) {
                    if (m_parents_ranks[element] > num_elements()) {
                        m_cluster_representatives.push_back(element);
                        f(element);
                    }
                }
                return;
            }

            size_t last = 0;
            for (size_t i = 0; i < m_cluster_representatives.size(); ++i) {
                auto element = m_cluster_representatives[i];
                if (m_parents_ranks[element] > num_elements()) {
                    m_cluster_representatives[last] = element;
                    last++;
                    f(element);
                }
            }
            m_cluster_representatives.erase(m_cluster_representatives.begin() + last, m_cluster_representatives.end());
        }

        [[nodiscard]] size_t num_elements() const {
            return m_parents_ranks.size();
        }

    private:
        std::vector<T> m_parents_ranks;

        mutable std::vector<T> m_cluster_representatives;
    };
}