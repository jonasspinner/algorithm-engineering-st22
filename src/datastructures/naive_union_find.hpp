#pragma once

#include <vector>
#include <numeric>

// This implementation of the Union-Find structure is deliberately _very_ bad. It is only here to be used with
// NaiveKruskal.
template<typename T>
class NaiveUnionFind {

    static_assert(std::is_integral_v<T> && std::is_unsigned_v<T>);

public:

    explicit NaiveUnionFind(const T& num_elements) {
        initialize(num_elements);
    }

    T find(T el) const {
        while (parent[el] != el) {
            el = parent[el];
        }
        return el;
    }

    void do_union(const T& el1, const T& el2) {
        assert(find(el1) != find(el2));
        parent[find(el1)] = find(el2);
    }

private:

    void initialize(const T& num_elements) {
        parent.resize(num_elements);
        std::iota(parent.begin(), parent.end(), 0);
    }

    std::vector<T> parent;
};

template<typename T>
using UnionFind = NaiveUnionFind<T>;
