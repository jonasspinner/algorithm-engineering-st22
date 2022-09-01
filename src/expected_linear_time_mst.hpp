#pragma once

#include <utility>
#include <random>
#include <algorithm>

#include "includes/definitions.hpp"
#include "includes/utils.hpp"

#include "datastructures/union_find.hpp"
#include "datastructures/robin_hood.h"

#include "datastructures/fast_reset_bitvector.hpp"

namespace js {
    template<typename T>
    void hash_combine(std::size_t &seed, T const &v) {
        // The Boost Software License, Version 1.0 applies to this function.
        // See https://www.boost.org/LICENSE_1_0.txt
        // https://www.boost.org/doc/libs/1_75_0/doc/html/hash/reference.html#boost.hash_combine
        seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    struct PairHash {
        template<typename A, typename B>
        std::size_t operator()(const std::pair<A, B> &pair) const {
            std::size_t seed = 0;
            hash_combine(seed, pair.first);
            hash_combine(seed, pair.second);
            return seed;
        }
    };


    class ExpectedLinearTimeMST {
        using i64 = int64_t;
        using u64 = uint64_t;
    public:
        algen::WEdgeList operator()(const algen::WEdgeList &edge_list,
                                    const algen::VertexId num_vertices) {
            using namespace algen;

            std::random_device rd;
            std::uniform_int_distribution<size_t> seed_dist;
            m_gen.seed(seed_dist(rd));
            //m_gen.seed(0);

            auto G = ContractedEdgeListGraph::from_algen_edge_list(edge_list, num_vertices);
            algen::WEdgeList mst;
            mst.reserve(num_vertices - 1);

            m_cheapest_edges_scratchpad.resize(G.num_nodes, INVALID_CONTRACTED_EDGE);
            m_active_nodes_scratchpad.reserve(G.num_nodes / 2);

            algorithm(std::move(G), mst);

            return mst;
        }

    private:
        struct ContractedEdge {
            u64 a;
            u64 b;
            u64 u;
            u64 v;
            algen::Weight weight;

            friend std::ostream &operator<<(std::ostream &os, const ContractedEdge &edge) {
                return os << "(" << edge.a << "," << edge.b << ") " << algen::WEdge{edge.u, edge.v, edge.weight};
            }
        };

        static constexpr ContractedEdge INVALID_CONTRACTED_EDGE = {
                algen::VERTEXID_UNDEFINED,
                algen::VERTEXID_UNDEFINED,
                algen::VERTEXID_UNDEFINED,
                algen::VERTEXID_UNDEFINED,
                algen::WEIGHT_UNDEFINED
        };

        struct ContractedEdgeListGraph {
            std::vector<ContractedEdge> edges;
            u64 num_nodes{};

            static ContractedEdgeListGraph from_algen_edge_list(const std::vector<algen::WEdge> &edges, u64 num_nodes) {
                ContractedEdgeListGraph G;
                G.edges.reserve(edges.size());
                G.num_nodes = num_nodes;
                for (const auto &[u, v, w]: edges) {
                    assert(u < num_nodes);
                    assert(v < num_nodes);
                    G.edges.push_back({u, v, u, v, w});
                }
                return G;
            }

            [[nodiscard]] bool check() const {
                std::unordered_map<std::pair<u64, u64>, std::pair<algen::WEdge, u64>, PairHash> inter_cluster_edges;
                for (auto [a, b, u, v, w]: edges) {
                    if (a == b) {
                        assert(false);
                        return false;
                    }
                    if (a > b) {
                        std::swap(a, b);
                        std::swap(u, v);
                    }
                    auto [it, inserted] = inter_cluster_edges.insert({{a,         b},
                                                                      {{u, v, w}, 1}});
                    if (!inserted) {
                        auto [a_, b_] = it->first;
                        auto [edge_, count] = it->second;
                        auto [u_, v_, w_] = edge_;
                        if (a_ > b_) {
                            std::swap(a_, b_);
                            std::swap(u_, v_);
                        }
                        if (count != 1) {
                            assert(false);
                            return false;
                        }
                        if (!(a == a_ && b == b_ && u == u_ && v == v_ && w == w_)) {
                            assert(false);
                            return false;
                        }
                        it->second.second++;
                    }
                }
                return true;
            }

            void compress(UnionFind<algen::VertexId> &union_find) {
                std::vector<u64> cc_ids(num_nodes, 0);
                u64 current_id = 0;
                for (u64 u = 0; u < num_nodes; ++u) {
                    if (union_find.find(u) == u && !union_find.is_singleton(u)) {
                        cc_ids[u] = current_id;
                        current_id++;
                    }
                }

                auto num_clusters = current_id;
                std::cout << std::setprecision(4) << std::setw(5)
                          << 100. * (1. - (double) num_clusters / (double) num_nodes) << "% ("
                          << std::setw(8) << num_nodes - num_clusters << "/"
                          << std::setw(8) << num_nodes << ") of nodes removed by compression\n";

                num_nodes = num_clusters;


                robin_hood::unordered_map<std::pair<u64, u64>, algen::WEdge, PairHash> inter_cluster_edges;
                inter_cluster_edges.reserve(num_nodes);

                for (auto &edge: edges) {
                    edge.a = union_find.find(edge.a);
                    edge.b = union_find.find(edge.b);
                    if (edge.a != edge.b) {
                        edge.a = cc_ids[edge.a];
                        edge.b = cc_ids[edge.b];

                        if (edge.a > edge.b) {
                            std::swap(edge.a, edge.b);
                        }
                        assert(edge.a != edge.b);
                        assert(edge.a < num_nodes);
                        assert(edge.b < num_nodes);
                        algen::WEdge uvw = {edge.u, edge.v, edge.weight};
                        auto [it, inserted] = inter_cluster_edges.emplace(std::pair{edge.a, edge.b}, uvw);
                        if (!inserted) {
                            if (edge.weight < it->second.weight) {
                                it->second = uvw;
                            }
                        }
                    }
                }

                std::cout << std::setprecision(4) << std::setw(5)
                          << 100. * (1. - (2. * (double) inter_cluster_edges.size()) / (double) edges.size()) << "% ("
                          << std::setw(8) << edges.size() / 2 - inter_cluster_edges.size() << "/"
                          << std::setw(8) << edges.size() / 2 << ") of edges removed as multi edges\n";

                edges.clear();
                edges.reserve(2 * inter_cluster_edges.size());
                for (const auto &[ab, edge]: inter_cluster_edges) {
                    auto [a, b] = ab;
                    assert(a < num_nodes);
                    assert(b < num_nodes);
                    assert(a != b);
                    auto [u, v, w] = edge;
                    edges.push_back({a, b, u, v, w});
                    edges.push_back({b, a, v, u, w});
                }
            }

            [[nodiscard]] std::vector<algen::WEdge> get_algen_edge_list() const {
                std::vector<algen::WEdge> result;
                result.reserve(edges.size());
                for (const auto &[a, b, _u, _v, w]: edges) {
                    result.emplace_back(a, b, w);
                }
                return result;
            }

            void remove_heavy_edges(const std::vector<bool> &is_light_edge) {
                assert(is_light_edge.size() == edges.size());
                i64 last = 0;
                for (size_t i = 0; i < edges.size(); i++) {
                    if (is_light_edge[i]) {
                        edges[last++] = edges[i];
                    }
                }
                assert(last % 2 == 0);
                assert(is_light_edge.size() % 2 == 0);
                assert(edges.size() % 2 == 0);
                std::cout << std::setprecision(4) << std::setw(5)
                          << 100. * (1. - ((double) last) / (double) edges.size()) << "% ("
                          << std::setw(8) << (edges.size() - last) / 2 << "/"
                          << std::setw(8) << edges.size() / 2 << ") of edges removed as heavy\n";
                edges.erase(edges.begin() + last, edges.end());
            }

            [[nodiscard]] ContractedEdgeListGraph get_subgraph_sample(std::mt19937_64 &gen) const {
                ContractedEdgeListGraph H;
                H.edges.reserve((size_t) (0.52 * (double) edges.size()));
                H.num_nodes = num_nodes;

                std::bernoulli_distribution dist;
                for (auto [a, b, _u, _v, w]: edges) {
                    // Make sure that the edge list H is symmetric, by ignoring edges with a > b and producing the edges
                    // in both directions if an edge with a < b is sampled.
                    if (a < b && dist(gen)) {
                        H.edges.push_back({a, b, a, b, w});
                        H.edges.push_back({b, a, b, a, w});
                    }
                }
                return H;
            }
        };

        void algorithm(ContractedEdgeListGraph &&G, algen::WEdgeList &mst, size_t level = 0) {
            std::cout << "level=" << level << "\n";

            // 1. If G is empty return an empty forest
            if (G.edges.size() <= 4) {
                for (const auto &edge: G.edges) {
                    mst.emplace_back(edge.u, edge.v, edge.weight);
                }
                return;
            }

            // 2. Create a contracted graph G' by running two successive Borůvka steps on G
            UnionFind<u64> union_find(G.num_nodes);

            //for (int i = 0; i < 3; ++i) {
            //    G = boruvka_step(std::move(G), union_find, mst);
            //    if (G.edges.empty()) return;
            //}
            auto G_prime = boruvka_steps(3, std::move(G), union_find, mst);
            if (G_prime.edges.empty()) return;

            G_prime.compress(union_find);

            if (G_prime.edges.size() <= 4) {
                for (const auto &edge: G_prime.edges) {
                    mst.emplace_back(edge.u, edge.v, edge.weight);
                }
                return;
            }

            // 3. Create a subgraph H by selecting each edge in G' with probability 1/2. Recursively apply the algorithm to H to get its minimum spanning forest F.
            auto H = G_prime.get_subgraph_sample(m_gen);

            algen::WEdgeList F;
            F.reserve(H.num_nodes - 1);

            algorithm(std::move(H), F, level + 1);


            // 4. Remove all F-heavy edges from G' (where F is the forest from step 3) using a linear time minimum spanning tree verification algorithm.
            auto is_light_edge = \
                algen::getEdgeClassifier().execute(G_prime.get_algen_edge_list(), F, G_prime.num_nodes);

            G_prime.remove_heavy_edges(is_light_edge);


            // 5. Recursively apply the algorithm to G' to get its minimum spanning forest.
            return algorithm(std::move(G_prime), mst, level + 1);

            // Output: The minimum spanning forest of G' and the contracted edges from the Borůvka steps
        }

        ContractedEdgeListGraph
        boruvka_step(ContractedEdgeListGraph &&G, UnionFind<algen::VertexId> &union_find, algen::WEdgeList &mst) {
            auto &cheapest_edges = m_cheapest_edges_scratchpad;
            assert(G.num_nodes == union_find.num_elements());
            assert(G.num_nodes <= cheapest_edges.size());

            auto is_better = [](const ContractedEdge &lhs, const ContractedEdge &rhs) {
                assert(rhs.weight == algen::WEIGHT_UNDEFINED || lhs.a == rhs.a);
                return std::tie(lhs.weight, lhs.b) < std::tie(rhs.weight, rhs.b);
            };

            for (auto &edge: G.edges) {
                auto &[a, b, _u, _v, _w] = edge;

                assert(a != b);
                assert(union_find.find(a) == a);
                assert(union_find.find(b) == b);

                if (is_better(edge, cheapest_edges[edge.a])) {
                    cheapest_edges[edge.a] = edge;
                }
            }

            u64 count = 0;
            for (u64 node = 0; node < G.num_nodes; ++node) {
                if (cheapest_edges[node].weight != algen::WEIGHT_UNDEFINED) {
                    auto &[a, b, u, v, w] = cheapest_edges[node];
                    if (union_find.do_union(a, b)) {
                        mst.emplace_back(u, v, w);
                        mst.emplace_back(v, u, w);
                    }
                    cheapest_edges[node].weight = algen::WEIGHT_UNDEFINED;
                } else {
                    count++;
                }
            }
            std::cout << std::setprecision(4) << std::setw(5) << 100. * (((double) count) / (double) G.num_nodes)
                      << "% ("
                      << std::setw(8) << (count) << "/"
                      << std::setw(8) << G.num_nodes << ") of nodes do not have a cheapest edge\n";

            i64 last = 0;
            for (size_t i = 0; i < G.edges.size(); i++) {
                auto &edge = G.edges[i];
                auto &[a, b, _u, _v, _w] = edge;
                a = union_find.find(a);
                b = union_find.find(b);
                if (G.edges[i].a != G.edges[i].b) {
                    G.edges[last++] = edge;
                }
            }
            std::cout << std::setprecision(4) << std::setw(5) << 100. * (1. - ((double) last) / (double) G.edges.size())
                      << "% ("
                      << std::setw(8) << (G.edges.size() - last) / 2 << "/"
                      << std::setw(8) << G.edges.size() / 2 << ") of edges removed as intra cluster edges\n";

            G.edges.erase(G.edges.begin() + last, G.edges.end());

            return G;
        }

        ContractedEdgeListGraph
        boruvka_steps(int num_steps, ContractedEdgeListGraph &&G, UnionFind<algen::VertexId> &union_find,
                      algen::WEdgeList &mst) {
            auto &cheapest_edges = m_cheapest_edges_scratchpad;
            auto &active_nodes = m_active_nodes_scratchpad;
            active_nodes.clear();

            assert(G.num_nodes == union_find.num_elements());
            assert(G.num_nodes <= cheapest_edges.size());

            auto is_better = [](const ContractedEdge &lhs, const ContractedEdge &rhs) {
                assert(rhs.weight == algen::WEIGHT_UNDEFINED || lhs.a == rhs.a);
                return std::tie(lhs.weight, lhs.b) < std::tie(rhs.weight, rhs.b);
            };

            /*
             * Call f for all inter cluster edges. Remove all intra cluster edges.
             */
            auto filter_and_map_edges = [&](auto f) {
                i64 last = 0;
                for (size_t i = 0; i < G.edges.size(); i++) {
                    auto &edge = G.edges[i];
                    auto &[a, b, _u, _v, _w] = edge;
                    a = union_find.find(a);
                    b = union_find.find(b);
                    if (G.edges[i].a != G.edges[i].b) {
                        G.edges[last++] = edge;

                        f(edge);
                    }
                }

                std::cout << std::setprecision(4) << std::setw(5)
                          << 100. * (1. - ((double) last) / (double) G.edges.size()) << "% ("
                          << std::setw(8) << (G.edges.size() - last) / 2 << "/"
                          << std::setw(8) << G.edges.size() / 2 << ") of edges removed as intra cluster edges"
                          << "\n";

                G.edges.erase(G.edges.begin() + last, G.edges.end());
            };

            /*
             * Call f for all nodes which have an edge.
             */
            auto map_all_nodes = [&](auto f) {
                u64 count = 0;
                for (u64 node = 0; node < G.num_nodes; ++node) {
                    if (cheapest_edges[node].weight != algen::WEIGHT_UNDEFINED) {
                        f(node);
                        cheapest_edges[node].weight = algen::WEIGHT_UNDEFINED;
                    } else {
                        count++;
                    }
                }
                std::cout << std::setprecision(4) << std::setw(5)
                          << 100. * ((double) (G.num_nodes - count) / (double) G.num_nodes)
                          << "% ("
                          << std::setw(8) << (G.num_nodes - count) << "/"
                          << std::setw(8) << G.num_nodes << ") of nodes were active\n";
            };

            /*
             * Call f for all nodes which have an edge and populate the active_nodes vector with them.
             */
            auto map_all_nodes_and_set_active = [&](auto f) {
                for (u64 node = 0; node < G.num_nodes; ++node) {
                    if (cheapest_edges[node].weight != algen::WEIGHT_UNDEFINED) {
                        f(node);
                        active_nodes.push_back(node);
                        cheapest_edges[node].weight = algen::WEIGHT_UNDEFINED;
                    }
                }
                std::cout << std::setprecision(4) << std::setw(5)
                          << 100. * (((double) active_nodes.size()) / (double) G.num_nodes)
                          << "% ("
                          << std::setw(8) << active_nodes.size() << "/"
                          << std::setw(8) << G.num_nodes << ") of nodes were active\n";
            };

            /*
             * Call f for all nodes which have an edge. Remove all other nodes from active_nodes.
             */
            auto filter_and_map_active_nodes = [&](auto f) {
                long long last = 0;
                for (u64 l = 0; l < active_nodes.size(); ++l) {
                    auto node = active_nodes[l];
                    if (cheapest_edges[node].weight != algen::WEIGHT_UNDEFINED) {
                        active_nodes[last++] = node;
                        f(node);
                        cheapest_edges[node].weight = algen::WEIGHT_UNDEFINED;
                    }
                }
                std::cout << std::setprecision(4) << std::setw(5)
                          << 100. * (((double) last) / (double) active_nodes.size())
                          << "% ("
                          << std::setw(8) << (last) << "/"
                          << std::setw(8) << active_nodes.size() << ") of nodes were active\n";

                active_nodes.erase(active_nodes.begin() + last, active_nodes.end());
            };


            // Cheapest edge selection.
            for (auto &edge: G.edges) {
                if (is_better(edge, cheapest_edges[edge.a])) {
                    cheapest_edges[edge.a] = edge;
                }
            }

            map_all_nodes([&](auto node) {
                auto &[a, b, u, v, w] = cheapest_edges[node];
                if (union_find.do_union(a, b)) {
                    mst.emplace_back(u, v, w);
                    mst.emplace_back(v, u, w);
                };
            });


            // Intra cluster edge filtering and cheapest edge selection.
            filter_and_map_edges([&](const auto &edge) {
                if (is_better(edge, cheapest_edges[edge.a])) {
                    cheapest_edges[edge.a] = edge;
                }
            });
            if (G.edges.empty()) return G;

            map_all_nodes_and_set_active([&](auto node) {
                auto &[a, b, u, v, w] = cheapest_edges[node];
                if (union_find.do_union(a, b)) {
                    mst.emplace_back(u, v, w);
                    mst.emplace_back(v, u, w);
                }
            });

            for (int k = 0; k < num_steps - 2; ++k) {
                // Intra cluster edge filtering and cheapest edge selection.
                filter_and_map_edges([&](const auto &edge) {
                    if (is_better(edge, cheapest_edges[edge.a])) {
                        cheapest_edges[edge.a] = edge;
                    }
                });
                if (G.edges.empty()) return G;

                filter_and_map_active_nodes([&](auto node) {
                    auto &[a, b, u, v, w] = cheapest_edges[node];
                    if (union_find.do_union(a, b)) {
                        mst.emplace_back(u, v, w);
                        mst.emplace_back(v, u, w);
                    }
                });
            }

            //filter_and_map_edges([](const auto&){});

            return G;
        }

        std::mt19937_64 m_gen;

        std::vector<ContractedEdge> m_cheapest_edges_scratchpad;
        std::vector<u64> m_active_nodes_scratchpad;
    };

}