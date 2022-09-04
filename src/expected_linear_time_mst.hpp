#pragma once

#include <utility>
#include <random>
#include <algorithm>
#include <iomanip>

#include "includes/definitions.hpp"
#include "includes/utils.hpp"

#include "datastructures/union_find.hpp"
#include "datastructures/robin_hood.h"

#include "datastructures/fast_reset_bitvector.hpp"

namespace js {
    class ExpectedLinearTimeMST {
        using i64 = int64_t;
        using u64 = uint64_t;
        using u32 = uint32_t;
        using VertexId = algen::VertexId;
        using EdgeIdx = algen::EdgeIdx;
        using Weight = algen::Weight;
        static constexpr VertexId VERTEXID_UNDEFINED = algen::VERTEXID_UNDEFINED;

        static constexpr bool LOG_INFO = false;

        struct PairHash {
            template<typename A, typename B>
            std::size_t operator()(const std::pair<A, B> &pair) const {
                return std::hash<A>{}(pair.first) << 32 ^ std::hash<B>{}(pair.second);
            }
        };

    public:

        explicit ExpectedLinearTimeMST(size_t num_boruvka_phases = 3) : m_gen(std::random_device{}()),
                                                                        m_num_boruvka_phases(num_boruvka_phases) {}

        algen::WEdgeList operator()(const algen::WEdgeList &edge_list,
                                    const algen::VertexId num_vertices) {
            using namespace algen;

            auto G = ContractedEdgeListGraph::from_algen_edge_list(edge_list, num_vertices);
            algen::WEdgeList mst;
            mst.reserve(num_vertices - 1);

            m_cheapest_edges_scratchpad.resize(G.num_nodes, INVALID_CONTRACTED_EDGE);
            m_active_nodes_scratchpad.reserve((size_t) (0.27 * (double) G.num_nodes));

            algorithm(std::move(G), mst);

            return mst;
        }

    private:
        struct ContractedEdge {
            VertexId a;
            VertexId b;
            VertexId u;
            VertexId v;
            Weight weight;

            friend std::ostream &operator<<(std::ostream &os, const ContractedEdge &edge) {
                return os << "(" << edge.a << "," << edge.b << ") " << algen::WEdge{edge.u, edge.v, edge.weight};
            }
        };

        static_assert(std::is_trivially_constructible_v<ContractedEdge>);
        static_assert(std::is_trivially_destructible_v<ContractedEdge>);
        static_assert(std::is_trivially_constructible_v<algen::WEdge>);
        static_assert(std::is_trivially_destructible_v<algen::WEdge>);

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

            void compress(UnionFind<algen::VertexId> &union_find) {
                robin_hood::unordered_map<VertexId, VertexId> new_vertex_ids;
                new_vertex_ids.reserve(64);

                auto new_num_nodes = 0;

                auto get_mapped_vertex_id = [&](VertexId u) -> VertexId {
                    assert(union_find.is_non_singleton_representative(u));
                    auto [it, inserted] = new_vertex_ids.try_emplace(u, new_num_nodes);
                    new_num_nodes += inserted ? 1 : 0;
                    return it->second;
                };

                robin_hood::unordered_map<std::pair<VertexId, VertexId>, algen::WEdge, PairHash> inter_cluster_edges;
                inter_cluster_edges.reserve(64);

                for (auto [a, b, u, v, w]: edges) {
                    a = union_find.find(a);
                    b = union_find.find(b);
                    if (a != b) {
                        assert(a == union_find.find(a));
                        assert(b == union_find.find(b));
                        assert(a != b);

                        a = get_mapped_vertex_id(a);
                        b = get_mapped_vertex_id(b);

                        if (a > b) std::swap(a, b);
                        assert(a != b);
                        assert(a < num_nodes);
                        assert(b < num_nodes);
                        algen::WEdge uvw = {u, v, w};
                        auto [it, inserted] = inter_cluster_edges.emplace(std::pair{a, b}, uvw);
                        if (!inserted) {
                            if (w < it->second.weight) {
                                it->second = uvw;
                            }
                        }
                    }
                }

                if constexpr (LOG_INFO)
                    std::cout << std::setprecision(4) << std::setw(5)
                              << 100. * (1. - (double) new_num_nodes / (double) num_nodes) << "% ("
                              << std::setw(8) << num_nodes - new_num_nodes << "/"
                              << std::setw(8) << num_nodes << ") of nodes removed by compression\n";

                num_nodes = new_num_nodes;

                if constexpr (LOG_INFO)
                    std::cout << std::setprecision(4) << std::setw(5)
                              << 100. * (1. - (2. * (double) inter_cluster_edges.size()) / (double) edges.size())
                              << "% ("
                              << std::setw(8) << edges.size() / 2 - inter_cluster_edges.size() << "/"
                              << std::setw(8) << edges.size() / 2
                              << ") of edges removed as intra cluster and multi edges\n";

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
                    edges[last] = edges[i];
                    last += is_light_edge[i] ? 1 : 0;
                }

                assert(last % 2 == 0);
                assert(is_light_edge.size() % 2 == 0);
                assert(edges.size() % 2 == 0);

                if constexpr (LOG_INFO)
                    std::cout << std::setprecision(4) << std::setw(5)
                              << 100. * (1. - ((double) last) / (double) edges.size()) << "% ("
                              << std::setw(8) << (edges.size() - last) / 2 << "/"
                              << std::setw(8) << edges.size() / 2 << ") of edges removed as heavy\n";

                edges.erase(edges.begin() + last, edges.end());
            }

            [[nodiscard]] ContractedEdgeListGraph get_subgraph_sample(std::minstd_rand &gen) const {
                auto t0 = std::chrono::steady_clock::now();
                ContractedEdgeListGraph H;
                H.edges.reserve((size_t) (0.52 * (double) edges.size()));
                H.num_nodes = num_nodes;

                std::uniform_int_distribution<u64> dist(0, std::numeric_limits<u64>::max());
                u64 r = dist(gen);
                int rw = 64;

                // Note: branchless version with
                // unconditioned writing to H.edges[last] and H.edges[last+1],
                // advancing last by 2 * ((a < b ? 1 : 0) & (r & 1)) and
                // doing 64 iterations until r = dist(gen)
                // was not faster than branching version.

                for (auto [a, b, _u, _v, w]: edges) {
                    // Make sure that the edge list H is symmetric, by ignoring edges with a > b and producing the edges
                    // in both directions if an edge with a < b is sampled.
                    if (a < b) {
                        if ((r & 1) != 0) {
                            H.edges.push_back({a, b, a, b, w});
                            H.edges.push_back({b, a, b, a, w});
                        }
                        r >>= 1;
                        rw = (rw + 1) & 63;
                        if (rw == 0) {
                            r = dist(gen);
                        }
                    }
                }
                return H;
            }
        };

        void algorithm(ContractedEdgeListGraph &&G, algen::WEdgeList &mst, size_t level = 0) {
            if constexpr (LOG_INFO)
                std::cout << "level=" << level << "\n";

            // 1. If G has at most 2 edges return all edges as forest
            auto tiny_graph_mst = [](const auto &G, auto &mst) {
                if (G.edges.size() <= 4) {
                    for (const auto &edge: G.edges) {
                        mst.emplace_back(edge.u, edge.v, edge.weight);
                    }
                    return true;
                }
                return false;
            };
            if (tiny_graph_mst(G, mst)) return;

            // 2. Create a contracted graph G' by running three successive Borůvka steps on G
            auto &union_find = m_union_find_scratchpad;
            union_find.reset(G.num_nodes);

            auto G_prime = boruvka_steps(m_num_boruvka_phases, std::move(G), union_find, mst);
            if (G_prime.edges.empty()) return;

            G_prime.compress(union_find);

            if (tiny_graph_mst(G_prime, mst)) return;

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
        boruvka_steps(int num_steps, ContractedEdgeListGraph &&G, UnionFind<algen::VertexId> &union_find,
                      algen::WEdgeList &mst) {
            auto &cheapest_edges = m_cheapest_edges_scratchpad;
            auto &active_nodes = m_active_nodes_scratchpad;
            active_nodes.clear();

            assert(G.num_nodes == union_find.num_elements());
            assert(G.num_nodes <= cheapest_edges.size());

            constexpr auto is_better = [](const ContractedEdge &lhs, const ContractedEdge &rhs) constexpr noexcept {
                assert(rhs.weight == algen::WEIGHT_UNDEFINED || lhs.a == rhs.a);
                return std::tie(lhs.weight, lhs.b) < std::tie(rhs.weight, rhs.b);
            };

            /*
             * Call f for all inter cluster edges. Remove all intra cluster edges.
             */
            auto filter_and_map_edges = [&](auto f) constexpr noexcept {
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

                if constexpr (LOG_INFO)
                    std::cout << std::setprecision(4) << std::setw(5)
                              << 100. * (1. - ((double) last) / (double) G.edges.size()) << "% ("
                              << std::setw(8) << (G.edges.size() - last) / 2 << "/"
                              << std::setw(8) << G.edges.size() / 2 << ") of edges removed as intra cluster edges"
                              << "\n";

                G.edges.erase(G.edges.begin() + last, G.edges.end());
            };

            /*
             * Call f for all nodes which have at least one incident edge.
             */
            auto map_all_nodes = [&](auto f) constexpr noexcept {
                size_t count = 0;
                for (VertexId node = 0; node < G.num_nodes; ++node) {
                    if (cheapest_edges[node].weight != algen::WEIGHT_UNDEFINED) {
                        f(node);
                        cheapest_edges[node].weight = algen::WEIGHT_UNDEFINED;
                    } else if constexpr (LOG_INFO) {
                        count++;
                    }
                }

                if constexpr (LOG_INFO)
                    std::cout << std::setprecision(4) << std::setw(5)
                              << 100. * ((double) (G.num_nodes - count) / (double) G.num_nodes)
                              << "% ("
                              << std::setw(8) << (G.num_nodes - count) << "/"
                              << std::setw(8) << G.num_nodes << ") of nodes were active\n";
            };

            /*
             * Call f for all nodes which have at least one incident edge and populate the active_nodes vector with them.
             */
            auto map_all_nodes_and_set_active = [&](auto f) constexpr noexcept {
                for (VertexId node = 0; node < G.num_nodes; ++node) {
                    if (cheapest_edges[node].weight != algen::WEIGHT_UNDEFINED) {
                        f(node);
                        active_nodes.push_back(node);
                        cheapest_edges[node].weight = algen::WEIGHT_UNDEFINED;
                    }
                }

                if constexpr (LOG_INFO)
                    std::cout << std::setprecision(4) << std::setw(5)
                              << 100. * (((double) active_nodes.size()) / (double) G.num_nodes)
                              << "% ("
                              << std::setw(8) << active_nodes.size() << "/"
                              << std::setw(8) << G.num_nodes << ") of nodes were active\n";
            };

            /*
             * Call f for all nodes which have at least one incident edge. Remove all other nodes from active_nodes.
             */
            auto filter_and_map_active_nodes = [&](auto f) constexpr noexcept {
                size_t last = 0;
                for (size_t i = 0; i < active_nodes.size(); ++i) {
                    auto node = active_nodes[i];
                    if (cheapest_edges[node].weight != algen::WEIGHT_UNDEFINED) {
                        active_nodes[last++] = node;
                        f(node);
                        cheapest_edges[node].weight = algen::WEIGHT_UNDEFINED;
                    }
                }

                if constexpr (LOG_INFO)
                    std::cout << std::setprecision(4) << std::setw(5)
                              << 100. * (((double) last) / (double) active_nodes.size())
                              << "% ("
                              << std::setw(8) << (last) << "/"
                              << std::setw(8) << active_nodes.size() << ") of nodes were active\n";

                active_nodes.erase(active_nodes.begin() + last, active_nodes.end());
            };

            /*
             * Call f for all nodes which have at least one incident edge. Don't update active_nodes.
             */
            auto map_active_nodes = [&](auto f) constexpr noexcept {
                size_t count = 0;
                for (auto node: active_nodes) {
                    if (cheapest_edges[node].weight != algen::WEIGHT_UNDEFINED) {
                        f(node);
                        cheapest_edges[node].weight = algen::WEIGHT_UNDEFINED;
                    } else if constexpr (LOG_INFO) {
                        count++;
                    }
                }

                if constexpr (LOG_INFO)
                    std::cout << std::setprecision(4) << std::setw(5)
                              << 100. * ((double) (active_nodes.size() - count) / (double) active_nodes.size())
                              << "% ("
                              << std::setw(8) << (active_nodes.size() - count) << "/"
                              << std::setw(8) << active_nodes.size() << ") of nodes were active\n";
            };


            // First Boruvka Step

            // Cheapest edge selection
            // All edges are guaranteed to be inter-cluster edges
            for (auto &edge: G.edges) {
                if (is_better(edge, cheapest_edges[edge.a])) {
                    cheapest_edges[edge.a] = edge;
                }
            }

            // Nearly all nodes have at least one incident edge for most input graphs
            map_all_nodes([&](auto node) constexpr noexcept {
                auto &[a, b, u, v, w] = cheapest_edges[node];
                if (union_find.do_union(a, b)) {
                    mst.emplace_back(u, v, w);
                    mst.emplace_back(v, u, w);
                };
            });


            // Second Boruvka Step

            // Intra cluster edge filtering and cheapest edge selection
            filter_and_map_edges([&](const auto &edge) constexpr noexcept {
                if (is_better(edge, cheapest_edges[edge.a])) {
                    cheapest_edges[edge.a] = edge;
                }
            });
            if (G.edges.empty()) return G;

            // About 25% of nodes have at least one incident edge after the first boruvka step
            map_all_nodes_and_set_active([&](auto node) constexpr noexcept {
                auto &[a, b, u, v, w] = cheapest_edges[node];
                if (union_find.do_union(a, b)) {
                    mst.emplace_back(u, v, w);
                    mst.emplace_back(v, u, w);
                }
            });

            // Third to second to last Boruvka Step
            for (int k = 0; k < num_steps - 3; ++k) {
                // Intra cluster edge filtering and cheapest edge selection
                filter_and_map_edges([&](const auto &edge) {
                    if (is_better(edge, cheapest_edges[edge.a])) {
                        cheapest_edges[edge.a] = edge;
                    }
                });
                if (G.edges.empty()) return G;

                filter_and_map_active_nodes([&](auto node) noexcept {
                    auto &[a, b, u, v, w] = cheapest_edges[node];
                    if (union_find.do_union(a, b)) {
                        mst.emplace_back(u, v, w);
                        mst.emplace_back(v, u, w);
                    }
                });
            }

            // Last Boruvka Step

            // Intra cluster edge filtering and cheapest edge selection
            filter_and_map_edges([&](const auto &edge) {
                if (is_better(edge, cheapest_edges[edge.a])) {
                    cheapest_edges[edge.a] = edge;
                }
            });
            if (G.edges.empty()) return G;

            // active_nodes won't be used later, so there is no need to update it
            map_active_nodes([&](auto node) noexcept {
                auto &[a, b, u, v, w] = cheapest_edges[node];
                if (union_find.do_union(a, b)) {
                    mst.emplace_back(u, v, w);
                    mst.emplace_back(v, u, w);
                }
            });

            return G;
        }

        std::minstd_rand m_gen;

        size_t m_num_boruvka_phases;

        UnionFind<VertexId> m_union_find_scratchpad;
        std::vector<ContractedEdge> m_cheapest_edges_scratchpad;
        std::vector<VertexId> m_active_nodes_scratchpad;
    };

}