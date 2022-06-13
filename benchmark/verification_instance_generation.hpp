#pragma once


#include <stack>

namespace benchmark {

    // Generates instances for an MST verification algorithm by accepting a correct MST and for a given number k
    // changing k edges s.t. the resulting spanning tree is no longer an MST.
    class CorruptedMSTGenerator {

    public:

        CorruptedMSTGenerator() = default;

        // Computes an MST for the given graph to be used by subsequent calls to generate_corrupted_mst().
        void preprocess(const algen::WEdgeList& graph_edge_list, const algen::WEdgeList& msf_edge_list, const algen::VertexId num_vertices, const std::size_t seed = 1) {

            num_vertices_ = num_vertices;
            seed_ = seed;
            graph_edges_ = graph_edge_list;
            make_inefficient_adjacency_structure(graph_edges_, first_graph_out_edge_, num_vertices);

            msf_edges_ = msf_edge_list;
            make_inefficient_adjacency_structure(msf_edges_, first_msf_out_edge_, num_vertices);
        }

        // Generate an ST by changing num_changed_edges edges in the MST generated in the last call to configure().
        // If num_changed_edges == 0, this returns a correct MST; otherwise it returns an ST that is not an MST.
        // Assumes that the graph is sufficiently dense to support randomly choosing edges until enough MSF-heavy edges
        // have been found to replace MSF edges.
        algen::WEdgeList generate_corrupted_mst(const std::size_t num_changed_edges, const bool print_status = false) const {

            if (print_status) std::cout << "editing " << num_changed_edges << " edges of the MSF ... " << std::flush;

            auto sf_edges = msf_edges_;
            auto first_sf_out_edge = first_msf_out_edge_;

            assert(algen::edge_list_format_check(sf_edges, num_vertices_).first);
            for (std::size_t i = 0; i < num_changed_edges; ++i) {
                find_random_cycle_and_edit(sf_edges, first_sf_out_edge);
            }
            assert(algen::edge_list_format_check(sf_edges, num_vertices_).first);

            if (print_status) std::cout << " done." << std::endl;

            return sf_edges;
        }

    private:

        bool mst_contains_edge(const algen::VertexId tail, const algen::VertexId head) const {
            for (auto i = first_msf_out_edge_[tail]; i < first_msf_out_edge_[tail + 1]; ++i) {
                assert(msf_edges_[i].tail == tail);
                if (msf_edges_[i].head == head)
                    return true;
            }
            return false;
        }

        void find_random_cycle_and_edit(algen::WEdgeList& sf_edges, std::vector<algen::EdgeIdx>& first_sf_out_edge) const {

            using namespace algen;

            auto sf_contains_edge = [&] (const algen::VertexId tail, const algen::VertexId head) {
                for (auto i = first_sf_out_edge[tail]; i < first_sf_out_edge[tail + 1]; ++i) {
                    assert(sf_edges[i].tail == tail);
                    if (sf_edges[i].head == head)
                        return true;
                }
                return false;
            };

            // Find a cycle using a DFS s.t. the non-SF edge is not an MSF edge that was previously removed and at
            // least one of the SF edges was initially part of the MSF.

            static constexpr EdgeIdx INVALID_EDGE_IDX = std::numeric_limits<VertexId>::max();
            std::vector<VertexId> idx_of_sf_edge_to(num_vertices_);
            std::vector<VertexId> depth(num_vertices_);
            std::vector<VertexId> next_sf_edge_offset(num_vertices_);
            std::stack<VertexId, std::vector<VertexId>> active;

            std::mt19937 gen(seed_);
            std::uniform_int_distribution<std::size_t> uniform_vertex_dist(0, num_vertices_ - 1);
            while (true) {
                // Draw random start vertex for DFS
                const VertexId s = uniform_vertex_dist(gen);

                std::vector<bool> marked(num_vertices_, false);
                std::fill(idx_of_sf_edge_to.begin(), idx_of_sf_edge_to.end(), INVALID_EDGE_IDX);
                std::fill(depth.begin(), depth.end(), VERTEXID_UNDEFINED);
                std::fill(next_sf_edge_offset.begin(), next_sf_edge_offset.end(), 0);
                while(!active.empty()) active.pop();

                // Find a graph edge e that is not in the SF and not in the MSF but that closes a cycle in the SF.
                active.push(s);
                marked[s] = true;
                depth[s] = 0;
                EdgeIdx graph_edge_idx_of_e = INVALID_EDGE_IDX;
                bool found_e = false;
                while (!active.empty()) {
                    const auto v = active.top();

                    for (auto graph_edge_idx = first_graph_out_edge_[v]; graph_edge_idx < first_graph_out_edge_[v + 1]; ++graph_edge_idx) {
                        const auto& e = graph_edges_[graph_edge_idx];
                        assert(e.tail == v);
                        if (marked[e.head] && !mst_contains_edge(e.tail, e.head) && !sf_contains_edge(e.tail, e.head)) {
                            graph_edge_idx_of_e = graph_edge_idx;
                            found_e = true;
                            break;
                        }
                    }
                    if (found_e) break;

                    const auto next_sf_edge = first_sf_out_edge[v] + next_sf_edge_offset[v];
                    if (next_sf_edge != first_sf_out_edge[v + 1]) {
                        const VertexId head = sf_edges[next_sf_edge].head;
                        ++next_sf_edge_offset[v];
                        if (!marked[head]) {
                            idx_of_sf_edge_to[head] = next_sf_edge;
                            active.push(head);
                            marked[head] = true;
                            depth[head] = depth[v] + 1;
                        }
                    } else {
                        active.pop();
                    }
                }
                assert(!active.empty());

                // Attempt to find an edge e' in the cycle closed by e s.t. e' is in the original MST and w(e') < w(e)

                // Wlog. let e.head have greater depth than e.tail. Then, all edges e' leading to e.head with
                // depth[e'.head] > e.tail are on the cycle.
                const WEdge& e = graph_edges_[graph_edge_idx_of_e];
                VertexId h = e.head;
                VertexId t = e.tail;

                auto check_edge = [&](const WEdge& edge){
                    return mst_contains_edge(edge.tail, edge.head) && edge.weight < e.weight;
                };

                while (depth[h] > depth[e.tail]) {
                    const auto& edge_to_h = sf_edges[idx_of_sf_edge_to[h]];
                    if (check_edge(edge_to_h)) {
                        replace_edge(edge_to_h, e, sf_edges, first_sf_out_edge);
                        return;
                    }
                    h = edge_to_h.tail;
                }

                while (depth[t] > depth[e.head]) {
                    const auto& edge_to_t = sf_edges[idx_of_sf_edge_to[t]];
                    if (check_edge(edge_to_t)) {
                        replace_edge(edge_to_t, e, sf_edges, first_sf_out_edge);
                        return;
                    }
                    t = edge_to_t.tail;
                }

                // All edges on the path to the deepest common ancestor of e.head and e.tail are on the cycle
                assert(depth[h] == depth[t] && depth[h] == std::min(depth[e.head], depth[e.tail]));
                while (h != t && depth[h] > 0) {
                    assert(depth[h] == depth[t]);
                    const auto& edge_to_h = sf_edges[idx_of_sf_edge_to[h]];
                    if (check_edge(edge_to_h)) {
                        replace_edge(edge_to_h, e, sf_edges, first_sf_out_edge);
                        return;
                    }
                    h = edge_to_h.tail;

                    const auto& edge_to_t = sf_edges[idx_of_sf_edge_to[t]];
                    if (check_edge(edge_to_t)) {
                        replace_edge(edge_to_t, e, sf_edges, first_sf_out_edge);
                        return;
                    }
                    t = edge_to_t.tail;
                }

                // If we reached the lowest common ancestor of e.head and e.tail in the DFS, then the cycle did not
                // contain an edge that was suitable for e to replace. In this case, we try again with a different
                // starting vertex.
            }

        }

        void replace_edge(const algen::WEdge to_remove, const algen::WEdge to_insert,
                          algen::WEdgeList& edges, std::vector<algen::EdgeIdx>& first_out_edge) const {

            // Remove to_remove in forward and reverse direction
            sorted_remove_edge(to_remove.tail, to_remove.head, edges, first_out_edge);
            sorted_remove_edge(to_remove.head, to_remove.tail, edges, first_out_edge);

            // Insert to_insert in forward and reverse direction
            sorted_insert_edge(to_insert.tail, to_insert.head, to_insert.weight, edges, first_out_edge);
            sorted_insert_edge(to_insert.head, to_insert.tail, to_insert.weight, edges, first_out_edge);

            assert(check_adjacency_structure_correct(edges, first_out_edge));
        }

        void sorted_remove_edge(const algen::VertexId tail, const algen::VertexId head,
                                algen::WEdgeList& edges, std::vector<algen::EdgeIdx>& first_out_edge) const {
            algen::EdgeIdx edge_idx = first_out_edge[tail];
            while (edge_idx < first_out_edge[tail + 1] && edges[edge_idx].head != head) {
                ++edge_idx;
            }
            assert(edge_idx < first_out_edge[tail + 1] && edges[edge_idx].head == head);
            edges.erase(edges.begin() + edge_idx);
            for (algen::VertexId v = tail + 1; v < num_vertices_ + 1; ++v) {
                --first_out_edge[v];
            }
        }

        void sorted_insert_edge(const algen::VertexId tail, const algen::VertexId head, const algen::Weight weight,
                                algen::WEdgeList& edges, std::vector<algen::EdgeIdx>& first_out_edge) const {
            algen::EdgeIdx edge_idx = first_out_edge[tail];
            while (edge_idx < first_out_edge[tail + 1] && head >= edges[edge_idx].head) {
                assert(head != edges[edge_idx].head);
                ++edge_idx;
            }
            edges.insert(edges.begin() + edge_idx, {tail, head, weight});
            for (algen::VertexId v = tail + 1; v < num_vertices_ + 1; ++v) {
                ++first_out_edge[v];
            }
        }

        bool check_adjacency_structure_correct(const algen::WEdgeList& edges, const std::vector<algen::EdgeIdx>& first_out_edge) const {

            if (first_out_edge.size() != num_vertices_ + 1 || first_out_edge[num_vertices_] != edges.size())
                return false;

            if (!std::is_sorted(edges.begin(), edges.end(), algen::TailHeadOrder<algen::WEdge>())) return false;
            for (algen::VertexId tail = 0; tail < num_vertices_; ++tail) {
                for (algen::VertexId edge_idx = first_out_edge[tail]; edge_idx < first_out_edge[tail + 1]; ++edge_idx) {
                    if (edges[edge_idx].tail != tail) return false;
                    if (edge_idx < first_out_edge[tail + 1] - 1 && edges[edge_idx].head > edges[edge_idx + 1].head)
                        return false;
                }
            }
            return true;
        }

        algen::WEdgeList graph_edges_;
        std::vector<algen::EdgeIdx> first_graph_out_edge_;

        algen::VertexId num_vertices_;
        algen::WEdgeList msf_edges_;
        std::vector<algen::EdgeIdx> first_msf_out_edge_;

        std::size_t seed_;
    };

//    // Generates instances for an MST verification algorithm by constructing spanning trees using a depth first search.
//    // In almost all cases, these STs will not be MSTs, i.e. the STs are almost always "no"-instances for the
//    // MST-verification algorithm.
//    class DFSSpanningTreeGenerator {
//
//        WEdgeList generate_st_using_dfs(const WEdgeList& graph_edge_list, const VertexId num_vertices) const {
//
//        }
//
//    };

} // end namespace benchmark
