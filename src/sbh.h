#ifndef SBH_H
#define SBH_H

#include "config.h"

struct Edge {
    size_t target;    // target k-mer index
    int weight;       // weight = k - overlap_length
    int overlap;      // overlap length

    Edge(size_t t, int w, int o) : target(t), weight(w), overlap(o) {}
};

struct SBHInstance {
    size_t n;                    // Target sequence length
    int k;                       // K-mer length
    std::string start_oligo;     // Starting k-mer
    int neg_errors;              // Number of negative errors
    bool has_repeats;            // Whether sequence has repeats
    int pos_errors;              // Number of positive errors
    std::vector<std::string> spectrum;  // K-mer spectrum
};

using AdjacencyList = std::vector<std::vector<Edge>>;
AdjacencyList build_adjacency_matrix(const std::vector<std::string>& kmers, int k);

std::string reconstruct_sequence(const std::vector<std::string>& kmers,
                                int k,
                                const std::vector<size_t>& trail,
                                size_t n,
                                const AdjacencyList& adj_matrix);

std::optional<int> edge_weight(const AdjacencyList& adj_matrix,
                               const std::vector<std::string>& kmers,
                               size_t u, size_t v, int k);

std::vector<size_t> find_eulerian_path(const AdjacencyList& adj_matrix);

bool find_hamiltonian_path_helper(const AdjacencyList &adj_matrix,
                                  std::vector<size_t> &path,
                                  std::vector<bool> &visited,
                                  size_t current,
                                  size_t target_length);

std::vector<size_t> find_hamiltonian_path(const AdjacencyList& adj_matrix);

bool is_connected(const AdjacencyList& adj_matrix);

std::string sequencing_by_hybridization(const std::vector<std::string>& kmers,
                                        int k,
                                        size_t target_length);

std::string sequencing_by_hybridization_with_start(const std::vector<std::string>& kmers,
                                        int k,
                                        size_t target_length,
                                        const std::string& start_oligo);

std::vector<size_t> find_path_from_start(const AdjacencyList& adj,
                                        size_t start,
                                        size_t total_nodes);

bool dfs_path(const AdjacencyList& adj,
              size_t current,
              std::vector<bool>& visited,
              std::vector<size_t>& path,
              size_t total_nodes);

std::vector<std::string> generate_kmers(const std::string& sequence, int k);

std::pair<size_t, std::vector<size_t>> dijkstra_shortest_path(const AdjacencyList& adj_matrix,
                                                              const std::vector<int>& used_count,
                                                              const std::vector<int>& repeat_limits,
                                                              size_t start,
                                                              const std::vector<size_t>& targets
                                                              );

#endif //SBH_H
