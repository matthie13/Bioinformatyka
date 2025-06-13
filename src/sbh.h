#ifndef SBH_H
#define SBH_H

#include "config.h"

struct Edge {
    int target;    // target k-mer index
    int weight;       // weight = k - overlap_length
    int overlap;      // overlap length

    Edge(int t, int w, int o) : target(t), weight(w), overlap(o) {}
};

struct SBHInstance {
    int n;                    // Target sequence length
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
                                const std::vector<int>& trail,
                                int n,
                                const AdjacencyList& adj_matrix);

std::optional<int> edge_weight(const AdjacencyList& adj_matrix,
                               const std::vector<std::string>& kmers,
                               int u, int v, int k);

std::vector<int> find_eulerian_path(const AdjacencyList& adj_matrix);

bool find_hamiltonian_path_helper(const AdjacencyList &adj_matrix,
                                  std::vector<int> &path,
                                  std::vector<bool> &visited,
                                  int current,
                                  int target_length);

std::vector<int> find_hamiltonian_path(const AdjacencyList& adj_matrix);

bool is_connected(const AdjacencyList& adj_matrix);

std::string sequencing_by_hybridization(const std::vector<std::string>& kmers,
                                        int k,
                                        int target_length);

std::string sequencing_by_hybridization_with_start(const std::vector<std::string>& kmers,
                                        int k,
                                        int target_length,
                                        const std::string& start_oligo);

std::vector<int> find_path_from_start(const AdjacencyList& adj,
                                        int start,
                                        int total_nodes);

bool dfs_path(const AdjacencyList& adj,
              int current,
              std::vector<bool>& visited,
              std::vector<int>& path,
              int total_nodes);

std::vector<std::string> generate_kmers(const std::string& sequence, int k);

std::pair<int, std::vector<int>> dijkstra_shortest_path(const AdjacencyList& adj_matrix,
                                                              const std::vector<int>& used_count,
                                                              const std::vector<int>& repeat_limits,
                                                              int start,
                                                              const std::vector<int>& targets
                                                              );

#endif //SBH_H
