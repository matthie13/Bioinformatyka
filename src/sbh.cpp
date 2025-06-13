#include "sbh.h"

#include <queue>
#include <unordered_set>


AdjacencyList build_adjacency_matrix(const std::vector<std::string>& kmers, int k) {
    for (const auto& kmer : kmers) {
        if (kmer.length() == static_cast<int>(k)) continue;
        std::cerr << "ERROR: Invalid k-mer length: " << kmer << "\n";
        exit(EXIT_FAILURE);
    }

    const int num_kmers = kmers.size();

    std::unordered_map<int, std::unordered_map<std::string, std::vector<int>>> prefix_map;

    for (int j = 0; j < num_kmers; j++) {
        const std::string& seq = kmers[j];
        for (int o = 1; o < k; o++) {
            std::string prefix = seq.substr(0, o);
            prefix_map[o][prefix].push_back(j);
        }
    }

    AdjacencyList adj_matrix(num_kmers);

    for (int i = 0; i < num_kmers; i++) {
        const std::string& ui = kmers[i];

        for (int o = k - 1; o >= 1; o--) {
            if (ui.length() >= static_cast<int>(o)) {
                std::string suffix = ui.substr(ui.length() - o);

                auto prefix_map_for_overlap = prefix_map.find(o);
                if (prefix_map_for_overlap != prefix_map.end()) {
                    auto suffix_entry = prefix_map_for_overlap->second.find(suffix);
                    if (suffix_entry != prefix_map_for_overlap->second.end()) {
                        for (int j : suffix_entry->second) {
                            if (j != i) {
                                int weight = k - o;
                                adj_matrix[i].emplace_back(j, weight, o);
                            }
                        }
                    }
                }
            }
        }
    }
    return adj_matrix;
}

std::string reconstruct_sequence(const std::vector<std::string>& kmers,
                                int k,
                                const std::vector<int>& trail,
                                int n,
                                const AdjacencyList& adj_matrix) {

    if (trail.empty() || kmers.empty()) return "";

    if (trail[0] >= kmers.size()) {
        std::cerr << "ERROR: Invalid index in trail: " << trail[0] << "\n";
        return "";
    }

    std::string seq = kmers[trail[0]];
    int seq_len = k;

    for (int i = 1; i < trail.size(); i++) {
        int u = trail[i - 1];
        int v = trail[i];

        if (u >= kmers.size() || v >= kmers.size()) {
            std::cerr << "ERROR: Invalid trail indices: " << u << ", " << v << "\n";
            break;
        }

        int overlap = 0;
        for (const auto edge : adj_matrix[u]) {
            if (edge.target != v) continue;
            overlap = edge.overlap; // Use stored overlap
            break;
        }

        int add_len = k - overlap;
        if (seq_len + add_len > n) {
            int to_add = n - seq_len;
            if (overlap + to_add <= kmers[v].size()) {
                seq += kmers[v].substr(overlap, to_add);
            } else {
                seq += kmers[v].substr(overlap);
            }
            break;
        }

        seq += kmers[v].substr(overlap);
        seq_len += add_len;

        if (seq_len >= n) break;
    }
    if (seq.length() > n) seq.resize(n);

    return seq;
}

std::optional<int> edge_weight(const AdjacencyList& adj_matrix,
                               const std::vector<std::string>& kmers,
                               int u, int v, int k) {

    for (const auto edge : adj_matrix[u]) {
        if (edge.target != v) continue;
        return edge.weight;
    }
    return std::nullopt;
}

std::vector<int> find_eulerian_path(const AdjacencyList& adj_matrix) {
    int n = adj_matrix.size();
    if (n == 0) return {};

    std::vector<int> out_degree(n, 0);
    std::vector<int> in_degree(n, 0);

    for (int i = 0; i < n; i++) {
        out_degree[i] = adj_matrix[i].size();
        for (const auto edge : adj_matrix[i]) {
            in_degree[edge.target]++;
        }
    }

    int start = 0;
    int start_nodes = 0;
    int end_nodes = 0;

    for (int i = 0; i < n; i++) {
        if (out_degree[i] == in_degree[i] + 1) {
            // Node with one more outgoing edge (start node)
            start = i;
            start_nodes++;
        } else if (in_degree[i] == out_degree[i] + 1) {
            // Node with one more incoming edge (end node)
            end_nodes++;
        } else if (out_degree[i] != in_degree[i]) {
            // Invalid difference - no Eulerian path exists
            return {};
        }
    }

    if (start_nodes == 0 && end_nodes == 0) {
        for (int i = 0; i < n; ++i) {
            if (out_degree[i] > 0) {
                start = i;
                break;
            }
        }
    } else if (start_nodes != 1 || end_nodes != 1) {
        return {};
    }

    std::vector<std::vector<Edge>> temp_adj_matrix = adj_matrix;
    std::vector<int> path;
    std::vector<int> circuit;
    circuit.push_back(start);

    int current = start;
    while (!circuit.empty()) {
        if (!temp_adj_matrix[current].empty()) {
            circuit.push_back(current);
            int next = temp_adj_matrix[current].back().target;
            temp_adj_matrix[current].pop_back();
            current = next;
        } else {
            path.push_back(current);
            current = circuit.back();
            circuit.pop_back();
        }
    }
    std::ranges::reverse(path);
    return path;
}

bool find_hamiltonian_path_helper(const AdjacencyList &adj_matrix, std::vector<int> &path,
                                  std::vector<bool> &visited, int current, int target_length) {

    if (path.size() == target_length) return true;

    for (const auto edge : adj_matrix[current]) {
        int next = edge.target;
        if (!visited[next]) {
            visited[next] = true;
            path.push_back(next);

            if (find_hamiltonian_path_helper(adj_matrix, path, visited, next, target_length)) return true;

            path.pop_back();
            visited[next] = false;
        }
    }
    return false;
}

std::vector<int> find_hamiltonian_path(const AdjacencyList& adj_matrix) {
    int n = adj_matrix.size();
    if (n == 0) return {};

    std::vector<bool> visited(n, false);
    std::vector<int> path;

    for (int start = 0; start < n; start++) {
        visited.assign(n, false);
        path.clear();
        path.push_back(start);
        visited[start] = true;

        if (find_hamiltonian_path_helper(adj_matrix, path, visited, start, n)) return path;
    }
    return {};
}

bool is_connected(const AdjacencyList& adj_matrix) {
    int n = adj_matrix.size();
    if (n <= 1) return true;

    std::vector<bool> visited(n, false);
    std::vector<int> stack;

    int start = 0;
    for (int i = 0; i < n; i++) {
        if (!adj_matrix[i].empty()) {
            start = i;
            break;
        }
    }

    stack.push_back(start);
    visited[start] = true;
    int visited_count = 1;
    while (!stack.empty()) {
        int current = stack.back();
        stack.pop_back();

        for (const auto edge : adj_matrix[current]) {
            if (!visited[edge.target]) {
                visited[edge.target] = true;
                stack.push_back(edge.target);
                visited_count++;
            }
        }
    }

    std::set<int> all_nodes;
    for (int i = 0; i < n; i++) {
        if (!adj_matrix[i].empty()) {
            all_nodes.insert(i);
            for (const auto edge : adj_matrix[i]) all_nodes.insert(edge.target);
        }
    }
    return visited_count >= all_nodes.size();
}

std::string sequencing_by_hybridization(const std::vector<std::string>& kmers,
                                        int k,
                                        int target_length) {

    if (kmers.empty()) return "";

    std::cout << "SBH: Building adjacency matrix...\n";
    AdjacencyList adj_matrix = build_adjacency_matrix(kmers, k);
    std::cout << "SBH: Building adjacency matrix finished.\n";

    std::cout << "SBH: Checking connectivity...\n";
    if (!is_connected(adj_matrix)) std::cout << "SBH: Graph is not connected.\n";

    std::cout << "SBH: Attempting to find an Eulerian path...\n";
    std::vector<int> path = find_eulerian_path(adj_matrix);
    if (path.empty()) {
        std::cout << "SBH: No Eulerian path found. Attempting to find a Hamiltonian path...\n";
        path = find_hamiltonian_path(adj_matrix);
    }

    if (path.empty()) {
        std::cout << "SBH: No valid path found in the graph.\n";
        return "";
    }

    std::cout << "SBH: Found path of length: " << path.size() << "\n";
    std::cout << "SBH: Reconstructing sequence...\n";

    return reconstruct_sequence(kmers, k, path, target_length, adj_matrix);
}

std::string sequencing_by_hybridization_with_start(const std::vector<std::string>& kmers,
                                        int k,
                                        int target_length,
                                        const std::string& start_oligo) {

    if (kmers.empty()) return "";

    std::cout << "SBH: Building adjacency matrix...\n";
    AdjacencyList adj_matrix = build_adjacency_matrix(kmers, k);
    std::cout << "SBH: Building adjacency matrix finished.\n";

    int start_idx = 0;
    bool found_start = false;
    for (int i = 0; i < kmers.size(); i++) {
        if (kmers[i] == start_oligo) {
            start_idx = i;
            found_start = true;
            break;
        }
    }

    if (!found_start) {
        std::cout << "SBH: WARN: Start oligo not found in spectrum. Attempting general approach...\n";
        return sequencing_by_hybridization(kmers, k, target_length);
    }

    std::cout << "SBH: Found starting oligo at index: " << start_idx << ".\n";

    std::vector<int> path = find_path_from_start(adj_matrix, start_idx, kmers.size());

    if (path.empty()) {
        std::cout << "SBH: WARN: No path found from starting oligo. Attempting general approach...\n";
        return sequencing_by_hybridization(kmers, k, target_length);
    }

    std::cout << "SBH: Found path of length: " << path.size() << ".\n";
    std::cout << "SBH: Reconstructing sequence...\n";
    return reconstruct_sequence(kmers, k, path, target_length, adj_matrix);
}

std::vector<int> find_path_from_start(const AdjacencyList& adj,
                                        int start,
                                        int total_nodes) {
    std::vector<int> path;
    std::vector<bool> visited(total_nodes, false);

    if (dfs_path(adj, start, visited, path, total_nodes)) {
        return path;
    }

    return {}; // No path found
}

bool dfs_path(const AdjacencyList& adj,
              int current,
              std::vector<bool>& visited,
              std::vector<int>& path,
              int total_nodes) {
    visited[current] = true;
    path.push_back(current);

    // If we've visited all nodes, we found a Hamiltonian path
    if (path.size() == total_nodes) {
        return true;
    }

    // Try all neighbors
    for (const auto edge : adj[current]) {
        if (!visited[edge.target]) {
            if (dfs_path(adj, edge.target, visited, path, total_nodes)) {
                return true;
            }
        }
    }

    // If we can't continue, check if we have a reasonable path
    if (path.size() * 10 >= total_nodes * 8) {  // Accept path covering 80% of nodes
        return true;
    }

    // Backtrack
    path.pop_back();
    visited[current] = false;
    return false;
}


std::vector<std::string> generate_kmers(const std::string& sequence, int k) {
    std::vector<std::string> kmers;
    if (sequence.length() < static_cast<int>(k)) return kmers;

    for (int i = 0; i <= sequence.length() - k; ++i) {
        kmers.push_back(sequence.substr(i, k));
    }
    return kmers;
}

std::pair<int, std::vector<int>> dijkstra_shortest_path(const AdjacencyList& adj_matrix,
                                                              const std::vector<int>& used_count,
                                                              const std::vector<int>& repeat_limits,
                                                              int start,
                                                              const std::vector<int>& targets
                                                              ) {

    int num_kmers = adj_matrix.size();

    constexpr int INF = std::numeric_limits<int>::max();

    std::vector dist(num_kmers, INF);
    std::vector prev(num_kmers, -1);
    dist[start] = 0;

    std::priority_queue<std::pair<int, int>,
                        std::vector<std::pair<int, int>>,
                        std::greater<>> heap;
    heap.push({0, start});

    std::unordered_set target_set(targets.begin(), targets.end());
    int best_target = -1;
    int best_dist = INF;

    while (!heap.empty()) {
        auto [d_u, u] = heap.top();
        heap.pop();

        if (d_u > dist[u]) continue;

        if (target_set.contains(u) && d_u < best_dist) {
            best_dist = d_u;
            best_target = u;
        }
        std::cout << "u = " << u << ", adj_matrix size = " << adj_matrix.size() << "\n";
        if (u >= adj_matrix.size()) {
            std::cerr << "ERROR: Invalid u index = " << u << "\n";
        }
        for (const auto edge : adj_matrix[u]) {
            int v = edge.target;
            if (v >= used_count.size() || v >= repeat_limits.size()) {
                std::cerr << "ERROR: Invalid target index v=" << v
                          << " (used_count.size()=" << used_count.size()
                          << ", repeat_limits.size()=" << repeat_limits.size() << ")\n";
                //exit(EXIT_FAILURE);
            }
            if (used_count[v] >= repeat_limits[v]) continue;

            int new_dist = d_u + edge.weight;
            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                prev[v] = u;
                heap.push({new_dist, v});
            }
        }
    }

    if (best_target == -1) return {-1, {}};

    std::vector<int> path;
    for (int node = best_target; node != -1; node = prev[node]) {
        path.push_back(node);
    }
    std::ranges::reverse(path);

    return {best_target, path};
}
