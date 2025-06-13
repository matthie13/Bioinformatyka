#include "aco.h"
#include "levenshtein.h"
#include "settings.h"

std::pair<int, int> choose_random(const std::vector<std::pair<int, int>>& candidates) {
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, candidates.size() - 1);
    return candidates[dist(rng)];
}

void two_opt_on_trail(Ant& ant,
                      const std::vector<std::string>& kmers,
                      int k,
                      int n,
                      const AdjacencyList& adj,
                      int max_swaps) {

    if (ant.trail_len < 4) return;

    double best_len = ant.cost;
    std::vector<int> best_trail(ant.trail.begin(), ant.trail.begin() + ant.trail_len);

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist_i(1, ant.trail_len - 3);

    for (int attempt = 0; attempt < max_swaps; ++attempt) {
        int i = dist_i(rng);
        std::uniform_int_distribution<int> dist_j(i + 1, ant.trail_len - 2);
        int j = dist_j(rng);

        // Build new trail with reversed segment
        std::vector<int> new_trail;
        new_trail.insert(new_trail.end(), ant.trail.begin(), ant.trail.begin() + i);
        for (int x = j; x >= i; --x) new_trail.push_back(ant.trail[x]);
        new_trail.insert(new_trail.end(), ant.trail.begin() + j + 1, ant.trail.begin() + ant.trail_len);

        // Validate new trail
        bool valid = true;
        double new_length = 0.0;
        int seq_len_acc = k;

        for (int idx = 1; idx < ant.trail_len; ++idx) {
            int u = new_trail[idx - 1];
            int v = new_trail[idx];

            std::optional<int> w_uv = edge_weight(adj, kmers, u, v, k);
            if (!w_uv.has_value()) {
                valid = false;
                break;
            }
            new_length += *w_uv;
            seq_len_acc += *w_uv;

            if (seq_len_acc > n + (k - 1)) {
                valid = false;
                break;
            }

            if (idx == ant.trail_len - 1 && seq_len_acc < n) {
                valid = false;
                break;
            }
        }

        if (valid && new_length < best_len) {
            best_len = new_length;
            best_trail = new_trail;
        }
    }

    // Apply improved trail
    if (best_len < ant.cost) {
        for (int i = 0; i < ant.trail_len; ++i)
            ant.trail[i] = best_trail[i];
        ant.cost = best_len;
        if (n >= 0) {
            ant.seq_len = std::max(ant.seq_len, static_cast<int>(n));
        }
    }
}

void update_pheromones(std::vector<std::vector<double>>& pheromone,
                       const std::vector<Ant>& ants,
                       int n,
                       double Q,
                       double rho,
                       const struct BestSolution& global_best) {

    int num_kmers = pheromone.size();
    double evaporation_factor = 1.0 - rho;

    for (int i = 0; i < num_kmers; ++i) {
        for (int j = 0; j < pheromone[i].size(); ++j) {
            pheromone[i][j] *= evaporation_factor;
        }
    }

    for (const Ant& ant : ants) {
        if (ant.seq_len < n || ant.cost <= 0) continue;
        double contribution = Q / ant.cost;
        for (int pos = ant.trail_len - 1; pos > 0; --pos) {
            int u = ant.trail[pos];
            int v = ant.trail[pos - 1];
            pheromone[u][v] += contribution;
        }
    }
    if (!global_best.trail.empty() && global_best.trail.size() > 1) {
        double contribution = Q / global_best.length;
        auto trail = global_best.trail;
        for (int pos = trail.size() - 1; pos > 0; --pos) {
            int u = trail[pos];
            int v = trail[pos - 1];
            pheromone[u][v] += contribution;
        }
    }
}

void update_ants(std::vector<Ant>& ants,
                 const std::vector<std::string>& kmers,
                 int k,
                 const AdjacencyList& adj_matrix,
                 const std::vector<std::vector<double>>& pheromone,
                 double alpha,
                 double beta,
                 int num_neg_errors,
                 int num_pos_errors,
                 bool has_repeats,
                 const std::vector<int>& repeat_limits,
                 int n) {

    const int num_kmers = kmers.size();
    const int max_steps = ants[0].trail.size();
    const bool ideal = (num_neg_errors == 0 && num_pos_errors == 0);
    const bool only_neg = (num_neg_errors > 0 && num_pos_errors == 0);
    const bool only_pos = (num_neg_errors == 0 && num_pos_errors > 0);
    const bool both = (num_neg_errors > 0 && num_pos_errors > 0);

    for (Ant& ant : ants) {
        auto start_idx = ant.trail[0];
        std::vector<int> fresh_trail(max_steps, -1);
        ant.trail = fresh_trail;
        ant.trail[0] = start_idx;
        ant.trail_len = 1;
        ant.seq_len = k;
        ant.cost = 0.0;
        for (int i = 0; i < num_kmers; i++) {
            ant.used_count[i] = 0;
            ant.used_count[start_idx] = 1;
        }

        while (ant.seq_len < n && ant.trail_len < max_steps) {
            auto u = ant.trail[ant.trail_len - 1];

            if (ideal) {
                // Krok 1: Overlap k-1
                std::vector<std::pair<int, int>> fresh_candidates;
                std::vector<std::pair<int, int>> reuse_candidates;
                for (const auto edge : adj_matrix[u]) {
                    if (edge.overlap == k - 1) {
                        if (ant.used_count[edge.target] == 0)
                            fresh_candidates.emplace_back(edge.target, edge.weight);
                        else if (ant.used_count[edge.target] < repeat_limits[edge.target])
                            reuse_candidates.emplace_back(edge.target, edge.weight);
                    }
                }

                if (!fresh_candidates.empty()) {
                    auto [chosen, w_chosen] = choose_random(fresh_candidates);
                    ant.trail[ant.trail_len] = chosen;
                    ant.used_count[chosen]++;
                    ant.cost += w_chosen;
                    ant.seq_len += w_chosen;
                    ant.trail_len++;
                    continue;
                }

                if (!reuse_candidates.empty()) {
                    auto [chosen, w_chosen] = choose_random(reuse_candidates);
                    ant.trail[ant.trail_len] = chosen;
                    ant.used_count[chosen]++;
                    ant.cost += w_chosen;
                    ant.seq_len += w_chosen;
                    ant.trail_len++;
                    continue;
                }

                // Krok 2: Overlap k-1, k-2, ..., 1
                bool found = false;
                for (int o = k - 2; o >= 1; --o) {
                    std::vector<std::pair<int, int>> new_fresh_candidates;
                    std::vector<std::pair<int, int>> new_reuse_candidates;

                    for (const auto edge : adj_matrix[u]) {
                        if (edge.overlap == o) {
                            if (ant.used_count[edge.target] == 0)
                                new_fresh_candidates.emplace_back(edge.target, edge.weight);
                            else if (ant.used_count[edge.target] < repeat_limits[edge.target])
                                new_reuse_candidates.emplace_back(edge.target, edge.weight);
                        }
                    }

                    if (!new_fresh_candidates.empty()) {
                        auto [chosen, w_chosen] = choose_random(new_fresh_candidates);
                        ant.trail[ant.trail_len] = chosen;
                        ant.used_count[chosen]++;
                        ant.cost += w_chosen;
                        ant.seq_len += w_chosen;
                        ant.trail_len++;
                        found = true;
                        break;
                    }

                    if (!new_reuse_candidates.empty()) {
                        auto [chosen, w_chosen] = choose_random(new_reuse_candidates);
                        ant.trail[ant.trail_len] = chosen;
                        ant.used_count[chosen]++;
                        ant.cost += w_chosen;
                        ant.seq_len += w_chosen;
                        ant.trail_len++;
                        found = true;
                        break;
                    }
                }
                if (found) continue;

                // Krok 3: Dijkstra
                std::vector<int> Y;
                for (int j = 0; j < num_kmers; ++j) {
                    if (ant.used_count[j] < repeat_limits[j]) {
                        Y.push_back(j);
                    }
                }

                if (!Y.empty()) {
                    auto [best_y, path] = dijkstra_shortest_path(adj_matrix, ant.used_count, repeat_limits, u, Y);

                    if (!path.empty() && path.size() >= 2) {
                        for (int idx = 1; idx < path.size(); ++idx) {
                            int v = path[idx];
                            int prev = path[idx - 1];

                            std::optional<int> weight_opt = edge_weight(adj_matrix, kmers, prev, v, k);
                            int w_uv = weight_opt.has_value() ? *weight_opt : k;

                            ant.trail[ant.trail_len++] = v;
                            ant.used_count[v]++;
                            ant.cost += w_uv;
                            ant.seq_len += w_uv;

                            if (ant.seq_len >= n || ant.trail_len >= max_steps)
                                break;
                        }
                        continue;
                    }
                }
                // Krok 4: virtual jump (overlap = 0)
                bool jumped = false;
                for (auto j = 0; j < num_kmers; ++j) {
                    if (ant.used_count[j] < repeat_limits[j]) {
                        ant.trail[ant.trail_len++] = j;
                        ant.used_count[j]++;
                        ant.cost += k;
                        ant.seq_len += k;
                        jumped = true;
                        break;
                    }
                }
                if (jumped) continue;
                break;
            }
            if (only_neg || only_pos) {
                // Krok 1: ACO
                std::vector<std::tuple<int, int, int>> candidates;
                std::vector<double> probs;
                for (const auto edge : adj_matrix[u]) {
                    int j = edge.target;
                    if (ant.used_count[j] < repeat_limits[j]) {
                        double pher = pheromone[u][j] > 0.0 ? pheromone[u][j] : 1e-12;
                        auto heur = static_cast<double>(edge.overlap);
                        double prob = std::pow(pher, alpha) * std::pow(heur, beta);
                        candidates.emplace_back(j, edge.weight, edge.overlap);
                        probs.emplace_back(prob);
                    }
                }

                if (!candidates.empty()) {
                    double total = std::accumulate(probs.begin(), probs.end(), 0.0);
                    std::uniform_real_distribution<double> dist(0.0, total);
                    static std::mt19937 rng(std::random_device{}());
                    double pick = dist(rng);
                    double cum = 0.0;

                    for (int i = 0; i < probs.size(); ++i) {
                        cum += probs[i];
                        if (pick <= cum) {
                            auto [chosen, w_chosen, _] = candidates[i];
                            ant.trail[ant.trail_len] = chosen;
                            ant.used_count[chosen]++;
                            ant.cost += w_chosen;
                            ant.seq_len += w_chosen;
                            ant.trail_len++;
                            break;
                        }
                    }
                }
                // Krok 2: Overlap k-1 ... 1
                bool found = false;
                for (int o = k - 1; o >= 1; --o) {
                    std::vector<std::pair<int, int>> fresh_candidates;
                    std::vector<std::pair<int, int>> reuse_candidates;

                    for (const auto edge : adj_matrix[u]) {
                        if (edge.overlap == o) {
                            if (ant.used_count[edge.target] == 0)
                                fresh_candidates.emplace_back(edge.target, edge.weight);
                            else if (ant.used_count[edge.target] < repeat_limits[edge.target])
                                reuse_candidates.emplace_back(edge.target, edge.weight);
                        }
                    }

                    if (!fresh_candidates.empty()) {
                        auto [chosen, w_chosen] = choose_random(fresh_candidates);
                        ant.trail[ant.trail_len] = chosen;
                        ant.used_count[chosen]++;
                        ant.cost += w_chosen;
                        ant.seq_len += w_chosen;
                        ant.trail_len++;
                        found = true;
                        break;
                    }

                    if (!reuse_candidates.empty()) {
                        auto [chosen, w_chosen] = choose_random(reuse_candidates);
                        ant.trail[ant.trail_len] = chosen;
                        ant.used_count[chosen]++;
                        ant.cost += w_chosen;
                        ant.seq_len += w_chosen;
                        ant.trail_len++;
                        found = true;
                        break;
                    }
                }
                if (found) continue;

                // Krok 3: Dijkstra
                std::vector<int> Y;
                for (int j = 0; j < num_kmers; ++j) {
                    if (ant.used_count[j] < repeat_limits[j]) {
                        Y.push_back(j);
                    }
                }

                if (!Y.empty()) {
                    auto [best_y, path] = dijkstra_shortest_path(adj_matrix, ant.used_count, repeat_limits, u, Y);
                    if (!path.empty() && path.size() >= 2) {
                        for (int idx = 1; idx < path.size(); ++idx) {
                            int v = path[idx];
                            int prev = path[idx - 1];

                            std::optional<int> weight_opt = edge_weight(adj_matrix, kmers, prev, v, k);
                            int w_uv = weight_opt.has_value() ? *weight_opt : k;

                            ant.trail[ant.trail_len++] = v;
                            ant.used_count[v]++;
                            ant.cost += w_uv;
                            ant.seq_len += w_uv;

                            if (ant.seq_len >= n || ant.trail_len >= max_steps)
                                break;
                        }
                        continue;
                    }
                }

                // Krok 4: virtual jump (overlap = 0)
                bool jumped = false;
                for (auto j = 0; j < num_kmers; ++j) {
                    if (ant.used_count[j] < repeat_limits[j]) {
                        ant.trail[ant.trail_len++] = j;
                        ant.used_count[j]++;
                        ant.cost += k;
                        ant.seq_len += k;
                        jumped = true;
                        break;
                    }
                }
                if (jumped) continue;
                break;
            }
            if (both) {
                // Krok 1: ACO
                std::vector<std::tuple<int, int, int>> candidates;
                std::vector<double> probs;
                for (const auto edge : adj_matrix[u]) {
                    int j = edge.target;
                    if (ant.used_count[j] < repeat_limits[j]) {
                        double pher = pheromone[u][j] > 0.0 ? pheromone[u][j] : 1e-12;
                        auto heur = static_cast<double>(edge.overlap);
                        double prob = std::pow(pher, alpha) * std::pow(heur, beta);
                        candidates.emplace_back(j, edge.weight, edge.overlap);
                        probs.emplace_back(prob);
                    }
                }
                if (!candidates.empty()) {
                    double total = std::accumulate(probs.begin(), probs.end(), 0.0);
                    std::uniform_real_distribution<double> dist(0.0, total);
                    static std::mt19937 rng(std::random_device{}());
                    double pick = dist(rng);
                    double cum = 0.0;

                    for (int i = 0; i < probs.size(); ++i) {
                        cum += probs[i];
                        if (pick <= cum) {
                            auto [chosen, w_chosen, _] = candidates[i];
                            ant.trail[ant.trail_len] = chosen;
                            ant.used_count[chosen]++;
                            ant.cost += w_chosen;
                            ant.seq_len += w_chosen;
                            ant.trail_len++;
                            break;
                        }
                    }
                }
                // Krok 2: Overlap k-1 ... 1
                bool found = false;
                for (int o = k - 1; o >= 1; --o) {
                    std::vector<std::pair<int, int>> fresh_candidates;
                    std::vector<std::pair<int, int>> reuse_candidates;

                    for (const auto edge : adj_matrix[u]) {
                        if (edge.overlap == o) {
                            if (ant.used_count[edge.target] == 0)
                                fresh_candidates.emplace_back(edge.target, edge.weight);
                            else if (ant.used_count[edge.target] < repeat_limits[edge.target])
                                reuse_candidates.emplace_back(edge.target, edge.weight);
                        }
                    }

                    if (!fresh_candidates.empty()) {
                        auto [chosen, w_chosen] = choose_random(fresh_candidates);
                        ant.trail[ant.trail_len] = chosen;
                        ant.used_count[chosen]++;
                        ant.cost += w_chosen;
                        ant.seq_len += w_chosen;
                        ant.trail_len++;
                        found = true;
                        break;
                    }

                    if (!reuse_candidates.empty()) {
                        auto [chosen, w_chosen] = choose_random(reuse_candidates);
                        ant.trail[ant.trail_len] = chosen;
                        ant.used_count[chosen]++;
                        ant.cost += w_chosen;
                        ant.seq_len += w_chosen;
                        ant.trail_len++;
                        found = true;
                        break;
                    }
                }
                if (found) continue;

                // Krok 3: Dijkstra
                std::vector<int> Y;
                for (int j = 0; j < num_kmers; ++j) {
                    if (ant.used_count[j] < repeat_limits[j]) {
                        Y.push_back(j);
                    }
                }

                if (!Y.empty()) {
                    auto [best_y, path] = dijkstra_shortest_path(adj_matrix, ant.used_count, repeat_limits, u, Y);
                    if (!path.empty() && path.size() >= 2) {
                        for (int idx = 1; idx < path.size(); ++idx) {
                            int v = path[idx];
                            int prev = path[idx - 1];

                            std::optional<int> weight_opt = edge_weight(adj_matrix, kmers, prev, v, k);
                            int w_uv = weight_opt.has_value() ? *weight_opt : k;

                            ant.trail[ant.trail_len++] = v;
                            ant.used_count[v]++;
                            ant.cost += w_uv;
                            ant.seq_len += w_uv;

                            if (ant.seq_len >= n || ant.trail_len >= max_steps)
                                break;
                        }
                        continue;
                    }
                }

                // Krok 4: virtual jump (overlap = 0)
                bool jumped = false;
                for (auto j = 0; j < num_kmers; ++j) {
                    if (ant.used_count[j] < repeat_limits[j]) {
                        ant.trail[ant.trail_len++] = j;
                        ant.used_count[j]++;
                        ant.cost += k;
                        ant.seq_len += k;
                        jumped = true;
                        break;
                    }
                }
                if (jumped) continue;
                break;
            }
        }
        if (ant.seq_len >= n) two_opt_on_trail(ant, kmers, k, n, adj_matrix, 20);
    }
}

std::string run_aco(const std::vector<std::string>& kmers,
                                        int k,
                                        int n,
                                        const std::string& start_oligo,
                                        int num_neg_errors,
                                        bool has_repeats,
                                        int num_pos_errors,
                                        int num_ants = settings::NUM_ANTS,
                                        double alpha = settings::ALPHA,
                                        double beta = settings::BETA,
                                        double rho = settings::RHO,
                                        double Q = settings::Q,
                                        double tau0 = settings::TAU0,
                                        int max_time = settings::MAX_TIME,
                                        int max_iter = settings::MAX_ITER) {
    int start_idx = 0;
    bool found_start = false;
    for (int i = 0; i < kmers.size(); i++) {
        if (kmers[i] == start_oligo) {
            start_idx = i;
            found_start = true;
            break;
        }
    }
    if (!found_start) std::cout << "ACO: WARN: Start oligo not found in spectrum. Attempting general approach...\n";

    AdjacencyList adj_matrix = build_adjacency_matrix(kmers, k);
    std::vector<std::vector<double>> pheromone(kmers.size(), std::vector<double>(kmers.size(), 0.0));
    for (int u = 0; u < kmers.size(); ++u)
        for (const auto [v, w, ov] : adj_matrix[u])
            pheromone[u][v] = tau0;

    int repeat_limit_value = has_repeats ? num_neg_errors + 1 : 1;
    std::vector<int> repeat_limits(kmers.size(), repeat_limit_value);

    std::vector<Ant> ants;
    int max_steps = kmers.size();
    for (int i = 0; i < num_ants; ++i) {
        Ant ant(max_steps, kmers.size());
        if (found_start) {
            ant.trail[0] = start_idx;
        } else {
            // Pick a random starting k-mer if the given one isn't found
            int random_start = std::rand() % kmers.size();
            ant.trail[0] = random_start;
        }
        ant.used_count[ant.trail[0]]++;
        ant.cost++;
        ant.trail_len = 1;
        ants.emplace_back(std::move(ant));
    }

    struct BestSolution global_best;

    auto start = std::chrono::steady_clock::now();
    int iter = 0;
    while (true) {
        iter++;

        update_ants(ants, kmers, k, adj_matrix, pheromone, alpha, beta,
                    num_neg_errors, num_pos_errors, has_repeats,
                    repeat_limits, n);

        for (const Ant& ant : ants) {
            if (ant.seq_len < n) continue;

            std::vector<int> trail_indices(ant.trail.begin(), ant.trail.end());
            std::string seq = reconstruct_sequence(kmers, k, trail_indices, n, adj_matrix);
            double length = ant.cost;
            int dist = levenshtein_score(seq, "");

            if (length < global_best.length || (length == global_best.length && dist < global_best.dist)) {
                global_best.length = length;
                global_best.trail = std::vector<int>(ant.trail.begin(), ant.trail.begin() + ant.trail_len);
                global_best.sequence = std::move(seq);
                global_best.dist = dist;
            }
        }
        update_pheromones(pheromone, ants, n, Q, rho, global_best);

        auto now = std::chrono::steady_clock::now();
        if ((max_time > 0 && std::chrono::duration_cast<std::chrono::seconds>(now - start).count() >= max_time) ||
            (max_iter > 0 && iter >= max_iter)) break;
    }
    return global_best.sequence;
}
