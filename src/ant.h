#ifndef ANT_H
#define ANT_H

#include "config.h"
#include "sbh.h"
#include "levenshtein.h"

class Ant {
public:
    std::vector<int> trail;         // indices of k-mers visited
    size_t trail_len = 0;           // how many steps taken
    std::vector<int> used_count;    // usage count for each k-mer
    double length = 0.0;            // total cost
    size_t seq_len = 0;             // reconstructed sequence length

    Ant(size_t max_steps, size_t num_kmers)
        : trail(max_steps, -1),
          used_count(num_kmers, 0) {}
};

struct BestSolution {
    double length = std::numeric_limits<double>::infinity();
    std::vector<int> trail;
    std::string sequence;
    int dist = std::numeric_limits<int>::max();

};

std::pair<size_t, int> choose_random(const std::vector<std::pair<size_t, int>>& candidates);

void two_opt_on_trail(Ant& ant,
                      const std::vector<std::string>& kmers,
                      int k,
                      int n,
                      const AdjacencyList& adj,
                      int max_swaps);

void update_pheromones(std::vector<std::vector<double>>& pheromone,
                       const Ant& ants,
                       int n,
                       double Q,
                       double rho,
                       const Ant& global_best);

void update_ants(const std::vector<Ant>& ants,
                 const std::vector<std::string>& kmers,
                 int k,
                 const AdjacencyList& adj_matrix,
                 std::vector<std::vector<double>>& pheromone,
                 double alpha,
                 double beta,
                 int num_neg_errors,
                 int num_pos_errors,
                 bool has_repeats,
                 const std::vector<int>& repeat_limits,
                 int n);

std::string run_aco(const std::vector<std::string>& kmers,
                                        int k,
                                        size_t n,
                                        const std::string& start_oligo,
                                        int num_neg_errors,
                                        bool has_repeats,
                                        int num_pos_errors,
                                        int num_ants,
                                        double alpha,
                                        double beta,
                                        double rho,
                                        double Q,
                                        double tau0,
                                        int max_time,
                                        int max_iter);

#endif //ANT_H
