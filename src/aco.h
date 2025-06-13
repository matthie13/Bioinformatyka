#ifndef ANT_H
#define ANT_H

#include "config.h"
#include "sbh.h"


struct Ant {
    std::vector<int> trail;         // indices of k-mers visited
    int trail_len;           // how many steps taken
    std::vector<int> used_count;    // usage count for each k-mer
    double cost;            // total cost
    int seq_len;             // reconstructed sequence length

    Ant(int max_steps, int num_kmers)
        : trail(max_steps, -1), trail_len(0),
          used_count(num_kmers, 0), cost(0.0), seq_len(0) {
    }
};

struct BestSolution {
    double length = std::numeric_limits<double>::infinity();
    std::vector<int> trail;
    std::string sequence;
    int dist = std::numeric_limits<int>::max();
};

std::pair<int, int> choose_random(const std::vector<std::pair<int, int>>& candidates);

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
                 int n);

std::string run_aco(const std::vector<std::string>& kmers,
                                        int k,
                                        int n,
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
