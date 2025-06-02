import sys
import time
import random
import argparse

def levenshtein(a: str, b: str) -> int:
    n, m = len(a), len(b)
    if n > m:
        a, b = b, a
        n, m = m, n

    current = list(range(n + 1))
    for i in range(1, m + 1):
        previous, current = current, [i] + [0] * n
        for j in range(1, n + 1):
            add = previous[j] + 1
            delete = current[j - 1] + 1
            change = previous[j - 1]
            if a[j - 1] != b[i - 1]:
                change += 1
            current[j] = min(add, delete, change)

    return current[n]

def read_instance(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    if len(lines) < 6:
        raise ValueError("Instance file must have at least 6 header lines and some k-mers.")

    n = int(lines[0])
    k = int(lines[1])
    start_oligo = lines[2]
    num_neg_errors = int(lines[3])
    has_repeats = (lines[4] == 'True' or lines[4] == '1')
    num_pos_errors = int(lines[5])

    kmers = lines[6:]
    return n, k, start_oligo, num_neg_errors, has_repeats, num_pos_errors, kmers

def compute_overlap(u: str, v: str, k: int) -> int:
    for o in range(k - 1, 0, -1):
        if u[-o:] == v[:o]:
            return o
    return 0

def init_pheromones(kmers, k, tau0):
    num_kmers = len(kmers)
    pher = [[0.0] * num_kmers for _ in range(num_kmers)]
    for i in range(num_kmers):
        for j in range(num_kmers):
            if i == j:
                continue
            if compute_overlap(kmers[i], kmers[j], k) > 0:
                pher[i][j] = tau0
    return pher

class Ant:
    def __init__(self, max_steps, num_kmers):
        self.trail = [-1] * max_steps       
        self.trail_len = 0                  
        self.used_count = [0] * num_kmers   
        self.length = 0.0                 
        self.seq_len = 0                    


def reconstruct_sequence(kmers, k, trail, trail_len, n):
    if trail_len == 0:
        return ""
    seq = kmers[trail[0]]
    seq_len = k
    for i in range(1, trail_len):
        u = trail[i - 1]
        v = trail[i]
        ov = compute_overlap(kmers[u], kmers[v], k)
        add_len = k - ov
        if seq_len + add_len > n:
            to_add = n - seq_len
            seq += kmers[v][ov:ov + to_add]
            seq_len = n
            break
        else:
            seq += kmers[v][ov:]
            seq_len += add_len
        if seq_len >= n:
            break
    return seq[:n]

def two_opt_on_trail(ant, kmers, k, n, max_swaps=20):
    best_len = ant.length
    best_trail = ant.trail[:ant.trail_len]
    for _ in range(max_swaps):
        if ant.trail_len < 4:
            break
        i = random.randint(1, ant.trail_len - 3)
        j = random.randint(i + 1, ant.trail_len - 2)
        new_trail = ant.trail[:i] + ant.trail[i:j+1][::-1] + ant.trail[j+1:ant.trail_len]
        new_length = 0.0
        seq_l = k
        valid = True
        for idx in range(1, ant.trail_len):
            u = new_trail[idx - 1]
            v = new_trail[idx]
            ov = compute_overlap(kmers[u], kmers[v], k)
            new_length += (k - ov)
            seq_l += (k - ov)
            if seq_l < n and idx == ant.trail_len - 1:
                valid = False
                break
            if seq_l > n + (k - 1):
                valid = False
                break
        if valid and seq_l >= n and new_length < best_len:
            best_len = new_length
            best_trail = new_trail[:]
    if best_len < ant.length:
        for idx in range(ant.trail_len):
            ant.trail[idx] = best_trail[idx]
        ant.length = best_len
        ant.seq_len = max(ant.seq_len, n)

def update_ants(ants, kmers, k, pher, alpha, beta, num_neg_errors, n):

    num_kmers = len(kmers)
    max_steps = len(ants[0].trail)

    for ant in ants:
        for i in range(num_kmers):
            ant.used_count[i] = 0
        ant.length = 0.0
        ant.trail_len = 1
        ant.seq_len = k

        start_idx = ant.trail[0]
        ant.used_count[start_idx] = 1

        for step in range(1, max_steps):
            u = ant.trail[step - 1]

            used_unique = sum(1 for x in ant.used_count if x > 0)
            remain_unique = num_kmers - used_unique
            max_add = remain_unique * (k - 1)
            to_go = n - ant.seq_len
            if to_go > max_add:
                break

            if ant.seq_len >= n:
                break
            if step < 10:
                candidates = []
                for j in range(num_kmers):
                    if ant.used_count[j] >= 1 + num_neg_errors:
                        continue
                    if compute_overlap(kmers[u], kmers[j], k) == k - 1:
                        candidates.append(j)
                if candidates:
                    chosen = random.choice(candidates)
                    ant.trail[step] = chosen
                    ant.used_count[chosen] += 1
                    ant.trail_len = step + 1
                    ant.length += 1  # weight = k - (k-1) = 1
                    ant.seq_len += 1
                    continue

            probs = [0.0] * num_kmers
            total_prob = 0.0
            for j in range(num_kmers):
                if ant.used_count[j] >= 1 + num_neg_errors:
                    continue
                ov = compute_overlap(kmers[u], kmers[j], k)
                if ov == 0:
                    continue
                tau = pher[u][j]
                if tau <= 0:
                    continue
                fer = tau ** alpha
                vis = (ov ** beta)
                p = fer * vis
                probs[j] = p
                total_prob += p

            if total_prob <= 0.0:
                break

            r = random.random() * total_prob
            cum = 0.0
            chosen = None
            for j in range(num_kmers):
                if probs[j] <= 0:
                    continue
                cum += probs[j]
                if cum >= r:
                    chosen = j
                    break

            if chosen is None:
                break

            ant.trail[step] = chosen
            ant.used_count[chosen] += 1
            ant.trail_len = step + 1
            ov = compute_overlap(kmers[u], kmers[chosen], k)
            w_uv = k - ov
            ant.length += w_uv
            ant.seq_len += (k - ov)

            if ant.seq_len >= n:
                break

        if ant.seq_len >= n:
            two_opt_on_trail(ant, kmers, k, n, max_swaps=20)

def update_pheromones(pher, ants, num_neg_errors, n, Q, rho):

    num_kmers = len(pher)
    for i in range(num_kmers):
        for j in range(num_kmers):
            pher[i][j] *= (1.0 - rho)
            if pher[i][j] < 1e-12:
                pher[i][j] = 0.0

    for ant in ants:
        if ant.seq_len < n or ant.length <= 0:
            continue
        contribution = Q / ant.length
        for pos in range(ant.trail_len - 1):
            u = ant.trail[pos]
            v = ant.trail[pos + 1]
            pher[u][v] += contribution

def run_aco_sbh(instance_file, original_file=None,
                num_ants=100, alpha=3.0, beta=8.0,
                rho=0.5, Q=20.0, tau0=1.0,
                max_time=60, max_iter=0):
    n, k, start_oligo, num_neg_errors, has_repeats, num_pos_errors, kmers = read_instance(instance_file)
    num_kmers = len(kmers)

    try:
        start_idx = kmers.index(start_oligo)
    except ValueError:
        print("Error: start_oligo not found in k-mers.")
        return

    original_seq = None
    if original_file:
        with open(original_file, 'r') as f:
            original_seq = f.read().strip()
            if len(original_seq) < n:
                original_seq = None

    pher = init_pheromones(kmers, k, tau0)

    # Initialize ants
    max_steps = (n - k + 1) + num_neg_errors + 5
    ants = [Ant(max_steps, num_kmers) for _ in range(num_ants)]
    for ant in ants:
        ant.trail[0] = start_idx

    best_sequence = None
    best_length = float('inf')
    best_dist = None

    print(f"Starting ACO: n={n}, k={k}, kmers={num_kmers}, neg_errors={num_neg_errors}")
    print(f"Params: ants={num_ants}, α={alpha}, β={beta}, ρ={rho}, Q={Q}, τ0={tau0}")
    start_time = time.time()
    iter_count = 0

    while True:
        iter_count += 1
        update_ants(ants, kmers, k, pher, alpha, beta, num_neg_errors, n)
        completed = sum(1 for ant in ants if ant.seq_len >= n)
        overlap_hist = [0] * (k)
        for ant in ants:
            for pos in range(ant.trail_len - 1):
                u = ant.trail[pos]
                v = ant.trail[pos + 1]
                ov = compute_overlap(kmers[u], kmers[v], k)
                overlap_hist[ov] += 1

        update_pheromones(pher, ants, num_neg_errors, n, Q, rho)

        improved_this_iter = False
        for ant in ants:
            if ant.seq_len >= n and ant.length < best_length:
                seq = reconstruct_sequence(kmers, k, ant.trail, ant.trail_len, n)
                dist = levenshtein(seq, original_seq) if original_seq else None
                best_length = ant.length
                best_sequence = seq
                best_dist = dist
                improved_this_iter = True

        # Logging
        elapsed = time.time() - start_time
        if improved_this_iter:
            log = (f"[Iter {iter_count:4d} | Time {elapsed:5.1f}s] "
                   f"New best length={best_length:.2f}")
            if best_dist is not None:
                log += f", Levenshtein={best_dist}"
            print(log)
        print(f"[Iter {iter_count}] Completed ants: {completed}/{num_ants} ({completed/num_ants:.1%})")
        hist_str = ", ".join(f"{o}:{overlap_hist[o]}" for o in range(1, k))
        print(f"[Iter {iter_count}] Overlap histogram: {{o:count}} = {hist_str}")

        # Stopping criteria
        if max_iter > 0 and iter_count >= max_iter:
            print(f"Reached max iterations: {max_iter}")
            break
        if max_time > 0 and elapsed >= max_time:
            print(f"Reached max time: {max_time}s")
            break

    print("\n=== Final Result ===")
    if best_sequence is None:
        print(f"No valid sequence of length >= {n} found.")
    else:
        print(f"Best reconstructed sequence (length {n}):")
        print(best_sequence)
        print(f"Suma wag (k-overlap): {best_length:.2f}")
        if best_dist is not None:
            print(f"Levenshtein distance to original: {best_dist}")

def parse_args():
    parser = argparse.ArgumentParser(description="ACO for SBH-DNA reconstruction")
    parser.add_argument("instance", help="Path to instance.txt file")
    parser.add_argument("-o", "--original", help="Path to original.txt for evaluation", default=None)
    parser.add_argument("-a", "--ants", type=int, help="Number of ants", default=100)
    parser.add_argument("--alpha", type=float, help="Alpha (pheromone influence)", default=3.0)
    parser.add_argument("--beta", type=float, help="Beta (visibility/influence overlap)", default=8.0)
    parser.add_argument("--rho", type=float, help="Rho (evaporation rate)", default=0.5)
    parser.add_argument("-q", "--Q", type=float, help="Q (pheromone deposit factor)", default=20.0)
    parser.add_argument("--tau0", type=float, help="Initial pheromone value", default=1.0)
    parser.add_argument("-t", "--time", type=int, help="Max time in seconds (0 = no limit)", default=60)
    parser.add_argument("-i", "--iter", type=int, help="Max iterations (0 = no limit)", default=0)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    random.seed(42)
    run_aco_sbh(
        instance_file=args.instance,
        original_file=args.original,
        num_ants=args.ants,
        alpha=args.alpha,
        beta=args.beta,
        rho=args.rho,
        Q=args.Q,
        tau0=args.tau0,
        max_time=args.time,
        max_iter=args.iter
    )
