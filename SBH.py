#!/usr/bin/env python3
import sys
import time
import random
import argparse
import heapq

def levenshtein(a: str, b: str) -> int:
    """Compute the Levenshtein distance between strings a and b."""
    if b is None:
        return None
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
    if len(start_oligo) != k:
        raise ValueError("start_oligo length must equal k.")

    num_neg_errors = int(lines[3])
    has_repeats = (lines[4] == '1')
    num_pos_errors = int(lines[5])

    kmers = lines[6:]
    return n, k, start_oligo, num_neg_errors, has_repeats, num_pos_errors, kmers

# Build adjacency list with weights using prefix-suffix hashing
def build_adjacency(kmers, k):
    """
    Build a list of neighbors, weights, and overlaps for each k-mer index,
    using hash maps to speed up overlap detection.
    Returns adj where adj[i] = list of (j, weight, overlap) for edges i->j.
    """
    num_kmers = len(kmers)
    # Build prefix_map: for each o in [1..k-1], prefix_map[o][prefix] = list of indices j with kmers[j][:o] = prefix
    prefix_map = {o: {} for o in range(1, k)}
    for j, seq in enumerate(kmers):
        for o in range(1, k):
            pref = seq[:o]
            prefix_map[o].setdefault(pref, []).append(j)

    adj = [[] for _ in range(num_kmers)]
    for i, ui in enumerate(kmers):
        # For each possible overlap length o from k-1 down to 1
        for o in range(k - 1, 0, -1):
            suffix = ui[-o:]
            if suffix in prefix_map[o]:
                for j in prefix_map[o][suffix]:
                    if i == j:
                        continue
                    w = k - o
                    adj[i].append((j, w, o))
    return adj

def init_pheromones(kmers, k, tau0, adj):
    """
    Initializes pheromone matrix for |kmers| x |kmers|.
    pher[i][j] = tau0 if edge i->j exists (i.e., overlap > 0), else 0.
    """
    num_kmers = len(kmers)
    pher = [[0.0] * num_kmers for _ in range(num_kmers)]
    for i in range(num_kmers):
        for (j, _, _) in adj[i]:
            pher[i][j] = tau0
    return pher

class Ant:
    def __init__(self, max_steps, num_kmers):
        self.trail = [-1] * max_steps       # indices of k-mers visited
        self.trail_len = 0                  # how many steps taken
        self.used_count = [0] * num_kmers   # how many times each k-mer was used
        self.length = 0.0                   # sum of weights (k - overlap)
        self.seq_len = 0                    # current reconstructed sequence length


# -----------------------------
# Dijkstra for weighted shortest path (fresh-first)
# -----------------------------
def dijkstra_shortest_path(adj, used_count, repeat_limits, start, targets, k):
    """
    Run Dijkstra from start to each target in 'targets',
    Allowing nodes v only if used_count[v] < repeat_limits[v].
    But priorytetyzując węzły z used_count[v]==0 (świeże) przed tymi, które mają used_count>0.
    Returns (best_target, path_to_best) with minimal sum of weights, or (None, None).
    """
    num_kmers = len(adj)
    INF = float('inf')
    dist = [INF] * num_kmers
    prev = [-1] * num_kmers

    dist[start] = 0
    heap = [(0, start)]

    target_set = set(targets)
    best_target = None
    best_dist = INF

    while heap:
        d_u, u = heapq.heappop(heap)
        if d_u > dist[u]:
            continue
        # Jeśli trafiliśmy na "świeże" węzeł z minimalnym kosztem
        if u in target_set and d_u < best_dist:
            best_dist = d_u
            best_target = u
        for (v, w_uv, _) in adj[u]:
            # Pozwól powrót do v dopóki used_count[v] < repeat_limits[v]
            if used_count[v] >= repeat_limits[v]:
                continue
            alt = d_u + w_uv
            if alt < dist[v]:
                dist[v] = alt
                prev[v] = u
                heapq.heappush(heap, (alt, v))

    if best_target is None:
        return None, None

    # Reconstruct path from start to best_target
    path = []
    node = best_target
    while node != -1:
        path.append(node)
        node = prev[node]
    path.reverse()
    return best_target, path

# Reconstruct sequence from trail
def reconstruct_sequence(kmers, k, trail, trail_len, n):
    """
    Given a trail of k-mer indices, reconstruct the DNA sequence of length n
    by overlapping appropriately.
    """
    if trail_len == 0:
        return ""
    seq = kmers[trail[0]]
    seq_len = k
    for i in range(1, trail_len):
        u = trail[i - 1]
        v = trail[i]
        # compute overlap
        ov = 0
        for o in range(k - 1, 0, -1):
            if kmers[u][-o:] == kmers[v][:o]:
                ov = o
                break
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

# Check if an edge exists from u to v
def edge_weight(adj, kmers, u, v, k):
    """
    Returns the weight (k - overlap) if an edge u->v exists in adj, else None.
    """
    for (jj, w_uv, ov) in adj[u]:
        if jj == v:
            return w_uv
    return None

# Aggressive 2-opt improvement (with validity checks)
def two_opt_on_trail(ant, kmers, k, n, adj, max_swaps=20):
    """
    Perform up to max_swaps random 2-opt attempts on ant's trail to reduce sum of weights,
    but only accept a new trail if every consecutive pair still forms a valid edge
    (overlap >= 1).
    """
    if ant.trail_len < 4:
        return
    best_len = ant.length
    best_trail = ant.trail[:ant.trail_len]
    for _ in range(max_swaps):
        i = random.randint(1, ant.trail_len - 3)
        j = random.randint(i + 1, ant.trail_len - 2)
        new_trail = ant.trail[:i] + ant.trail[i:j+1][::-1] + ant.trail[j+1:ant.trail_len]
        # Check validity of new_trail: every consecutive pair must exist in adj
        valid = True
        new_length = 0.0
        seq_len_acc = k
        for idx in range(1, ant.trail_len):
            u = new_trail[idx - 1]
            v = new_trail[idx]
            w_uv = edge_weight(adj, kmers, u, v, k)
            if w_uv is None:
                valid = False
                break
            new_length += w_uv
            seq_len_acc += w_uv
            if seq_len_acc < n and idx == ant.trail_len - 1:
                valid = False
                break
            if seq_len_acc > n + (k - 1):  # too long
                valid = False
                break
        if valid and seq_len_acc >= n and new_length < best_len:
            best_len = new_length
            best_trail = new_trail[:]
    if best_len < ant.length:
        for idx in range(ant.trail_len):
            ant.trail[idx] = best_trail[idx]
        ant.length = best_len
        ant.seq_len = max(ant.seq_len, n)


# ACO: Build trails for all ants
def update_ants(ants, kmers, k, adj, pher, alpha, beta,
                num_neg_errors, num_pos_errors, has_repeats,
                repeat_limits, n):
    """
    For each ant, build a trail (sequence of k-mer indices) until
    seq_len >= n or max_steps reached, using:
      - TRYB IDEALNY: if num_neg_errors==0 and num_pos_errors==0, pure overlap
      - TRYB PÓŁ-IDEALNY (tylko negatywy lub tylko pozytywy): 
        proste przeszukiwanie overlapy w kolejności o=k-1..1,
        potem Dijkstra, potem virtual jump
      - TRYB PEŁNY (błędy mieszane): ACO → overlapy → Dijkstra → virtual jump
    """
    num_kmers = len(kmers)
    max_steps = len(ants[0].trail)

    ideal = (num_neg_errors == 0 and num_pos_errors == 0)
    only_neg = (num_pos_errors == 0 and num_neg_errors > 0)
    only_pos = (num_neg_errors == 0 and num_pos_errors > 0)

    for ant in ants:
        # Start from the given start_oligo
        start_idx = ant.trail[0]
        # Reset ant state
# w update_ants, na samym początku pętli for ant in ants:
        ant.trail     = [-1] * max_steps
        ant.trail[0]  = start_idx
        ant.trail_len = 1
        ant.seq_len   = k
        ant.length    = 0.0
        for i in range(num_kmers):
            ant.used_count[i] = 0
            ant.used_count[start_idx] = 1

        # Build the trail
        while ant.seq_len < n and ant.trail_len < max_steps:
            u = ant.trail[ant.trail_len - 1]

            # === TRYB IDEALNY ===
            if ideal:
                # Step 1: try overlap = k-1, but w pierwszej kolejności tylko "świeże" węzły
                fresh_candidates = []
                reuse_candidates = []
                for (j, w_uv, ov) in adj[u]:
                    if ov == k - 1:
                        if ant.used_count[j] == 0:
                            fresh_candidates.append((j, w_uv))
                        elif ant.used_count[j] < repeat_limits[j]:
                            reuse_candidates.append((j, w_uv))
                if fresh_candidates:
                    chosen, w_chosen = random.choice(fresh_candidates)
                    ant.trail[ant.trail_len] = chosen
                    ant.used_count[chosen] += 1
                    ant.length += w_chosen
                    ant.seq_len += w_chosen
                    ant.trail_len += 1
                    continue
                if reuse_candidates:
                    chosen, w_chosen = random.choice(reuse_candidates)
                    ant.trail[ant.trail_len] = chosen
                    ant.used_count[chosen] += 1
                    ant.length += w_chosen
                    ant.seq_len += w_chosen
                    ant.trail_len += 1
                    continue

                # Step 2: try overlap = k-1, k-2, ..., 1 (fresh first, then reuse)
                found = False
                for o in range(k - 2, 0, -1):
                    fresh2 = []
                    reuse2 = []
                    for (j, w_uv, ov) in adj[u]:
                        if ov == o:
                            if ant.used_count[j] == 0:
                                fresh2.append((j, w_uv))
                            elif ant.used_count[j] < repeat_limits[j]:
                                reuse2.append((j, w_uv))
                    if fresh2:
                        chosen, w_chosen = random.choice(fresh2)
                        ant.trail[ant.trail_len] = chosen
                        ant.used_count[chosen] += 1
                        ant.length += w_chosen
                        ant.seq_len += w_chosen
                        ant.trail_len += 1
                        found = True
                        break
                    if reuse2:
                        chosen, w_chosen = random.choice(reuse2)
                        ant.trail[ant.trail_len] = chosen
                        ant.used_count[chosen] += 1
                        ant.length += w_chosen
                        ant.seq_len += w_chosen
                        ant.trail_len += 1
                        found = True
                        break
                    # jeśli nic w overlap=o, idź dalej o-1
                if found:
                    continue

                # Step 3: Dijkstra do jakiegokolwiek j, gdzie used_count[j] < repeat_limits[j]
                Y = [j for j in range(num_kmers) if ant.used_count[j] < repeat_limits[j]]
                if Y:
                    best_y, path = dijkstra_shortest_path(adj, ant.used_count, repeat_limits, u, Y, k)
                    if path and len(path) >= 2:
                        for idx in range(1, len(path)):
                            v = path[idx]
                            prev = path[idx - 1]
                            w_uv = edge_weight(adj, kmers, prev, v, k)
                            if w_uv is None:
                                w_uv = k
                            ant.trail[ant.trail_len] = v
                            ant.trail_len += 1
                            ant.used_count[v] += 1
                            ant.length += w_uv
                            ant.seq_len += w_uv
                            if ant.seq_len >= n or ant.trail_len >= max_steps:
                                break
                        continue

                # Step 4: Virtual jump (ov = 0)
                jumped = False
                for j in range(num_kmers):
                    if ant.used_count[j] < repeat_limits[j]:
                        ant.trail[ant.trail_len] = j
                        ant.trail_len += 1
                        ant.used_count[j] += 1
                        ant.length += k
                        ant.seq_len += k
                        jumped = True
                        break
                if jumped:
                    continue

                # Dead end
                break

            # === TRYB PÓŁ-IDEALNY (Tylko negatywy lub tylko pozytywy), ale z wykorzystaniem feromonów ===
            if only_neg or only_pos:
                # Phase 1: ACO-biased selection spośród overlap>=1
                candidates = []
                probs = []
                for (j, w_uv, ov) in adj[u]:
                    if ant.used_count[j] < repeat_limits[j]:
                        pheromone = pher[u][j] if pher[u][j] > 0 else 1e-12
                        heur = ov
                        prob = (pheromone ** alpha) * (heur ** beta)
                        candidates.append((j, w_uv, ov))
                        probs.append(prob)
                if candidates:
                    total = sum(probs)
                    pick = random.uniform(0, total)
                    cum = 0.0
                    chosen, w_chosen, _ = None, None, None
                    for idx_cand, p in enumerate(probs):
                        cum += p
                        if pick <= cum:
                            chosen, w_chosen, _ = candidates[idx_cand]
                            break
                    ant.trail[ant.trail_len] = chosen
                    ant.used_count[chosen] += 1
                    ant.length += w_chosen
                    ant.seq_len += w_chosen
                    ant.trail_len += 1
                    continue

                # Phase 2: try overlap = k-1..1 (fresh first, then reuse)
                found = False
                for o in range(k - 1, 0, -1):
                    fresh = []
                    reuse = []
                    for (j, w_uv, ov) in adj[u]:
                        if ov == o:
                            if ant.used_count[j] == 0:
                                fresh.append((j, w_uv))
                            elif ant.used_count[j] < repeat_limits[j]:
                                reuse.append((j, w_uv))
                    if fresh:
                        chosen, w_chosen = random.choice(fresh)
                        ant.trail[ant.trail_len] = chosen
                        ant.used_count[chosen] += 1
                        ant.length += w_chosen
                        ant.seq_len += w_chosen
                        ant.trail_len += 1
                        found = True
                        break
                    if reuse:
                        chosen, w_chosen = random.choice(reuse)
                        ant.trail[ant.trail_len] = chosen
                        ant.used_count[chosen] += 1
                        ant.length += w_chosen
                        ant.seq_len += w_chosen
                        ant.trail_len += 1
                        found = True
                        break
                if found:
                    continue

                # Phase 3: Dijkstra (fresh first, potem reuse)
                Y_fresh = [j for j in range(num_kmers) if ant.used_count[j] == 0]
                if Y_fresh:
                    best_y, path = dijkstra_shortest_path(adj, ant.used_count, repeat_limits, u, Y_fresh, k)
                else:
                    Y_reuse = [j for j in range(num_kmers) if 0 < ant.used_count[j] < repeat_limits[j]]
                    best_y, path = dijkstra_shortest_path(adj, ant.used_count, repeat_limits, u, Y_reuse, k)

                if path and len(path) >= 2:
                    for idx in range(1, len(path)):
                        v = path[idx]
                        prev = path[idx - 1]
                        w_uv = edge_weight(adj, kmers, prev, v, k)
                        if w_uv is None:
                            w_uv = k
                        ant.trail[ant.trail_len] = v
                        ant.trail_len += 1
                        ant.used_count[v] += 1
                        ant.length += w_uv
                        ant.seq_len += w_uv
                        if ant.seq_len >= n or ant.trail_len >= max_steps:
                            break
                    continue

                # Phase 4: Virtual jump (ov = 0)
                jumped = False
                for j in range(num_kmers):
                    if ant.used_count[j] < repeat_limits[j]:
                        ant.trail[ant.trail_len] = j
                        ant.trail_len += 1
                        ant.used_count[j] += 1
                        ant.length += k
                        ant.seq_len += k
                        jumped = True
                        break
                if jumped:
                    continue

                # Dead end
                break


            # === TRYB PEŁNY (Mieszane błędy) ===
            # Phase 1: ACO-biased selection among overlap >=1
            candidates = []
            probs = []
            for (j, w_uv, ov) in adj[u]:
                if ant.used_count[j] < repeat_limits[j]:
                    pheromone = pher[u][j] if pher[u][j] > 0 else 1e-12
                    heur = ov  # prefer larger overlap
                    prob = (pheromone ** alpha) * (heur ** beta)
                    candidates.append((j, w_uv, ov))
                    probs.append(prob)

            if candidates:
                total = sum(probs)
                pick = random.uniform(0, total)
                cum = 0.0
                chosen, w_chosen, _ = None, None, None
                for idx_cand, p in enumerate(probs):
                    cum += p
                    if pick <= cum:
                        chosen, w_chosen, _ = candidates[idx_cand]
                        break
                ant.trail[ant.trail_len] = chosen
                ant.used_count[chosen] += 1
                ant.length += w_chosen
                ant.seq_len += w_chosen
                ant.trail_len += 1
                continue

            # Phase 2: try overlap = k-1..1 (fresh first, then reuse)
            found = False
            for o in range(k - 1, 0, -1):
                fresh = []
                reuse = []
                for (j, w_uv, ov) in adj[u]:
                    if ov == o:
                        if ant.used_count[j] == 0:
                            fresh.append((j, w_uv))
                        elif ant.used_count[j] < repeat_limits[j]:
                            reuse.append((j, w_uv))
                if fresh:
                    chosen, w_chosen = random.choice(fresh)
                    ant.trail[ant.trail_len] = chosen
                    ant.used_count[chosen] += 1
                    ant.length += w_chosen
                    ant.seq_len += w_chosen
                    ant.trail_len += 1
                    found = True
                    break
                if reuse:
                    chosen, w_chosen = random.choice(reuse)
                    ant.trail[ant.trail_len] = chosen
                    ant.used_count[chosen] += 1
                    ant.length += w_chosen
                    ant.seq_len += w_chosen
                    ant.trail_len += 1
                    found = True
                    break
                # jeśli nic overlap=o, przeglądaj dalej o-1
            if found:
                continue

            # Phase 3: Dijkstra (fresh first, then reuse)
            Y_fresh = [j for j in range(num_kmers) if ant.used_count[j] == 0]
            if Y_fresh:
                best_y, path = dijkstra_shortest_path(adj, ant.used_count, repeat_limits, u, Y_fresh, k)
            else:
                Y_reuse = [j for j in range(num_kmers) if 0 < ant.used_count[j] < repeat_limits[j]]
                best_y, path = dijkstra_shortest_path(adj, ant.used_count, repeat_limits, u, Y_reuse, k)

            if path and len(path) >= 2:
                for idx in range(1, len(path)):
                    v = path[idx]
                    prev = path[idx - 1]
                    w_uv = edge_weight(adj, kmers, prev, v, k)
                    if w_uv is None:
                        w_uv = k
                    ant.trail[ant.trail_len] = v
                    ant.trail_len += 1
                    ant.used_count[v] += 1
                    ant.length += w_uv
                    ant.seq_len += w_uv
                    if ant.seq_len >= n or ant.trail_len >= max_steps:
                        break
                continue

            # Phase 4: Virtual jump (ov = 0)
            jumped = False
            for j in range(num_kmers):
                if ant.used_count[j] < repeat_limits[j]:
                    ant.trail[ant.trail_len] = j
                    ant.trail_len += 1
                    ant.used_count[j] += 1
                    ant.length += k
                    ant.seq_len += k
                    jumped = True
                    break
            if jumped:
                continue

            # Dead end
            break

        # After building a path, attempt local 2-opt improvement if valid
        if ant.seq_len >= n:
            two_opt_on_trail(ant, kmers, k, n, adj, max_swaps=20)

# ACO: Update pheromones
def update_pheromones(pher, ants, n, Q, rho, global_best):
    num_kmers = len(pher)
    # Evaporation
    for i in range(num_kmers):
        for j in range(num_kmers):
            old = pher[i][j]
            pher[i][j] *= (1.0 - rho)
           # if abs(pher[i][j] - old) > 1e-8:
                # Wypisz tylko kilkanaście “zmian”, by nie spamować konsoli
            #   print(f"DEBUG: Pho[{i}][{j}] {old:.4e} → {pher[i][j]:.4e}")

    # Deposit from each qualifying ant
    for ant in ants:
        if ant.seq_len < n or ant.length <= 0:
            continue
        contribution = Q / ant.length
        for pos in range(ant.trail_len - 1):
            u = ant.trail[pos]
            v = ant.trail[pos + 1]
            pher[u][v] += contribution

    # Elitist deposit (global best)
    if global_best.get('trail') is not None and len(global_best['trail']) > 1:
        contribution = Q / global_best['length']
        trail = global_best['trail']
        for idx in range(len(trail) - 1):
            u = trail[idx]
            v = trail[idx + 1]
            pher[u][v] += contribution

# Main ACO routine for SBH-DNA
def run_aco_sbh(instance_file, original_file=None,
                num_ants=100, alpha=1.0, beta=2.0,
                rho=0.1, Q=20.0, tau0=1.0,
                max_time=60, max_iter=0):
    # Read instance
    n, k, start_oligo, num_neg_errors, has_repeats, num_pos_errors, kmers = read_instance(instance_file)
    num_kmers = len(kmers)

    print(f"DEBUG: n={n}, k={k}, start_oligo='{start_oligo}', "
          f"num_neg_errors={num_neg_errors}, has_repeats={has_repeats}, num_pos_errors={num_pos_errors}")
    print(f"DEBUG: pierwsze 5 k-merów: {kmers[:5]}")

    # Find index of start_oligo
    try:
        start_idx = kmers.index(start_oligo)
    except ValueError:
        print("Error: start_oligo not found in k-mers.")
        return

    # Build adjacency
    adj = build_adjacency(kmers, k)
    neighbors_kminus1 = [j for (j, _, ov) in adj[start_idx] if ov == k - 1]
    print(f"DEBUG: z start_oligo wychodzi {len(neighbors_kminus1)} łuków overlap=k-1")

    # Read original sequence if provided
    original_seq = None
    if original_file:
        with open(original_file, 'r') as f:
            original_seq = f.read().strip()
            if len(original_seq) < n:
                original_seq = None

    pher = init_pheromones(kmers, k, tau0, adj)
    # Jeśli powtórzenia, limit = num_neg_errors+1; w przeciwnym razie 1
    repeat_limit_value = (num_neg_errors + 1) if has_repeats else 1
    repeat_limits = [repeat_limit_value] * num_kmers

    max_steps = len(kmers)
    ants = [Ant(max_steps, num_kmers) for _ in range(num_ants)]
    for ant in ants:
        ant.trail[0] = start_idx

    global_best = {
        'length': float('inf'),
        'trail': None,
        'dist': None,
        'seq': None
    }

    print(f"Starting ACO: n={n}, k={k}, kmers={num_kmers}, neg_errors={num_neg_errors}, pos_errors={num_pos_errors}")
    print(f"Params: ants={num_ants}, α={alpha}, β={beta}, ρ={rho}, Q={Q}, τ0={tau0}")
    start_time = time.time()
    iter_count = 0

    while True:
        iter_count += 1

        # Build trails
        update_ants(ants, kmers, k, adj, pher, alpha, beta,
                    num_neg_errors, num_pos_errors, has_repeats,
                    repeat_limits, n)

        # Count completed ants and overlap histogram
        completed = sum(1 for ant in ants if ant.seq_len >= n)
        overlap_hist = [0] * k
        for ant in ants:
            for pos in range(ant.trail_len - 1):
                u = ant.trail[pos]
                v = ant.trail[pos + 1]
                ov = 0
                for o in range(k - 1, 0, -1):
                    if kmers[u][-o:] == kmers[v][:o]:
                        ov = o
                        break
                overlap_hist[ov] += 1

        # Check each ant for a valid solution and update best
        improved_this_iter = False
        for ant in ants:
            if ant.seq_len < n:
                continue

            # zawsze przeliczamy odległość Levenshteina
            seq = reconstruct_sequence(kmers, k, ant.trail, ant.trail_len, n)
            dist = levenshtein(seq, original_seq)

            # Jeśli wcześniej nie było żadnego best
            if global_best['trail'] is None:
                global_best['length'] = ant.length
                global_best['dist'] = dist
                global_best['trail'] = ant.trail[:ant.trail_len]
                global_best['seq'] = seq
                improved_this_iter = True

            else:
                # strategia (B): najpierw porównujemy length, gdy równe → Levenshtein
                if ant.length < global_best['length'] or \
                   (ant.length == global_best['length'] and dist < global_best['dist']):
                    global_best['length'] = ant.length
                    global_best['dist'] = dist
                    global_best['trail'] = ant.trail[:ant.trail_len]
                    global_best['seq'] = seq
                    improved_this_iter = True

                # alternatywnie: jeśli priorytet ma Levenshtein → użyj:
                #if dist < global_best['dist']:
                #     global_best['length'] = ant.length
                #     global_best['dist'] = dist
                 #    global_best['trail'] = ant.trail[:ant.trail_len]
                 #    global_best['seq'] = seq
                 #    improved_this_iter = True


        # Update pheromones with elitist deposit, ale tylko jeśli chociaż jedna mrówka ukończyła
        if completed > 0:
            update_pheromones(pher, ants, n, Q, rho, global_best)

        # Logging
        elapsed = time.time() - start_time
        if improved_this_iter:
            log = (f"[Iter {iter_count:4d} | Time {elapsed:5.1f}s] "
                   f"New best length={global_best['length']:.2f}")
            if global_best['dist'] is not None:
                log += f", Levenshtein={global_best['dist']}"
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

    # Output results
    print("\n=== Final Result ===")
    if global_best['trail'] is None:
        print(f"No valid sequence of length >= {n} found.")
    else:
        print(f"Best reconstructed sequence (length {n}):")
        print(global_best['seq'])
        print(f"Suma wag (k - overlap): {global_best['length']:.2f}")
        if global_best['dist'] is not None:
            print(f"Levenshtein distance to original: {global_best['dist']}")

# -----------------------------
# Command-line interface
# -----------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="ACO for SBH-DNA reconstruction")
    parser.add_argument("instance", help="Path to instance.txt file")
    parser.add_argument("-o", "--original", help="Path to original.txt for evaluation", default=None)
    parser.add_argument("-a", "--ants", type=int, help="Number of ants", default=100)
    parser.add_argument("--alpha", type=float, help="Alpha (pheromone influence)", default=1.0)
    parser.add_argument("--beta", type=float, help="Beta (heuristic influence)", default=2.0)
    parser.add_argument("--rho", type=float, help="Rho (evaporation rate)", default=0.1)
    parser.add_argument("-q", "--Q", type=float, help="Q (pheromone deposit factor)", default=20.0)
    parser.add_argument("--tau0", type=float, help="Initial pheromone value", default=1.0)
    parser.add_argument("-t", "--time", type=int, help="Max time in seconds (0 = no limit)", default=60)
    parser.add_argument("-i", "--iter", type=int, help="Max iterations (0 = no limit)", default=0)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    random.seed()
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
