import random
from collections import Counter

def generate_single_instance(n, k, num_neg_errors, num_pos_errors, has_repeats):
    if has_repeats:
        n0 = int(0.6 * n)
        base_seq = ''.join(random.choices('ATCG', k=n0))
        repeats = []
        num_repeat_kmers = max(1, n // 100)
        for _ in range(num_repeat_kmers):
            i = random.randrange(0, n0 - k + 1)
            segment = base_seq[i:i+k]
            repeats.append(segment)
        dna_list = list(base_seq)
        for segment in repeats:
            if len(dna_list) < k:
                break
            insert_pos = random.randrange(0, len(dna_list) - k + 1)
            dna_list[insert_pos:insert_pos+k] = list(segment)
        while len(dna_list) < n:
            dna_list.append(random.choice('ATCG'))
        dna = ''.join(dna_list[:n])
    else:
        dna = ''.join(random.choices('ATCG', k=n))

    all_kmers = [dna[i : i+k] for i in range(n - k + 1)]
    ideal_set = set(all_kmers)
    spectrum = list(ideal_set)
    to_remove = []

    if num_neg_errors > 0:
        if has_repeats:
            counts = Counter(all_kmers)
            repeat_kmers = [kmer for kmer, c in counts.items() if c >= 2]
            random.shuffle(repeat_kmers)
            r = min(len(repeat_kmers), num_neg_errors)
            to_remove = repeat_kmers[:r]
            if r < num_neg_errors:
                others = list(ideal_set - set(repeat_kmers))
                random.shuffle(others)
                to_remove += others[: (num_neg_errors - r) ]
        else:
            removable = spectrum.copy()
            random.shuffle(removable)
            to_remove = removable[: min(len(removable), num_neg_errors) ]
        for kmer in to_remove:
            if kmer in spectrum:
                spectrum.remove(kmer)

    existing = set(spectrum)
    ideal_list_for_pos = list(ideal_set)
    pos_errors_added = 0
    idx = 0
    while pos_errors_added < num_pos_errors and idx < len(ideal_list_for_pos):
        base = ideal_list_for_pos[idx]
        idx += 1
        last = base[-1]
        for alt in 'ACGT':
            if alt != last:
                fake1 = base[:-1] + alt
                if (fake1 not in ideal_set) and (fake1 not in existing):
                    spectrum.append(fake1)
                    existing.add(fake1)
                    pos_errors_added += 1
                    break
        if pos_errors_added >= num_pos_errors:
            break
        mid_index = (k // 2)
        mid = base[mid_index]
        for alt in 'ACGT':
            if alt != mid:
                fake2 = base[:mid_index] + alt + base[mid_index+1:]
                if (fake2 not in ideal_set) and (fake2 not in existing):
                    spectrum.append(fake2)
                    existing.add(fake2)
                    pos_errors_added += 1
                    break

    random.shuffle(spectrum)
    return dna, spectrum, dna[:k], len(to_remove), pos_errors_added

if __name__ == "__main__":
    random.seed()
    N_INSTANCES = 1
    n = 700
    k = 8
    num_neg_errors = 70
    num_pos_errors = 0
    has_repeats = True

    with open(f"instance-{N_INSTANCES}-{n}.txt", "w") as inst_file, open(f"original-{N_INSTANCES}-{n}.txt", "w") as orig_file:
        for idx in range(N_INSTANCES):
            dna, spectrum, start_oligo, neg_errors_actual, pos_errors_actual = generate_single_instance(
                n, k, num_neg_errors, num_pos_errors, has_repeats
            )
            # Zapis instancji do wspólnego pliku
            inst_file.write(f"# Instance {idx+1}\n")
            inst_file.write(f"{n}\n")
            inst_file.write(f"{k}\n")
            inst_file.write(f"{start_oligo}\n")
            inst_file.write(f"{neg_errors_actual}\n")
            inst_file.write(f"{int(has_repeats)}\n")
            inst_file.write(f"{pos_errors_actual}\n")
            for kmer in spectrum:
                inst_file.write(kmer + "\n")
            inst_file.write("\n")

            # Zapis DNA do osobnego pliku
            orig_file.write(f"# Instance {idx+1}\n{dna}\n")

    print(f"Wygenerowano {N_INSTANCES} instancji.")
    print(f"Pliki: instance-{N_INSTANCES}-{n}.txt, original-{N_INSTANCES}-{n}.txt")
