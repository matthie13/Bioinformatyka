import random
from collections import Counter

def generate_sbh_instance(
    n,
    k,
    num_neg_errors,
    num_pos_errors,
    has_repeats,
    instance_filename="instance.txt",
    original_filename="original.txt"
):
    """
    Generuje instancję problemu SBH z błędami negatywnymi i pozytywnymi (+ mieszaną)
    i zapisuje do plików:
      - instance_filename: parametry oraz spektrum (spectrum)
      - original_filename: oryginalną sekwencję DNA

    Parametry:
    - n: długość DNA
    - k: długość oligonukleotydów (k-merów)
    - num_neg_errors: liczba k-merów do usunięcia (błędy negatywne)
    - num_pos_errors: liczba k-merów do dodania (błędy pozytywne)
    - has_repeats: bool, czy negatywne błędy mogą wynikać z powtórzeń w DNA
    - instance_filename: nazwa pliku z instancją
    - original_filename: nazwa pliku z oryginalnym DNA
    """
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

        #Jeżeli długość jeszcze < n, dopełnij losowymi nukleotydami
        while len(dna_list) < n:
            dna_list.append(random.choice('ATCG'))

        dna = ''.join(dna_list[:n])

    else:
        dna = ''.join(random.choices('ATCG', k=n))
    all_kmers = [dna[i : i+k] for i in range(n - k + 1)]
    ideal_set = set(all_kmers)  # unikalne k-mery
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

        #Usuń wybrane k-mery
        for kmer in to_remove:
            if kmer in spectrum:
                spectrum.remove(kmer)

    #Generowanie błędów pozytywnych
    existing = set(spectrum)
    ideal_list_for_pos = list(ideal_set)  #baza do tworzenia fałszywek
    pos_errors_added = 0
    idx = 0

    while pos_errors_added < num_pos_errors and idx < len(ideal_list_for_pos):
        base = ideal_list_for_pos[idx]
        idx += 1

        #zamiana ostatniego nukleotydu
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

        #zamiana nukleotydu środkowego
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
    #Mieszamy spektrum (żeby symulować, że nie znamy kolejności k-merów)
    random.shuffle(spectrum)

    start_oligo = dna[:k]
    with open(instance_filename, "w") as f:
        f.write(f"{n}\n")
        f.write(f"{k}\n")
        f.write(f"{start_oligo}\n")
        f.write(f"{len(to_remove)}\n")
        f.write(f"{int(has_repeats)}\n")
        f.write(f"{pos_errors_added}\n")
        for kmer in spectrum:
            f.write(kmer + "\n")

    with open(original_filename, "w") as f:
        f.write(dna)

    print(f"Instancja zapisana do '{instance_filename}'")
    print(f"Oryginalne DNA zapisane do '{original_filename}'")

if __name__ == "__main__":
    random.seed(42)  # opcjonalnie dla powtarzalności
    generate_sbh_instance(
        n=700,
        k=8,
        num_neg_errors=5,
        num_pos_errors=0,
        has_repeats=True,
        instance_filename="instance.txt",
        original_filename="original.txt"
    )
