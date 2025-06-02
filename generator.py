import random

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
    Generuje instancję problemu SBH z błędami negatywnymi i pozytywnymi (może być mieszaną)
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
    # 1. Wygeneruj losową sekwencję DNA długości n
    dna = ''.join(random.choices('ATCG', k=n))

    # 2. Utwórz idealne spektrum (kolejne k-mery)
    ideal_kmers = [dna[i : i + k] for i in range(n - k + 1)]
    # Dla szybkiego sprawdzania przynależności k-merów:
    ideal_set = set(ideal_kmers)

    # 3. Skopiuj idealne k-mery i usuń num_neg_errors (błędy negatywne)
    spectrum = ideal_kmers.copy()
    # Nie usuwaj więcej, niż jest dostępnych:
    to_remove = random.sample(spectrum, min(num_neg_errors, len(spectrum)))
    for kmer in to_remove:
        spectrum.remove(kmer)

    # 4. Dodaj num_pos_errors k-merów fałszywych (błędy pozytywne)
    #    Generujemy losowe k-mery aż uzyskamy odpowiednią liczbę, upewniając się,
    #    że nie należą do ideal_set i nie powtarzają się w 'spectrum' (unika duplikatów).
    existing = set(spectrum)
    while num_pos_errors > 0:
        candidate = ''.join(random.choices('ATCG', k=k))
        if (candidate not in ideal_set) and (candidate not in existing):
            spectrum.append(candidate)
            existing.add(candidate)
            num_pos_errors -= 1

    # 5. Wymieszaj kolejność k-merów (SBH nie zwraca uporządkowanego spektrum)
    random.shuffle(spectrum)

    # 6. Przygotuj pozostałe parametry instancji
    start_oligo = dna[:k]
    # Po usunięciu negatywnych i dodaniu pozytywnych, spectrum to lista, a:
    # num_neg_errors pozostało takie, jakie było na wejściu,
    # num_pos_errors został wyzerowany powyżej w pętli.

    # 7. Zapisz instancję do pliku instance_filename
    with open(instance_filename, "w") as f:
        f.write(f"{n}\n")
        f.write(f"{k}\n")
        f.write(f"{start_oligo}\n")
        f.write(f"{len(to_remove)}\n")      # faktycznie usuniętych negatywnych
        f.write(f"{has_repeats}\n")
        f.write(f"{len(spectrum) - (n - k + 1 - len(to_remove))}\n")  # liczba pozytywnych błędów: finalne_size - (ideal - neg)
        for kmer in spectrum:
            f.write(kmer + "\n")

    # 8. Zapisz oryginalne DNA do pliku original_filename
    with open(original_filename, "w") as f:
        f.write(dna)

    print(f"Instancja zapisana do '{instance_filename}'")
    print(f"Oryginalne DNA zapisane do '{original_filename}'")


# Przykładowe wywołanie generatora:
random.seed(42)  # Ustawienie ziarna dla powtarzalności wyników
generate_sbh_instance(
    n=500,
    k=8,
    num_neg_errors=50,
    num_pos_errors=0,
    has_repeats=True,
    instance_filename="instance.txt",
    original_filename="original.txt"
)
