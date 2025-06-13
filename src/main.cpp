#include "config.h"
#include "instance_runner.h"
#include "sbh.h"
#include "levenshtein.h"
#include "ant.h"


int main(const int argc, char** argv) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <instance_file> <original_file>\n";
        std::cout << "Example: " << argv[0] << " instance-50-300.txt original-50-300.txt\n";
        return EXIT_FAILURE;
    }

    auto instance_file = std::string(argv[1]);
    auto original_file = std::string(argv[2]);

    std::cout << "INFO: Reading instances from: " << instance_file << "\n";
    auto instances = read_instances(instance_file);

    std::cout << "INFO: Reading original sequences from: " << original_file << "\n";
    auto originals = read_original_sequences(original_file);

    if (instances.size() != originals.size()) {
        std::cerr << "ERROR: Number of instances (" << instances.size()
                  << ") doesn't match number of original sequences (" << originals.size() << ")\n";
        return EXIT_FAILURE;
    }

    std::cout << "INFO: Processing " << instances.size() << " instances...\n\n";

    int total_levenshtein = 0;

    for (size_t i = 0; i < instances.size(); ++i) {
        const auto& instance = instances[i];
        const auto& original = originals[i];

        std::cout << "=== Instance " << (i + 1) << " ===\n";
        std::cout << "Target length: " << instance.n << "\n";
        std::cout << "K-mer length: " << instance.k << "\n";
        std::cout << "Spectrum size: " << instance.spectrum.size() << "\n";
        std::cout << "Start oligo: " << instance.start_oligo << "\n";
        std::cout << "Negative errors: " << instance.neg_errors << "\n";
        std::cout << "Positive errors: " << instance.pos_errors << "\n";
        std::cout << "Has repeats: " << (instance.has_repeats ? "Yes" : "No") << "\n";

        std::string reconstructed = run_aco(
            instance.spectrum,
            instance.k,
            instance.n,
            instance.start_oligo,
            instance.neg_errors,
            instance.has_repeats,
            instance.pos_errors,
            100,
            1.0,
            2.0,
            0.1,
            20.0,
            1.0,
            0,
            100
        );
        int levenshtein = levenshtein_score(reconstructed, original);
        total_levenshtein += levenshtein;

        std::cout << "Original:      " << original.substr(0, 50) << "...\n";
        std::cout << "Reconstructed: " << reconstructed.substr(0, 50) << "...\n";
        std::cout << "Levenshtein score: " << levenshtein << '\n';
    }

    std::cout << "=== SUMMARY ===\n";
    std::cout << "Total instances: " << instances.size() << "\n";
    std::cout << "Average Levenshtein score: " << (total_levenshtein / instances.size()) << '\n';

    return EXIT_SUCCESS;
}
