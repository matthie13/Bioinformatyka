#include "instance_runner.h"

std::vector<SBHInstance> read_instances(const std::string& filename) {
    std::vector<SBHInstance> instances;
    std::ifstream file{filename};

    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        SBHInstance instance;

        // Read n
        instance.n = std::stoull(line);

        // Read k
        if (!std::getline(file, line)) {
            std::cerr << "ERROR: Unexpected EOF reading k" << std::endl;
            exit(EXIT_FAILURE);
        }
        instance.k = std::stoi(line);

        // Read start_oligo
        if (!std::getline(file, line)) {
            std::cerr << "ERROR: Unexpected EOF reading start_oligo" << std::endl;
            exit(EXIT_FAILURE);
        }
        instance.start_oligo = line;

        // Read neg_errors
        if (!std::getline(file, line)) {
            std::cerr << "ERROR: Unexpected EOF reading neg_errors" << std::endl;
            exit(EXIT_FAILURE);
        }
        instance.neg_errors = std::stoi(line);

        // Read has_repeats
        if (!std::getline(file, line)) {
            std::cerr << "ERROR: Unexpected EOF reading has_repeats" << std::endl;
            exit(EXIT_FAILURE);
        }
        instance.has_repeats = (std::stoi(line) == 1);

        // Read pos_errors
        if (!std::getline(file, line)) {
            std::cerr << "ERROR: Unexpected EOF reading pos_errors" << std::endl;
            exit(EXIT_FAILURE);
        }
        instance.pos_errors = std::stoi(line);

        // Read spectrum (k-mers until empty line or EOF)
        instance.spectrum.clear();
        while (std::getline(file, line)) {
            if (line.empty()) {
                break;
            }
            instance.spectrum.push_back(line);
        }
        instances.push_back(instance);
    }

    file.close();
    return instances;
}

std::vector<std::string> read_original_sequences(const std::string& filename) {
    std::vector<std::string> sequences;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip comments
        if (line.empty() || line[0] == '#') {
            continue;
        }
        sequences.push_back(line);
    }

    file.close();
    return sequences;
}
