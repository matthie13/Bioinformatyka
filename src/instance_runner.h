#ifndef INSTANCE_RUNNER_H
#define INSTANCE_RUNNER_H

#include "config.h"
#include "sbh.h"


std::vector<SBHInstance> read_instances(const std::string& filename);
std::vector<std::string> read_original_sequences(const std::string& filename);

#endif //INSTANCE_RUNNER_H
