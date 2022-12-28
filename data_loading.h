#pragma once

#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <array>
#include <fstream>

std::pair<std::vector<std::pair<int, int>>, std::string> load_data(
        const std::string& islands_path, const std::string& sequence_path
) {
    std::ifstream islands_file(islands_path);
    if (!islands_file.good()) {
        std::cerr << "Path to islands " << islands_path << " is not valid" << std::endl;
        return {{},{}};
    }

    std::string line;

    std::vector<std::pair<int, int>> islands;
    while (std::getline(islands_file, line)) {
        auto comma = line.find(',');
        islands.emplace_back(
                std::stoi(line.substr(0, comma)),
                std::stoi(line.substr(comma + 1))
        );
    }

    std::ifstream sequence_file(sequence_path);
    if (!sequence_file.good()) {
        std::cerr << "Path to sequence " << sequence_path << " is not valid" << std::endl;
        return {{},{}};
    }

    std::string sequence;

    while (std::getline(sequence_file, line)) {
        std::transform(line.begin(), line.end(), line.begin(),
                       [](auto c) { return std::toupper(c); });
        sequence.append(line);
    }

    return std::make_pair(islands, sequence);
}

std::string from_islands_to_str(
        const std::vector<std::pair<int,int>>& islands, const std::string& emissions, bool simple = true
) {
    // assumption: all pairs are ascending

    std::string states;

    int current_idx = 0;
    for (auto island: islands) {
        if (simple) {
            states.append(island.first - current_idx, '-');
            states.append(island.second - island.first + 1, '+');
        } else {
            std::transform(emissions.begin() + current_idx, emissions.begin() + island.first,
                           std::back_inserter(states),
                           [](auto c) { return std::tolower(c); });
            states.insert(island.first, emissions, island.first, island.second - island.first + 1);
        }
        current_idx = island.second + 1;
    }
    if (simple)
        states.append(emissions.length() - current_idx, '-');
    else
        std::transform(emissions.begin() + current_idx, emissions.end(), states.begin() + current_idx,
                       [](auto c) { return std::tolower(c); });
    return states;
}

std::pair<std::string, std::string> parse_data(
        const std::string& islands_path, const std::string& sequence_path, bool simple = true
) {
    std::vector<std::pair<int,int>> islands;
    std::string emissions;

    std::tie(islands, emissions) = load_data(islands_path, sequence_path);
    return std::make_pair(from_islands_to_str(islands,emissions, simple), emissions);
}