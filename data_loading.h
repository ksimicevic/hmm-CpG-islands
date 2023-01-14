#pragma once

#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <array>
#include <fstream>

std::vector<std::pair<int,int>> load_island(const std::string& islands_path) {
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

    return islands;
}

std::string load_sequence(const std::string& sequence_path) {
    std::ifstream sequence_file(sequence_path);
    if (!sequence_file.good()) {
        std::cerr << "Path to sequence " << sequence_path << " is not valid" << std::endl;
        return {{},{}};
    }

    std::string sequence, line;
    while (std::getline(sequence_file, line)) {
        std::transform(line.begin(), line.end(), line.begin(),
                       [](auto c) { return std::toupper(c); });
        sequence.append(line);
    }

    return sequence;
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
        std::transform(emissions.begin() + current_idx, emissions.end(), std::back_inserter(states),
                       [](auto c) { return std::tolower(c); });
    return states;
}

std::string from_str_to_islands(const std::string& str, bool simple = true) {
    std::string ret;

    if (simple) {
        for (size_t pos = 0; pos != std::string::npos && pos < str.length(); pos++) {
            auto start = str.find('+', pos);
            if (start != std::string::npos) {
                auto end = str.find('-', start + 1);
                if (end == std::string::npos)
                    end = str.length() - 1;

                ret.append(std::to_string(start))
                        .append(",")
                        .append(std::to_string(end))
                        .append("\n");
                pos = end;
            } else {
                break;
            }
        }
    } else {
        for (size_t pos = 0; pos != std::string::npos && pos < str.length(); pos++) {
            if (std::isupper(str.at(pos))){ // start of an island
                auto start = pos;
                while (std::isupper(str.at(pos)) && pos < str.length() - 1) pos++;
                auto end = pos;

                ret.append(std::to_string(start))
                        .append(",")
                        .append(std::to_string(end))
                        .append("\n");
                pos = end;
            }
        }
    }

    return ret;
}

std::pair<std::string, std::string> parse_data(
        const std::string& islands_path, const std::string& sequence_path, bool simple = true
) {
    auto islands = load_island(islands_path);
    auto sequence = load_sequence(sequence_path);

    return std::make_pair(from_islands_to_str(islands,sequence, simple), sequence);
}
