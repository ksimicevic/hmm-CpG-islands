#pragma once

#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <array>
#include <fstream>

// forward declaration
class Test;

// N = len(states), M = len(symbols)
template<int N, int M>
class hidden_markov_chain {
    friend class Test;

// known from the beginning
    std::vector<char> _states;
    std::vector<char> _symbols;

// parameters that we need to figure out
    double _transition_probabilities[N][N]; // N x N
    double _emission_probabilities[N][M]; // N x M
    double _states_probabilities[N]; // N
public:
    hidden_markov_chain(const std::vector<char>& states, const std::vector<char>& symbols) :
            _states(states), _symbols(symbols) {
    }

    void fit(const std::string& states, const std::string& emissions) {
        estimate_initial_probabilities(states, emissions);
    }

    //predict
    void decode(const std::string& data) {
        // call viterbi algorithm
    }

private:
    //TODO: convert std::vector<std::pair<int,int>> to string representing states
    static std::pair<std::vector<std::pair<int, int>>, std::string> load_data(
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

    //TODO: this fun behaves differently depending on simple or complicated model
    static std::string from_islands_to_str(
            const std::vector<std::pair<int,int>>& islands, const std::string& emissions, bool simple = true) {
        // assumption: all pairs are ascending

        std::string states;

        //TODO: improve?
        int current_idx = 0;
        for (auto island: islands) {
            states.append(island.first - current_idx, '-');
            states.append(island.second - island.first + 1, '+');
            current_idx = island.second + 1;
        }
        states.append(emissions.length() - current_idx, '-');

        return states;
    }

    void estimate_initial_probabilities(const std::string& states, const std::string& emissions) {
        //estimate initial probabilities for transitions, emissions and states based on frequencies, slide 4 (HMM2)
        if (states.length() != emissions.length()) {
            std::cerr << "Given states and corresponding emissions MUST match in length." << std::endl;
            return;
        }

        // INITIAL STATES PROBABILITIES

        // count all states occurrences
        std::unordered_map<char, int> states_freqs;
        for (char state : states)
            states_freqs[state] += 1;

        // estimate initial states probabilities as freq(state) / number_of_emissions
        for (auto i = 0; i < _states.size(); ++i)
            _states_probabilities[i] = (double) states_freqs[_states[i]] / (double) states.length();


        // INITIAL TRANSITIONS PROBABILITIES

        // count all transitions
        std::unordered_map<std::string, int> transition_freqs;
        for (auto i = 0; i < states.length() - 1; ++i)
            transition_freqs[states.substr(i, 2)] += 1;

        // put frequencies of transitions into matrix
        for (auto i = 0; i < _states.size(); ++i) {
            for (auto j = 0; j < _states.size(); ++j) {
                auto from_state = _states[i];
                auto to_state = _states[j];
                std::string transition {from_state, to_state};
                _transition_probabilities[i][j] = (float) transition_freqs[transition];
            }
        }

        // make probabilities of frequencies found in matrix
        for (auto i = 0; i < _states.size(); ++i) {
            double col_sum = 0;
            for (auto j = 0; j < _states.size(); ++j) col_sum += _transition_probabilities[j][i];
            for (auto j = 0; j < _states.size(); ++j) _transition_probabilities[j][i] /= col_sum;
        }

        // INITIAL EMISSIONS PROBABILITIES
        std::unordered_map<char, std::vector<int>> emission_freqs;
        std::unordered_map<char, int> state_to_idx;
        for (char _symbol : _symbols)
            for (auto j = 0; j < _states.size(); ++j) {
                emission_freqs[_symbol].resize(_states.size());
                state_to_idx[_states[j]] = j;
            }

        for (auto i = 0; i < emissions.length(); ++i) {
            auto state = states[i];
            auto emission = emissions[i];

            emission_freqs[emission][state_to_idx[state]] += 1;
        }

        for (auto i = 0; i < _states.size(); ++i) {
            for (auto j = 0; j < _symbols.size(); ++j)
                _emission_probabilities[i][j] =
                        (double) emission_freqs[_symbols[j]][i] / (double) states_freqs[_states[i]];
        }
    }

    float forward_parameter(int t, int i) {
        // probability that model at time t is in state i, slide 6 (HMM2)
    }

    float backward_parameter(int t, int i) {
        // probability that model at time t is in state i and that it emitted T-t symbols, slide 4 (HMM2)
    }

    void viterbi_algorithm() {

    }

    void baum_welch_algorithm() {
        // baum welch algorithm, parameters and placement pending (it might accept the entire hmm as whole)
    }
};