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
        std::vector<int> DNA = {}; // DNA sequence

        double V[DNA.size()][N]; // Dynamic programming table
        int backtrack[DNA.size()][N]; // Backtracking table (stores pointers to previous state)


        // Initialize first column of dynamic programming table
        for (int i = 0; i < N; i++){
            V[0][i] = _states_probabilities[i] * _emission_probabilities[i][DNA[0]]
            backtrack[0][i] = -1;
        }

        // Fill in rest of dynamic programming table
        for (int t = 1; t < DNA.size(); t++){
            for (int j = 0; j < N; j++){
                double max_prob = 0;
                int max_state = -1;
                for (int i = 0; i < N; i++){
                    double prob = V[t-1][i] * _transition_probabilities[i][j] * _emission_probabilities[j][DNA[t]];
                    if (prob > max_prob){
                        max_prob = prob;
                        max_state = i;
                    }
                }
                V[t][j] = max_prob;
                backtrack[t][j] = max_state;
            }
        }

        // Find state with highest probability in final column of dynamic programming table
        double max_prob = 0;
        int max_state = -1;
        for (int i = 0; i < N; i++){
            if (V[DNA.size()-1][i] > max_prob){
                max_prob = V[DNA.size()-1][i];
                max_state = i;
            }
        }

        // Backtrack through dynamic programming table to reconstruct most likely sequence of states
        std::vector<int> path;
        int state = max_state;
        for (int t = DNA.size()-1; t >= 0; t--){
            path.push_back(state);
            state = backtrack[t][state];
        }

        // Reverse the sequence to get the most likely sequence in the correct order
        return std::reverse(path.begin(), path.end());
    }

    void baum_welch_algorithm() {
        // baum welch algorithm, parameters and placement pending (it might accept the entire hmm as whole)
    }
};
