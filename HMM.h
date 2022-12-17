#pragma once

#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <array>

const static int N = 2;
const static int M = 4;

// forward declaration
class Test;

class hidden_markov_chain {
    friend class Test;

// known from the beginning
    std::vector<char> _states;
    std::vector<char> _symbols;

// parameters that we need to figure out
    double _transition_probabilities[N][N]; // N x N
    double _emission_probabilities[N][M]; // N x M, N = len(states), M = len(emission)
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
        for (char state : states) {
            states_freqs[state] += 1;
        }

//        std::for_each(states.cbegin(), states.cbegin(),
//              [&states_freqs](auto state) {
//                    states_freqs[state] += 1;
//        });

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
            double row_sum = 0;
            for (auto j = 0; j < _states.size(); ++j) row_sum += _transition_probabilities[i][j];
            for (auto j = 0; j < _states.size(); ++j) _transition_probabilities[i][j] /= row_sum;
        }

        // INITIAL EMISSIONS PROBABILITIES
        std::unordered_map<char, int> emission_freqs;
        std::for_each(emissions.cbegin(), emissions.cend(),
              [&emission_freqs] (auto emission) {
                    emission_freqs[emission] += 1;
        });

        for (auto i = 0; i < _states.size(); ++i) {
            for (auto j = 0; j < emission_freqs.size(); ++j)
                _emission_probabilities[i][j] = (double) emission_freqs[emissions[j]] / (double) states_freqs[_states[i]];
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
