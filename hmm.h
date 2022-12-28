#pragma once

#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <array>
#include <fstream>
#include <set>   

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
    double _transition_probabilities[N][N]; // N x N    // a
    double _emission_probabilities[N][M]; // N x M      // b
    double _states_probabilities[N]; // N               // initial_distribution

// emission to emission_index_map
    std::unordered_map<char, int> emission_to_idx;
    std::unordered_map<char, int> state_to_idx;

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

    //TODO: these are never used?
    double **get_transition_probabilities() {
        return _transition_probabilities;
    }

    double **get_emission_probabilities() {
        return _emission_probabilities;
    }

    const double *get_states_probabilities() const {
        return _states_probabilities;
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

    void create_default_emission_to_idx_map(){
        emission_to_idx['A'] = 0;
        emission_to_idx['C'] = 1;
        emission_to_idx['G'] = 2;
        emission_to_idx['T'] = 3;
    }

    static std::unordered_map<char, int> create_emission_to_idx_map(const std::string& emissions) {
        std::set<char> set_of_emmisions;
        for(char c : emissions) 
            set_of_emmisions.insert(c);
        
        //for(char c : set_of_emmisions) 
            //std::cout << c << " ";

        int index = 0;
        std::unordered_map<char, int> emission_to_idx;
        for(char c : set_of_emmisions) 
            emission_to_idx[c] = index++;
        
        for(char c : set_of_emmisions) 
            std::cout << c << " " << emission_to_idx[c] << std::endl;

        return emission_to_idx;
    }

    static std::unordered_map<char, int> create_state_to_idx_map(const std::string& states) {
        std::set<char> set_of_states;
        for(char c : states) 
            set_of_states.insert(c);
        
        //for(char c : set_of_emmisions) 
            //std::cout << c << " ";

        int index = 0;
        std::unordered_map<char, int> state_to_idx;
        for(char c : set_of_states) 
            state_to_idx[c] = index++;
        
        for(char c : set_of_states) 
            std::cout << c << " " << state_to_idx[c] << std::endl;

        return state_to_idx;
    }

    float **forward(const std::string& data) {
        int data_length = data.length(); // T
        // data // V
        // _transition_probabilities // a
        // _emission_probabilities // b
        // _states_probabilities // initial_distribution
        float** alpha = 0;
        alpha = new float*[data_length]; // [data_length][N]
        for(int i = 0; i < data_length; i++) {
            alpha[i] = new float[N];
            for(int j = 0; j < N; j++) 
                alpha[i][j] = 0;
        }

        int V[data_length];
        for(int i = 0; i < data_length; i++) V[i] = emission_to_idx[data[i]];

        for(int j = 0; j < N; j++) {
            alpha[0][j] = _states_probabilities[j] * _emission_probabilities[j][V[0]];
            // std::cout << _states_probabilities[j] * _emission_probabilities[j][V[0]] << std::endl;
        }

        for(int i = 1; i < data_length; i++) {
            for(int j = 0; j < N; j++) {
                float result = 0;
                for(int k = 0; k < N; k++) 
                    result += alpha[i-1][k] * _transition_probabilities[k][j];
                
                result *= _emission_probabilities[j][V[i]];
                alpha[i][j] = result;
            }
        }
        std::cout << "alpha" << std::endl;
        for(int i = 0; i < data_length; i++) {
            for(int j = 0; j < N; j++) {
                std::cout << alpha[i][j] << " ";
            }
            std::cout << std::endl;
        }

        return alpha;
    }

    float **backward(const std::string& data) {
        int data_length = data.length(); // T
        // data // V
        // _transition_probabilities // a
        // _emission_probabilities // b

        float** beta = 0;
        beta = new float*[data_length]; // [data_length][N]
        for(int i = 0; i < data_length; i++) {
            beta[i] = new float[N];
            for(int j = 0; j < N; j++) 
                beta[i][j] = 0;
        }

        int V[data_length];
        for(int i = 0; i < data_length; i++) V[i] = emission_to_idx[data[i]];

        for(int j = 0; j < N; j++) {
            beta[data_length-1][j] = 1;
            // wstd::cout << _states_probabilities[j] * _emission_probabilities[j][V[0]] << std::endl;
        }

        //std::cout << "beta calc" << std::endl;
        for(int i = data_length-2; i >= 0; i--) {
            for(int j = 0; j < N; j++) {
                float result = 0;
                float* result_matrix = new float[N];
                for(int k = 0; k < N; k++) {
                    result_matrix[k] = beta[i+1][k] * _emission_probabilities[k][V[i+1]];
                    //std::cout << beta[i+1][k] << " * " << _emission_probabilities[k][V[i+1]] << " = " << result_matrix[k] << ", ";
                }
                //std::cout << std::endl;
                
                for(int k = 0; k < N; k++) 
                    result += result_matrix[k] * _transition_probabilities[j][k];

                beta[i][j] = result;
            }
        }
        std::cout << "beta" << std::endl;
        for(int i = 0; i < data_length; i++) {
            for(int j = 0; j < N; j++) {
                std::cout << beta[i][j] << " ";
            }
            std::cout << std::endl;
        }

        return beta;
    }

    float forward_parameter(int t, int i) {
        // probability that model at time t is in state i, slide 6 (HMM2)
        return 0;
    }

    float backward_parameter(int t, int i) {
        // probability that model at time t is in state i and that it emitted T-t symbols, slide 4 (HMM2)
        return 0;
    }

    void viterbi_algorithm() {
        std::vector<int> DNA = {}; // DNA sequence

        double V[DNA.size()][N]; // Dynamic programming table
        int backtrack[DNA.size()][N]; // Backtracking table (stores pointers to previous state)


        // Initialize first column of dynamic programming table
        for (int i = 0; i < N; i++){
            V[0][i] = _states_probabilities[i] * _emission_probabilities[i][DNA[0]];
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
