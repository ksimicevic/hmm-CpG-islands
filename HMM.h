#pragma once

#include <unordered_map>
#include <string>
#include <vector>

const static int N = 2;
const static int M = 4;

class hidden_markov_chain {
// known from the beginning
    std::vector<std::string> _states;
    std::vector<char> _observations;

// parameters that we need to figure out
    float _transition_probabilities[N][N]; // N x N
    float _emission_probabilities[N][M]; // N x M, N = len(states), M = len(emission)
    float _states_probabilities[N]; // N
public:
    hidden_markov_chain(const std::vector<std::string>& states, const std::vector<char>& observations) :
        _states(states), _observations(observations) {
    }

    void fit(const std::string& states, const std::string& data) {
        estimate_initial_probabilities();
    }

    //predict
    void decode(const std::string& data) {
        // call viterbi algorithm
    }

private:
    void estimate_initial_probabilities() {
        //estimate initial probabilities for transitions, emissions and states based on frequencies, slide 4 (HMM2)
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