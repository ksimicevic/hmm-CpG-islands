#include <iostream>
#include "HMM.h"

class Test {

public:
    static void test_estimate_initial_probabilities() {
        std::cout << ">>> Test estimate initial probabilities begin. <<<" << std::endl;

        const std::string states = "FBBFFBBFBFFBFFB";
        const std::string emissions = "562364166532654";

        hidden_markov_chain hmm({'F', 'B'}, {'1', '2', '3', '4', '5', '6'});
        hmm.estimate_initial_probabilities(states, emissions);

        double states_initial_probs[] = {8/15, 7/15};
        double transition_initial_probs[][2] = {
                {3/7, 5/7}, //F->F, F->B
                {4/7, 2/7} //B->F, B->B
        };
        double emission_initial_probs[][6] = {
                {0/8, 0/8, 2/8, 0/8, 3/8,3/8},
                {1/7, 2/7, 0/7, 2/7, 0/7, 2/7}
        };

        //TODO: recheck ground truth

        if (!is_equal_arrays(states_initial_probs, hmm._states_probabilities, 2))
            std::cerr << "Initial states probabilities are not equal!" << std::endl;

//        if (!is_equal_matrices(transition_initial_probs, hmm._transition_probabilities, 2, 2))
//            std::cerr << "Transition states probabilities are not equal!" << std::endl;
//
//        if (!is_equal_matrices(emission_initial_probs, hmm._emission_probabilities, 2, 6))
//            std::cerr << "Emission states probabilities are not equal!" << std::endl;

        std::cout << ">>> Test estimate initial probabilities end. <<<" << std::endl;
    }

private:
    static bool is_equal_arrays(double* a, double* b, int n) {
        for (auto i = 0; i < n; i++) {
            if (a[i] != b[i]) return false;
        }
        return true;
    }

    static bool is_equal_matrices(double** a, double** b, int n, int m) {
        //TODO: implementation
    }
};

int main() {
    Test::test_estimate_initial_probabilities();
    return 0;
}
