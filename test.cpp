#include <iostream>
#include <cmath>

#include "hmm.h"
#include "data_loading.h"

class Test {

public:
    static void test_create_emission_to_idx_map() {
        const auto islands_path = "data\\relaxed-examples\\mm39_1_islands.txt";
        const auto sequence_path = "data\\relaxed-examples\\mm39_1.txt";

        const int N = 2;
        const int M = 4;

        std::vector<std::pair<int, int>> islands;
        std::string sequence;
        std::tie(islands, sequence) = load_data(islands_path, sequence_path);

        hidden_markov_chain<N, M>::create_emission_to_idx_map(sequence);

        std::cout << ">>> Test create emission to idx map done. <<<" << std::endl;
    }

    static void test_forward() {
        std::cout << ">>> Test forward start. <<<" << std::endl;
        const std::string states = "FBBFFBBFBFFBFFB";
        const std::string emissions = "562364166532654";

        const int N = 2;
        const int M = 6;

        hidden_markov_chain<N, M> hmm({'F', 'B'}, {'1', '2', '3', '4', '5', '6'});
        hmm.estimate_initial_probabilities(states, emissions);

        hmm.emission_to_idx = hidden_markov_chain<N, M>::create_emission_to_idx_map(emissions);
        double ** alpha = hmm.forward(emissions);

        std::cout << "alpha" << std::endl;
        for(int i = 0; i < emissions.length(); i++) {
            for(int j = 0; j < N; j++) {
                std::cout << alpha[i][j] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << ">>> Test forward done. <<<" << std::endl;
    }

    static void test_backward() {
        std::cout << ">>> Test backward start. <<<" << std::endl;
        const std::string states = "FBBFFBBFBFFBFFB";
        const std::string emissions = "562364166532654";

        const int N = 2;
        const int M = 6;

        hidden_markov_chain<N, M> hmm({'F', 'B'}, {'1', '2', '3', '4', '5', '6'});
        hmm.estimate_initial_probabilities(states, emissions);

        hmm.emission_to_idx = hidden_markov_chain<N, M>::create_emission_to_idx_map(emissions);
        double ** beta = hmm.backward(emissions);

        std::cout << "beta" << std::endl;
        for(int i = 0; i < emissions.length(); i++) {
            for(int j = 0; j < N; j++) {
                std::cout << beta[i][j] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << ">>> Test backward done. <<<" << std::endl;
    }

    static void test_baum_welch_algorithm() {
        std::cout << ">>> Test baum_welch_algorithm start. <<<" << std::endl;
        const std::string states = "FBBFFBBFBFFBFFB";
        const std::string emissions = "562364166532654";

        const int N = 2;
        const int M = 6;

        hidden_markov_chain<N, M> hmm({'F', 'B'}, {'1', '2', '3', '4', '5', '6'});
        hmm.estimate_initial_probabilities(states, emissions);

        hmm.emission_to_idx = hidden_markov_chain<N, M>::create_emission_to_idx_map(emissions);
        hmm.baum_welch_algorithm(emissions, 1000);

        std::cout << ">>> Test baum_welch_algorithm done. <<<" << std::endl;
    }

    static void test_estimate_initial_probabilities() {
        const std::string states = "FBBFFBBFBFFBFFB";
        const std::string emissions = "562364166532654";

        const int N = 2;
        const int M = 6;

        hidden_markov_chain<N, M> hmm({'F', 'B'}, {'1', '2', '3', '4', '5', '6'});
        hmm.estimate_initial_probabilities(states, emissions);

        double states_initial_probs[] = {8/15., 7/15.};
        double transition_initial_probs[N][N] = {
                {3/7., 5/7.}, //F->F, F->B
                {4/7., 2/7.} //B->F, B->B
        };
        double emission_initial_probs[N][M] = {
                {0/8., 0/8., 2/8., 0/8., 3/8., 3/8.},
                {1/7., 2/7., 0/7., 2/7., 0/7., 2/7.}
        };

        if (!is_equal_arrays<N>(states_initial_probs, hmm._states_probabilities))
            std::cerr << "Initial states probabilities are not equal!" << std::endl;

        if (!is_equal_matrices<N, N>(transition_initial_probs, hmm._transition_probabilities))
            std::cerr << "Transition probabilities are not equal!" << std::endl;

        if (!is_equal_matrices<N, M>(emission_initial_probs, hmm._emission_probabilities))
            std::cerr << "Emission probabilities are not equal!" << std::endl;

        std::cout << ">>> Test estimate initial probabilities done. <<<" << std::endl;
    }

    static void test_load_data() {
        const auto islands_path = "data\\relaxed-examples\\mm39_1_islands.txt";
        const auto sequence_path = "data\\relaxed-examples\\mm39_1.txt";

        const int N = 2;
        const int M = 4;

        std::vector<std::pair<int, int>> islands;
        std::string sequence;
        std::tie(islands, sequence) = load_data(islands_path, sequence_path);

        std::cout << "islands" << std::endl;
        for (int i = 0; i < islands.size(); i++)
            std::cout << islands[i].first << ", " << islands[i].second << std::endl;


        std::tie(islands, sequence) = load_data(islands_path, sequence_path);
    }

    static void test_convert_islands_to_string() {
        const auto emissions = "ACTGCGCGCATTTGCGCTGCA";
        std::vector<std::pair<int, int>> islands;
        islands.emplace_back(3, 8); // inclusive borders
        islands.emplace_back(13, 20);

        const auto states_simple = "---++++++----++++++++";
        auto states_simple_ = from_islands_to_str(islands, emissions);

        if (states_simple != states_simple_)
            std::cerr << "Conversion from islands to string (simple) failed! Strings do not match." << std::endl;

        const auto states = "actGCGCGCatttGCGCTGCA";
        auto states_ = from_islands_to_str(islands, emissions, false);

        if (states != states_)
            std::cerr << "Conversion from islands to string failed! Strings do not match." << std::endl;

        std::cout << ">>> Test convert islands to string done. <<<" << std::endl;
    }

private:
    template<int N>
    static bool is_equal_arrays(double a[N], double b[N]) {
        for (auto i = 0; i < N; ++i) {
            if (std::abs(a[i] - b[i]) > std::numeric_limits<double>::epsilon()) return false;
        }
        return true;
    }

    template<int N, int M>
    static bool is_equal_matrices(double a[N][M], double b[N][M]) {
        for (auto i = 0; i < N; ++i) {
            if (!is_equal_arrays<M>(a[i], b[i])) return  false;
        }

        return true;
    }
};

int main() {
    Test::test_estimate_initial_probabilities();
    Test::test_load_data();
    Test::test_convert_islands_to_string();
    Test::test_create_emission_to_idx_map();
    //Test::test_forward();
    //Test::test_backward();
    Test::test_baum_welch_algorithm();
    return 0;
}