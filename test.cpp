#include <iostream>
#include <cmath>

#include "hmm.h"
#include "data_loading.h"

class Test {

public:
//    static void test_create_emission_to_idx_map() {
//        const auto islands_path = "data\\relaxed-examples\\mm39_1_islands.txt";
//        const auto sequence_path = "data\\relaxed-examples\\mm39_1.txt";
//
//        const int N = 2;
//        const int M = 4;
//
//        std::string sequence = load_sequence(sequence_path);
//
//        hidden_markov_chain<N, M>::create_emission_to_idx_map(sequence);
//
//        std::cout << ">>> Test create emission to idx map done. <<<" << std::endl;
//    }

    static void test_forward() {
        std::cout << ">>> Test forward start. <<<" << std::endl;
        const std::string states = "FBBFFBBFBFFBFFB";
        const std::string emissions = "562364166532654";

        const int N = 2;
        const int M = 6;

        hidden_markov_chain<N, M> hmm({'F', 'B'}, {'1', '2', '3', '4', '5', '6'});
        hmm.estimate_initial_probabilities(states, emissions);

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

        // hmm.baum_welch_algorithm(emissions, 1000);

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

    static void test_viterbi(){
        std::cout << ">>> Test viterbi start. <<<" << std::endl;
        
        const std::string states = "FBBFFBBFBFFBFFB";
        const std::string emissions = "562364166532654";

        const int N = 2;
        const int M = 6;

        hidden_markov_chain<N, M> hmm({'F', 'B'}, {'1', '2', '3', '4', '5', '6'});

        hmm.estimate_initial_probabilities(states, emissions);

        std::string result = hmm.viterbi_algorithm(states);

        std::cout << result << std::endl;

        std::cout << ">>> Test viterbi done. <<<" << std::endl;
    }

    static void test_mm39_relaxed_simple_without_baum_welch_algorithm() {
        std::cout << ">>> Test mm39_relaxed_simple start. <<<" << std::endl;
        const std::string train_seq_path = "data\\sequences\\mm39_1.txt";
        const std::string train_path = "data\\relaxed-examples\\mm39_1_islands.txt";

        const std::string test_seq_path = "data\\sequences\\mm39_2.txt";
        const std::string test_path = "data\\relaxed-examples\\mm39_2_islands.txt";

        std::string train_sequence, train_emissions, test_sequence, test_emissions;
        std::tie(train_sequence, train_emissions) = parse_data(train_path, train_seq_path, true);
        std::tie(test_sequence, test_emissions) = parse_data(test_path, test_seq_path, true);

        hidden_markov_chain<2, 4> hmm({'+', '-'}, {'A', 'C', 'G', 'T'});
        hmm.estimate_initial_probabilities(train_sequence, train_emissions);

        auto predicted_emissions = hmm.predict(test_sequence);

        std::string hit_or_miss;
        double accuracy;
        std::tie(accuracy, hit_or_miss) = hmm.evaluate(test_sequence, predicted_emissions);
        std::cout << "Accuracy: " << accuracy << std::endl;
        //std::cout << "Hit-or-miss: " << hit_or_miss << std::endl;
        
        std::cout << ">>> Test mm39_relaxed_simple done. <<<" << std::endl;
    }

    static void test_mm39_relaxed_simple() {
        std::cout << ">>> Test mm39_relaxed_simple start. <<<" << std::endl;
        const std::string train_seq_path = "data\\relaxed-examples\\mm39_1.txt";
        const std::string train_path = "data\\relaxed-examples\\mm39_1_islands.txt";

        const std::string test_seq_path = "data\\relaxed-examples\\mm39_2.txt";
        const std::string test_path = "data\\relaxed-examples\\mm39_2_islands.txt";

        std::string train_sequence, train_emissions, test_sequence, test_emissions;
        std::tie(train_sequence, train_emissions) = parse_data(train_path, train_seq_path, true);
        std::tie(test_sequence, test_emissions) = parse_data(test_path, test_seq_path, true);

        hidden_markov_chain<2, 4> hmm({'+', '-'}, {'A', 'C', 'G', 'T'});
        hmm.fit(train_sequence, train_emissions, 300, 1);

        auto predicted_emissions = hmm.predict(test_sequence);

        std::string hit_or_miss;
        double accuracy;
        std::tie(accuracy, hit_or_miss) = hmm.evaluate(test_sequence, predicted_emissions);

        int plus = 0;
        int minus = 0;

        for(int i = 0; i < train_sequence.length(); i++) {
            if(train_sequence.at(i) == '+') plus++;
            else minus++;
        }

        std::cout << "Plus and minus counts:" << std::endl;
        std::cout << plus << " " << minus << std::endl;
        
        // std::cout << hit_or_miss << std::endl;
        int c_matrix[2][2] = {{0, 0}, 0, 0};

        for(int i = 0; i < test_sequence.length(); i++) {
            char real_c = test_sequence.at(i);
            char predicted_c = predicted_emissions.at(i);

            if(real_c == '+') {
                if(predicted_c == '+') c_matrix[0][0]++;
                else c_matrix[1][0]++;
            }
            else {
                if(predicted_c == '-') c_matrix[1][1]++;
                else c_matrix[0][1]++;
            }
        }
        std::cout << "Accuracy: " << accuracy << std::endl;
        //std::cout << "Hit-or-miss: " << hit_or_miss << std::endl;

        std::cout << "Confusion matrix:" << std::endl;
        std::cout << c_matrix[0][0] << " " << c_matrix[0][1] << std::endl;
        std::cout << c_matrix[1][0] << " " << c_matrix[1][1] << std::endl;
        
        std::cout << std::endl;
        
        std::cout << ">>> Test mm39_relaxed_simple done. <<<" << std::endl;
    }

    static void test_chr19_0_relaxed_simple() {
        std::cout << ">>> Test chr19_0_relaxed_simple start. <<<" << std::endl;
        const std::string train_seq_path = "even\\sequences\\chr19_0.txt";
        const std::string train_path = "even\\islands\\chr19_0.txt";

        const std::string test_seq_path = "even\\sequences\\chr19_3.txt";
        const std::string test_path = "even\\islands\\chr19_3.txt";

        std::string train_sequence, train_emissions, test_sequence, test_emissions;
        std::tie(train_sequence, train_emissions) = parse_data(train_path, train_seq_path, true);
        std::tie(test_sequence, test_emissions) = parse_data(test_path, test_seq_path, true);

        hidden_markov_chain<2, 4> hmm({'+', '-'}, {'A', 'C', 'G', 'T'});
        hmm.fit(train_sequence, train_emissions, 300, 10);

        auto predicted_emissions = hmm.predict(test_sequence);

        std::string hit_or_miss;
        double accuracy;
        std::tie(accuracy, hit_or_miss) = hmm.evaluate(test_sequence, predicted_emissions);

        int plus = 0;
        int minus = 0;

        for(int i = 0; i < train_sequence.length(); i++) {
            if(train_sequence.at(i) == '+') plus++;
            else minus++;
        }

        std::cout << "Plus and minus counts:" << std::endl;
        std::cout << plus << " " << minus << std::endl;
        
        // std::cout << hit_or_miss << std::endl;
        int c_matrix[2][2] = {{0, 0}, 0, 0};

        for(int i = 0; i < test_sequence.length(); i++) {
            char real_c = test_sequence.at(i);
            char predicted_c = predicted_emissions.at(i);

            if(real_c == '+') {
                if(predicted_c == '+') c_matrix[0][0]++;
                else c_matrix[1][0]++;
            }
            else {
                if(predicted_c == '-') c_matrix[1][1]++;
                else c_matrix[0][1]++;
            }
        }
        std::cout << "Accuracy: " << accuracy << std::endl;
        //std::cout << "Hit-or-miss: " << hit_or_miss << std::endl;

        std::cout << "Confusion matrix:" << std::endl;
        std::cout << c_matrix[0][0] << " " << c_matrix[0][1] << std::endl;
        std::cout << c_matrix[1][0] << " " << c_matrix[1][1] << std::endl;
        
        std::cout << std::endl;
        
        std::cout << ">>> Test chr19_0_relaxed_simple done. <<<" << std::endl;
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
//    Test::test_estimate_initial_probabilities();
//    Test::test_convert_islands_to_string();
//    Test::test_create_emission_to_idx_map();
//    Test::test_forward();
//    Test::test_backward();
//    Test::test_baum_welch_algorithm();
//    Test::test_viterbi();
    //   Test::test_mm39_relaxed_simple_without_baum_welch_algorithm();
    // Test::test_mm39_relaxed_simple();
    Test::test_chr19_0_relaxed_simple();
    return 0;
}