#include <iostream>
#include <cmath>

#include "hmm.h"
#include "data_loading.h"

class Test {

public:
    static void test_forward() {
        std::cout << ">>> Test forward start. <<<" << std::endl;
        const std::string states = "FBBFFBBFBFFBFFB";
        const std::string emissions = "562364166532654";

        const int N = 2;
        const int M = 6;

        hidden_markov_chain<N, M> hmm({'F', 'B'}, {'1', '2', '3', '4', '5', '6'});
        hmm.estimate_initial_probabilities(states, emissions);

        double** alpha = hmm.forward(emissions);

        std::cout << "alpha" << std::endl;
        for(int i = 0; i < emissions.length(); i++) {
            for(int j = 0; j < N; j++) {
                std::cout << alpha[i][j] << " ";
            }
            std::cout << std::endl;
        }

        delete[] alpha;
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

        double** beta = hmm.backward(emissions);

        std::cout << "beta" << std::endl;
        for(int i = 0; i < emissions.length(); i++) {
            for(int j = 0; j < N; j++) {
                std::cout << beta[i][j] << " ";
            }
            std::cout << std::endl;
        }

        delete[] beta;
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

    static double test_simple(const std::string& train_seq_path,const std::string& train_path,
                const std::string& test_seq_path, const std::string& test_path) {
        std::cout << ">>> Test simple start. <<<" << std::endl;

        std::string train_sequence, train_emissions, test_sequence, test_emissions;
        std::tie(train_sequence, train_emissions) = parse_data(train_path, train_seq_path, true);
        std::tie(test_sequence, test_emissions) = parse_data(test_path, test_seq_path, true);

        hidden_markov_chain<2, 4> hmm({'+', '-'}, {'A', 'C', 'G', 'T'});
        hmm.fit(train_sequence, train_emissions, 1);
        hmm.print_transition_probabilities_and_emission_probabilities();

        auto predicted_emissions = hmm.predict(test_emissions);

        std::string hit_or_miss;
        double accuracy;
        std::tie(accuracy, hit_or_miss) = evaluate(test_sequence, predicted_emissions);
        std::cout << "Accuracy: " << accuracy << std::endl;

        std::cout << predicted_emissions << std::endl;

        std::cout << from_str_to_islands(predicted_emissions) << std::endl;

        std::cout << ">>> Test simple done. <<<" << std::endl;

        return accuracy;
    }

    static void test_complex(const std::string& train_seq_path,const std::string& train_path,
                             const std::string& test_seq_path, const std::string& test_path) {
        std::cout << ">>> Test complex start. <<<" << std::endl;

        std::string train_sequence, train_emissions, test_sequence, test_emissions;
        std::tie(train_sequence, train_emissions) = parse_data(train_path, train_seq_path, false);
        std::tie(test_sequence, test_emissions) = parse_data(test_path, test_seq_path, false);

        hidden_markov_chain<8, 4> hmm({'a', 'c', 'g', 't', 'A', 'C', 'G', 'T'}, {'A', 'C', 'G', 'T'});
        hmm.fit(train_sequence, train_emissions, 1);
        hmm.print_transition_probabilities_and_emission_probabilities();

        auto predicted_emissions = hmm.predict(test_emissions);

        std::string hit_or_miss;
        double accuracy;
        std::tie(accuracy, hit_or_miss) = evaluate(test_sequence, predicted_emissions);
        std::cout << "Accuracy: " << accuracy << std::endl;

        std::cout << predicted_emissions << std::endl;

        std::cout << from_str_to_islands(predicted_emissions, false) << std::endl;

        std::cout << ">>> Test complex done. <<<" << std::endl;
    }

    static void convert_from_str_to_islands(const std::string& island_path, bool simple=true) {
        auto sequence = load_raw(island_path);

        auto islands = from_str_to_islands(sequence, simple);
        std::cout << islands << std::endl;
    }

    static void evaluate_islands(const std::string& seq_path, const std::string& true_path,
                                 const std::string& predict_path, bool simple=true) {
        // true observations are always saved as indexes, don't tell me twice it's a bad design
        std::string sequence, true_obs;
        std::tie(true_obs, sequence) = parse_data(true_path, seq_path);

        std::string predicted_obs = load_raw(predict_path);

        std::string hit_or_miss;
        double accuracy;
        std::tie(accuracy, hit_or_miss) = evaluate(true_obs, predicted_obs, simple);

        std::cout << "Accuracy: " << accuracy << std::endl;
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
            if (!is_equal_arrays<M>(a[i], b[i])) return false;
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
//    Test::test_mm39_relaxed_simple();
//    Test::evaluate_islands(R"(..\data\seven-thrice-nonislands\sequences\chr19_14.txt)", R"(..\data\seven-thrice-nonislands\islands\chr19_14.txt)",
//           R"(..\data\seven-thrice-nonislands\predicted\chr19_14_viterbi.txt)", true);
//    Test::convert_from_str_to_islands(R"(..\data\seven-thrice-nonislands\predicted\chr19_14_viterbi.txt)", true);
    Test::test_simple(R"(..\data\seven-even-islands\sequences\chr19_0.txt)", R"(..\data\seven-even-islands\islands\chr19_0.txt)",
      R"(..\data\seven-thrice-nonislands\sequences\chr19_14.txt)", R"(..\data\seven-thrice-nonislands\islands\chr19_14.txt)");
//    Test::test_complex(R"(..\data\even\sequences\chr19_0.txt)", R"(..\data\even\islands\chr19_0.txt)",
//    R"(..\data\even\sequences\chr19_3.txt)", R"(..\data\even\islands\chr19_3.txt)");
    return 0;
}