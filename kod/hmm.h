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
            _states(states), _symbols(symbols), state_to_idx(create_state_to_idx_map(states)),
            emission_to_idx(create_emission_to_idx_map(symbols)) {
    }

    void fit(const std::string& states, const std::string& emissions, int n_iter) {
        create_default_emission_to_idx_map();
        create_default_state_to_idx_map();
        estimate_initial_probabilities(states, emissions);
        baum_welch_algorithm(emissions, n_iter);
    }

    std::string predict(const std::string& data) {
        return viterbi_algorithm(data);
    }

   std::pair<double, std::string> evaluate(const std::string& emissions, const std::string& predicted) {
        if (emissions.length() != predicted.length()) {
            std::cerr << "True emissions and predicted emissions differ in length!" << std::endl;
            return {{},{}};
        }

        int correct = 0;
        std::string hit_or_miss;
        for (auto i = 0; i < emissions.length(); ++i) {
            if (emissions[i] == predicted[i]) {
                correct++;
                hit_or_miss.append("+");
            } else {
                hit_or_miss.append("-");
            }
        }

        return {correct / (emissions.length() * 1.0), hit_or_miss};
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
        for (char state: states)
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
                std::string transition{from_state, to_state};
                _transition_probabilities[i][j] = (double) transition_freqs[transition];
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
        for (char _symbol: _symbols)
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

    void create_default_emission_to_idx_map() {
        emission_to_idx['A'] = 0;
        emission_to_idx['C'] = 1;
        emission_to_idx['G'] = 2;
        emission_to_idx['T'] = 3;
    }

    void create_default_state_to_idx_map() {
        state_to_idx['+'] = 0;
        state_to_idx['-'] = 1;
    }

    static std::unordered_map<char, int> create_emission_to_idx_map(const std::vector<char>& emissions) {
        int index = 0;
        std::unordered_map<char, int> emission_to_idx;
        for (char c: emissions)
            emission_to_idx[c] = index++;

        return emission_to_idx;
    }

    
    static std::unordered_map<char, int> create_state_to_idx_map(const std::vector<char>& states) {
        int index = 0;
        std::unordered_map<char, int> state_to_idx;
        for (char c: states)
            state_to_idx[c] = index++;

        //for(char c : set_of_states) 
        //    std::cout << c << " " << state_to_idx[c] << std::endl;

        return state_to_idx;
    }

    double** forward(const std::string& data) {
        int data_length = data.length(); // T
        // data // V
        // _transition_probabilities // a
        // _emission_probabilities // b
        // _states_probabilities // initial_distribution
        double** alpha = 0;
        alpha = new double* [data_length]; // [data_length][N]
        for (int i = 0; i < data_length; i++) {
            alpha[i] = new double[N];
            for (int j = 0; j < N; j++)
                alpha[i][j] = 0;
        }

        int V[data_length];
        for (int i = 0; i < data_length; i++) V[i] = emission_to_idx[data[i]];

        for (int j = 0; j < N; j++) {
            alpha[0][j] = _states_probabilities[j] * _emission_probabilities[j][V[0]];
            // std::cout << _states_probabilities[j] * _emission_probabilities[j][V[0]] << std::endl;
        }

        for (int i = 1; i < data_length; i++) {
            for (int j = 0; j < N; j++) {
                double result = 0;
                for (int k = 0; k < N; k++)
                    result += alpha[i - 1][k] * _transition_probabilities[k][j];

                result *= _emission_probabilities[j][V[i]];
                alpha[i][j] = result;
            }
        }

        return alpha;
    }

    double** backward(const std::string& data) {
        int data_length = data.length(); // T
        // data // V
        // _transition_probabilities // a
        // _emission_probabilities // b

        double** beta = new double* [data_length]; // [data_length][N]
        for (int i = 0; i < data_length; i++) {
            beta[i] = new double[N];
            for (int j = 0; j < N; j++)
                beta[i][j] = 0;
        }

        int V[data_length];
        for (int i = 0; i < data_length; i++) V[i] = emission_to_idx[data[i]];

        for (int j = 0; j < N; j++) {
            beta[data_length - 1][j] = 1;
            // wstd::cout << _states_probabilities[j] * _emission_probabilities[j][V[0]] << std::endl;
        }

        //std::cout << "beta calc" << std::endl;
        for (int i = data_length - 2; i >= 0; i--) {
            for (int j = 0; j < N; j++) {
                double result = 0;
                double* result_matrix = new double[N];
                for (int k = 0; k < N; k++) {
                    result_matrix[k] = beta[i + 1][k] * _emission_probabilities[k][V[i + 1]];
                    //std::cout << beta[i+1][k] << " * " << _emission_probabilities[k][V[i+1]] << " = " << result_matrix[k] << ", ";
                }
                //std::cout << std::endl;

                for (int k = 0; k < N; k++)
                    result += result_matrix[k] * _transition_probabilities[j][k];

                beta[i][j] = result;
                //if(result == 0 && data_length - 1000 < i) std::cout << "i: " << i << ", j:" << j << std::endl;
            }
        }

        return beta;
    }

    std::string viterbi_algorithm(const std::string& data) {
        std::vector<int> DNA;

        std::unordered_map<int, char> idx_to_state;
        for (auto& it: state_to_idx) idx_to_state[it.second] = it.first;

        for (char c: data) DNA.push_back(state_to_idx[c]);

        int data_length = DNA.size();
        int transition_length = sizeof _transition_probabilities / sizeof _transition_probabilities[0];

        double** omega = nullptr;
        omega = new double* [data_length]; // [T][M]
        for (int i = 0; i < data_length; i++) {
            omega[i] = new double[transition_length];
            for (int j = 0; j < transition_length; j++)
                omega[i][j] = 0;
        }

        for (int i = 0; i < transition_length; i++) {
            omega[0][i] = log(_states_probabilities[i] * _emission_probabilities[i][DNA[0]]);
        }

        int** prev = nullptr;
        prev = new int* [data_length - 1]; // [T-1][M]
        for (int i = 0; i < data_length - 1; i++) {
            prev[i] = new int[transition_length];
            for (int j = 0; j < transition_length; j++)
                prev[i][j] = 0;
        }

        for (int t = 1; t < data_length; t++) {
            for (int j = 0; j < transition_length; j++) {
                double max = -INFINITY;
                int max_idx = 0;
                for (int i = 0; i < transition_length; i++) {
                    double probability = omega[t - 1][i] + log(_transition_probabilities[i][j]) +
                                         log(_emission_probabilities[j][DNA[t]]);
                    if (probability > max) {
                        max = probability;
                        max_idx = i;
                    }
                }
                omega[t][j] = max;
                prev[t - 1][j] = max_idx;

            }
        }

        int* S = new int[data_length];
        for (int i = 0; i < data_length; i++) S[i] = 0;

        double max = -INFINITY;
        int max_idx = 0;
        for (int i = 0; i < transition_length; i++) {
            if (omega[data_length - 1][i] > max) {
                max = omega[data_length - 1][i];
                max_idx = i;
            }
        }
        int last_state = max_idx;

        S[0] = last_state;

        int backtrack_index = 1;
        for (int i = data_length - 2; i >= 0; i--) {
            S[backtrack_index] = prev[i][last_state];
            last_state = prev[i][last_state];
            backtrack_index++;
        }

        std::string result;
        for (int i = data_length - 1; i >= 0; i--) {
            result += idx_to_state[S[i]];
        }

        return result;

    }

    void print_transition_probabilities_and__emission_probabilities() {
        std::cout << "_transition_probabilities" << std::endl;
        for (int t_alpha = 0; t_alpha < N; t_alpha++) {
            for (int t_transition = 0; t_transition < N; t_transition++)
                std::cout << _transition_probabilities[t_alpha][t_transition] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << "_transition_probabilities end" << std::endl << std::endl;

        std::cout << "_emission_probabilities" << std::endl;
        for (int vpom1 = 0; vpom1 < N; vpom1++) {
            for (int vpom2 = 0; vpom2 < M; vpom2++)
                std::cout << _emission_probabilities[vpom1][vpom2] << " ";
            std::cout << std::endl;
        }
        std::cout << "_emission_probabilities end" << std::endl << std::endl;
    }

    void baum_welch_algorithm(const std::string& data, int n_iter) {
        int data_length = data.length(); // T
        // data // V
        // _transition_probabilities // a
        // _emission_probabilities // b
        // _states_probabilities // initial_distribution
        // M na internetu je ovdje N
        // print_transition_probabilities_and__emission_probabilities();

        int *V = new int[data_length];
        for (int i = 0; i < data_length; i++) V[i] = emission_to_idx[data[i]];

        for (int n = 0; n < n_iter; n++) {
            double** alpha = forward(data);
            double** beta = backward(data);

            // Inicijalizacija xi
            double*** xi = new double** [N];
            for (int j = 0; j < N; j++) {
                xi[j] = new double* [N];
                for (int k = 0; k < N; k++) {
                    xi[j][k] = new double[data_length - 1];
                    for (int m = 0; m < data_length - 1; m++) xi[j][k][m] = 0;
                }
            }

            for (int t = 0; t < data_length - 1; t++) {
                double* first_dot = new double[N]; // np.dot(alpha2[t, :].T, a2)
                for (int t_alpha = 0; t_alpha < N; t_alpha++) {
                    double first_dot_r = 0;
                    for (int t_transition = 0; t_transition < N; t_transition++)
                        first_dot_r += alpha[t][t_transition] * _transition_probabilities[t_transition][t_alpha];
                    first_dot[t_alpha] = first_dot_r;
                }

                double* second_result = new double[N]; // b2[:, V2[t + 1]].T
                for (int vpom = 0; vpom < N; vpom++)
                    second_result[vpom] = _emission_probabilities[vpom][V[t + 1]];

                double* third_result = new double[N]; // np.dot(alpha2[t, :].T, a2) * b2[:, V2[t + 1]].T
                for (int vpom = 0; vpom < N; vpom++)
                    third_result[vpom] = first_dot[vpom] * second_result[vpom];

                double* fourth_result = new double[N]; // beta2[t + 1, :]
                for (int vpom = 0; vpom < N; vpom++)
                    fourth_result[vpom] = beta[t + 1][vpom];

                // for (int vpom = 0; vpom < N; vpom++) {
                //     if(t == 0) {
                //         std::cout << "third_result[vpom]: " << third_result[vpom] << std::endl;
                //         std::cout << "fourth_result[vpom]: " << fourth_result[vpom] << std::endl;
                //     }
                // }


                double denominator = 0;
                for (int vpom = 0; vpom < N; vpom++)
                    denominator += third_result[vpom] * fourth_result[vpom];
                // if(t == 0)std::cout << "denominator: " << denominator << std::endl;


                for (int i = 0; i < N; i++) {
                    double* numerator = new double[N]; // alpha2[t, i] * a2[i, :] * b2[:, V2[t + 1]].T * beta2[t + 1, :].T
                    for (int vpom = 0; vpom < N; vpom++)
                        numerator[vpom] = alpha[t][i] * _transition_probabilities[i][vpom] *
                                          _emission_probabilities[vpom][V[t + 1]] * beta[t + 1][vpom];

                    for (int vpom = 0; vpom < N; vpom++) 
                        xi[i][vpom][t] = numerator[vpom] / denominator;

                }
            }

            // std::cout << "1" << std::endl;
            // print_transition_probabilities_and__emission_probabilities();


            double** gamma = new double* [N]; // np.sum(xi, axis=1)
            double* gamma_sum_by_1_axis = new double[N]; // np.sum(gamma, axis=1).reshape((-1, 1))
            for (int vpom1 = 0; vpom1 < N; vpom1++) {
                gamma[vpom1] = new double[data_length - 1];
                gamma_sum_by_1_axis[vpom1] = 0;
                for (int vpom2 = 0; vpom2 < data_length - 1; vpom2++) {
                    gamma[vpom1][vpom2] = 0;
                    for (int vpom3 = 0; vpom3 < N; vpom3++)
                        gamma[vpom1][vpom2] += xi[vpom1][vpom3][vpom2];
                    gamma_sum_by_1_axis[vpom1] += gamma[vpom1][vpom2];
                    // std::cout << gamma_sum_by_1_axis[vpom1] << " ";
                }
                // std::cout << std::endl;
            }

            double** pom_var_gamma = new double* [N];// np.sum(xi, 2)
            for (int vpom1 = 0; vpom1 < N; vpom1++) {
                pom_var_gamma[vpom1] = new double[N];
                for (int vpom2 = 0; vpom2 < N; vpom2++) {
                    pom_var_gamma[vpom1][vpom2] = 0;
                    for (int vpom3 = 0; vpom3 < data_length - 1; vpom3++)
                        pom_var_gamma[vpom1][vpom2] += xi[vpom1][vpom2][vpom3];

                    _transition_probabilities[vpom1][vpom2] = pom_var_gamma[vpom1][vpom2] / gamma_sum_by_1_axis[vpom1];
                }
            }

            // std::cout << "2" << std::endl;
            // print_transition_probabilities_and__emission_probabilities();


            // Add additional T'th element in gamma
            double* xi_sum_by_axis_2 = new double[N]; // np.sum(xi[:, :, T - 2], axis=0)
            for (int vpom1 = 0; vpom1 < N; vpom1++) {
                xi_sum_by_axis_2[vpom1] = 0;
                for (int vpom2 = 0; vpom2 < N; vpom2++)
                    xi_sum_by_axis_2[vpom1] += xi[vpom2][vpom1][data_length - 2];
            }

            double** new_gamma = new double* [N];// gamma = np.hstack((gamma, np.sum(xi[:, :, T - 2], axis=0).reshape((-1, 1))))

            for (int vpom1 = 0; vpom1 < N; vpom1++) {
                new_gamma[vpom1] = new double[data_length];
                for (int vpom2 = 0; vpom2 < data_length - 1; vpom2++)
                    new_gamma[vpom1][vpom2] = gamma[vpom1][vpom2];

                new_gamma[vpom1][data_length - 1] = xi_sum_by_axis_2[vpom1];
            }

            // M je K
            double* denominator_2 = new double[N];
            for (int vpom1 = 0; vpom1 < N; vpom1++)
                denominator_2[vpom1] = gamma_sum_by_1_axis[vpom1] + xi_sum_by_axis_2[vpom1];

            for (int vpom1 = 0; vpom1 < M; vpom1++) {
                for (int vpom2 = 0; vpom2 < N; vpom2++) {
                    double pom_new_gamma_sum = 0;
                    for (int vpom3 = 0; vpom3 < data_length; vpom3++) {
                        if (vpom1 == V[vpom3])
                            pom_new_gamma_sum += new_gamma[vpom2][vpom3];
                    }
                    _emission_probabilities[vpom2][vpom1] = pom_new_gamma_sum / denominator_2[vpom2];
                }
            }
        }
        
        // print_transition_probabilities_and__emission_probabilities();
        free(V);

    }
};