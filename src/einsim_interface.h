/**
 * @file einsim_interface.h
 *
 * @brief Routines that interface BEER with the EINSim simulator
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#ifndef EINSIM_INTERFACE_H
#define EINSIM_INTERFACE_H

#include <string>
#include <set>

/* external libraries */
#include "Eigen/Eigen"

// generate a new ECC code using EINSim
std::string einsim_generate_ecc_code(std::string &einsim_binary, uint64_t desired_k, uint64_t desired_p);

// different methods for generating miscorrections for given test pattern sets
// TODO: fix documentation and discrepancy between 'miscorrections' and 'observed_errors'
void einsim_compute_observed_errors_for_test_pattern(
      const std::string &einsim_binary
    , const uint64_t K_gold, const uint64_t N_gold 
    , const std::string &ecc_code_json_cfg_fname
    , const float noise_ratio
    , const std::string &true_anti_cell_layout
    , const Eigen::Matrix< bool, 1, Eigen::Dynamic > &n_charged_pattern
    , Eigen::Matrix< long, 1, Eigen::Dynamic > &observed_errors);

void compute_experimental_observed_errors_using_einsim(
      const std::string &einsim_binary
    , const uint64_t K_gold, const uint64_t N_gold 
    , const std::string &ecc_code_json_cfg_fname
    , const float noise_ratio
    , const std::string &true_anti_cell_layout
    , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &test_pattern_list
    , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &observed_errors_per_test_pattern);

void einsim_compute_golden_observed_errors_using_einsim(
      const std::string &einsim_binary
    , const uint64_t K_gold, const uint64_t N_gold 
    , const std::string &ecc_code_json_cfg_fname
    , const std::string &true_anti_cell_layout
    , const std::set< int > &n_hot_patterns
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &test_pattern_list
    , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &observed_errors_per_test_pattern);

void einsim_extract_miscorrection_profile_from_file_and_generate_missing_entries(
      const std::string &einsim_binary
    , const rapidjson::Value &json_obj
    , const uint64_t K_gold
    , const uint64_t N_gold 
    , const std::string &ecc_code_json_cfg_fname
    , const std::string &true_anti_cell_layout
    , const std::set< int > &n_hot_patterns
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
    , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections);

#endif /* EINSIM_INTERFACE_H */