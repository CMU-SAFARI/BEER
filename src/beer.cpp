/**
 * @file beer.cpp
 *
 * @brief Main entry point and supporting routines for running BEER on arbitrary
 * input data.
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <array>
#include <chrono>
#include <limits>
#include <cctype>
#include <sstream>
#include <set>
#include <algorithm>
#include <random>
#include <fstream>

/* libraries */
#define CXXOPTS_VECTOR_DELIMITER ';' /**< hideous hack to allow backwards compatability */
#include "cxxopts/cxxopts.h"
#include "Eigen/Eigen"
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/writer.h"
#include "rapidjson/error/en.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"
#include "z3++.h"

/* project includes */
#include "unk_hamming_code.h"
#include "einsim_interface.h"
#include "beer_utils.h"

/** tracks the output verbosity, with larger values being more verbose */
int g_verbosity = 0;

/**
 * @brief Generates a miscorrection profile for the supplied code using the SAT
 * solver itself
 *
 * Note that this method is often far slower than using EINSim to numerically
 * brute-force the miscorrection profile because the SAT solver must be run
 * multiple times
 *
 * @param K_gold dataword length of the ECC code to simulate
 * @param N_gold codeword length of the ECC code to simulate 
 * @param ecc_code_json_cfg_fname JSON configuration file for the ECC code to simulate
 * @param n_hot_patterns set of 'n' values for the n-hot test patterns to simulate
 * @param mc_profile_test_patterns generated test patterns
 * @param mc_profile_miscorrections generate miscorrection profile entry
 */
void compute_golden_miscorrection_profile_using_SAT_solver
(
      const uint64_t K_gold, const uint64_t N_gold 
    , const std::string &ecc_code_json_cfg_fname
    , const std::set< int > &n_hot_patterns
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
    , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections
)
{
    std::cout << "[INFO] generating miscorrection profile for code: " << ecc_code_json_cfg_fname << std::endl;
    std::cout << "[INFO] using n_hot_patterns:";
    for(const auto &c : n_hot_patterns) 
        std::cout << " " << c;
    std::cout << std::endl;

    // reset the matrices just in case
    mc_profile_test_patterns.resize(0, K_gold);

    std::vector< bool > test_pattern;
    for(uint64_t k_ones : n_hot_patterns)
    {
        test_pattern_generator_k_ones_in_n_bits_initialize(test_pattern, K_gold, k_ones);
        do
        {
            Eigen::Matrix< bool, 1, Eigen::Dynamic > test_pattern_mat(K_gold);
            for(uint64_t i = 0; i < K_gold; i++) 
                test_pattern_mat(i) = test_pattern[i];

            assert((uint64_t)mc_profile_test_patterns.cols() == K_gold);
            mc_profile_test_patterns.conservativeResize(mc_profile_test_patterns.rows() + 1, K_gold);
            mc_profile_test_patterns.row(mc_profile_test_patterns.rows() - 1) = test_pattern_mat;
        } while(test_pattern_generator_k_ones_in_n_bits_advance(test_pattern, K_gold, k_ones));
    }

    // init a Hamming code and figure out where the miscorrections are
    z3::context z3_ctx;
    z3::solver s(z3_ctx);
    gf2::init_with_z3_context(z3_ctx);
    unk_hamming_code hc_known(z3_ctx, s, ecc_code_json_cfg_fname);
    
    Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > mc_profile_miscorrections_binarized;
    hc_known.compute_miscorrection_profile(mc_profile_test_patterns, mc_profile_miscorrections_binarized, EO_ALL_T);
    mc_profile_miscorrections.resize(mc_profile_miscorrections_binarized.rows(), mc_profile_miscorrections_binarized.cols());
    for(Eigen::Index r = 0; r < mc_profile_miscorrections_binarized.rows(); r++)
        for(Eigen::Index c = 0; c < mc_profile_miscorrections_binarized.cols(); c++)
            mc_profile_miscorrections(r, c) = mc_profile_miscorrections_binarized(r, c);
}


/**
 * @brief Extracts the ECC code parameters out of a JSON configuration file
 *
 * If the JSON configuration file contains the golden miscorrection profile,
 * this function will use the provided miscorrection profile entries and will
 * generate any missing entries using EINSim.
 *
 * @param einsim_binary path to EINSim executable binary
 * @param ecc_code_json_cfg_fname configuration file containing ECC code parameters to read
 * @param true_anti_cell_layout 
 * @param K ECC dataword length
 * @param N ECC codeword length
 * @param NmK number of ECC parity-check bits 
 * @param n_hot_patterns set of the 'n' values for the n-hot patterns to use
 * @param mc_profile_test_patterns read and/or generated test patterns
 * @param mc_profile_miscorrections read and/or generate miscorrection profile entry
 *
 * @return true the JSON configuration file contains the golden miscorrection
 * profile
 * @return false the JSON configuration file does NOT contain the golden
 * miscorrection profile
 */
bool extract_k_n_nmk_and_mcp_from_ecc_code_cfg_file
(
      const std::string &einsim_binary
    , const std::string &ecc_code_json_cfg_fname
    , const std::string &true_anti_cell_layout
    , uint64_t &K, uint64_t &N, uint64_t &NmK
    , const std::set< int > &n_hot_patterns
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
    , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections
)
{
    std::cout << "[INFO] reading Hamming code from file: " << ecc_code_json_cfg_fname << std::endl;

    // read the JSON configuration file in
    rapidjson::Document d = read_json_cfg_file(ecc_code_json_cfg_fname);

    // extract the details we need
    assert(d.HasMember("s") && "malformed JSON file");
    const std::string &ecc_scheme_str = d["s"].GetString();
    if(ecc_scheme_str != "HSC") 
    {
        std::cout << "[ERROR] BEER only supports HSC codes for now. bad configuration file: " << ecc_code_json_cfg_fname << std::endl;
        exit(-1);
    }

    // extract the code details we are interested in
    assert(d.HasMember("k") && "malformed JSON file");
    K = d["k"].GetUint64();
    NmK = unk_hamming_code::compute_hamming_code_n_parity_bits(K);
    N = K + NmK;

    if(d.HasMember("miscorrection_profile"))
    {
        einsim_extract_miscorrection_profile_from_file_and_generate_missing_entries(einsim_binary, d["miscorrection_profile"], K, N
            , ecc_code_json_cfg_fname, true_anti_cell_layout, n_hot_patterns, mc_profile_test_patterns, mc_profile_miscorrections); 
        return true;
    }

    // no profile found
    return false;
}

/**
 * @brief Adds a miscorrection profile to a JSON ECC configuration file in a
 * manner ready for output
 *
 * @param ecc_code_json_cfg_fname configuration file containing ECC code parameters to read
 * @param mc_profile_test_patterns list of test patterns corresponding to miscorrection profile entries
 * @param mc_profile_miscorrections_binarized list of miscorrections possible corresponding to the test patterns
 * 
 * @return std::string JSON configuration for ECC code parameters containing the miscorrection profile
 */
std::string add_miscorrection_profile_to_json_cfg_file
(
      const std::string &ecc_code_json_cfg_fname
    , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
    , const Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections_binarized
)
{
    rapidjson::Document d = read_json_cfg_file(ecc_code_json_cfg_fname);

    // create the array type
    rapidjson::Value mc_profile(rapidjson::kArrayType);
    rapidjson::Document::AllocatorType &allocator = d.GetAllocator();
    assert(mc_profile_test_patterns.rows() == mc_profile_miscorrections_binarized.rows());
    assert(mc_profile_test_patterns.cols() == mc_profile_miscorrections_binarized.cols());
    for(int r = 0; r < mc_profile_test_patterns.rows(); r++)
    {
        rapidjson::Value row_tp(rapidjson::kArrayType);
        rapidjson::Value row_mc(rapidjson::kArrayType);
        for(int c = 0; c < mc_profile_test_patterns.cols(); c++)
        {
            row_tp.PushBack(int(mc_profile_test_patterns(r, c)), allocator);
            row_mc.PushBack(int(mc_profile_miscorrections_binarized(r, c)), allocator);
        }
        rapidjson::Value row(rapidjson::kArrayType);
        row.PushBack(row_tp, allocator);
        row.PushBack(row_mc, allocator);
        mc_profile.PushBack(row, allocator);
    }
    d.AddMember("miscorrection_profile", mc_profile, allocator);

    // no profile
    rapidjson::StringBuffer strbuf;
    rapidjson::PrettyWriter< rapidjson::StringBuffer > writer(strbuf);
    writer.SetFormatOptions(rapidjson::kFormatSingleLineArray);
    d.Accept(writer);
    
    std::string json_str;
    json_str = strbuf.GetString();
    json_str.erase(std::remove_if(json_str.begin(), json_str.end(), ::isspace), json_str.end());

    // REGEX STACK OVERFLOW FOR CODES k120 and greater
    // json_str = std::regex_replace(json_str, std::regex("\\[\\[\\[(.*)]]]"), "[\n          [[$1]]\n    ]"); // format the mc profile
    // json_str = std::regex_replace(json_str, std::regex("], \\["), "]\n        , ["); // format arrays nicely
    // json_str = std::regex_replace(json_str, std::regex("(\\[\\[[01 ,]+])\n        (, \\[[01 ,]+]])"), "$1$2"); // undo line breaks in mc profile
    return json_str;
}

/**
 * @brief Solves for a set of G/H/R matrices that can produce the a given
 * miscorrection profile
 *
 * Returns one possible SAT solution. If is_rerun is set, will incrementally
 * solve for the next set of solutions with the additional constraints of
 * excluding the exclude_Hs set of matrices (along with any of their equivalent codes)
 *
 * @param hc_unk the unknown Hamming code to solve for
 * @param mc_profile_test_patterns list of test patterens in the miscorrection profile
 * @param mc_profile_miscorrections_binarized list of miscorrections in the miscorrection profile
 * @param G solved-for generator matrix
 * @param H solved-for parity-check matrix
 * @param R solved-for degenerator matrix
 * @param exclude_Hs list of H-matrices to exclude in the current solve
 * @param is_rerun whether to incrementally solve over the last invocation
 * @param constrain_mc_exists whether to apply the constraint that a miscorrection must exist 
 * @param constrain_mc_not_exists whether to apply the constraint that a miscorrection does not exist
 * @param eo descriptor for the error model to consider when solving for the ECC function
 * 
 * @return true successsful solve - G/H/R matrices are valid
 * @return false unsuccessful solve - no solution found
 */
bool solve_for_GHR
(
      unk_hamming_code &hc_unk
    , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &mc_profile_test_patterns
    , const Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic > &mc_profile_miscorrections_binarized
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &G
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &H
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &R
    , std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > exclude_Hs
    , const bool is_rerun
    , const bool constrain_mc_exists
    , const bool constrain_mc_not_exists
    , const enum error_operator eo
)
{
    uint64_t K = hc_unk.get_n_data_bits();
    assert((uint64_t)mc_profile_test_patterns.cols() == K);
    assert((uint64_t)mc_profile_miscorrections_binarized.cols() == K);
        
    // establish the profile as constraints and run the SAT solver
    hc_unk.solver_push();
    bool solved = hc_unk.determine_ecc_function_from_binarized_miscorrection_profile(mc_profile_test_patterns
        , mc_profile_miscorrections_binarized, G, H, R, exclude_Hs, is_rerun, constrain_mc_exists, constrain_mc_not_exists, eo);
    hc_unk.solver_pop();

    return solved;
}

/**
 * @brief Primary routine to solve for all possible G/H/R matrices that satisfy
 * the provided miscorrection profile
 *
 * @param K_gold ECC dataword length
 * @param N_gold ECC codeword length
 * @param NmK_gold number of ECC parity-check bits 
 * @param mc_profile_gold_test_patterns list of test patterens in the miscorrection profile
 * @param mc_profile_gold_miscorrections_binarized list of miscorrections in the miscorrection profile
 * @param Gs list of solved-for generator matrices
 * @param Hs list of solved-for parity-check matrices
 * @param Rs list of solved-for degenerator matrices
 * @param constrain_mc_exists whether to apply the constraint that a miscorrection must exist 
 * @param constrain_mc_not_exists whether to apply the constraint that a miscorrection does not exist
 * @param eo descriptor for the error model to consider when solving for the ECC function 
 * @param early_exit whether to immediately exit after the first discovered solution
 * 
 * @return int 
 */
int find_all_solutions_for_mc_profile(
      const uint64_t K_gold
    , const uint64_t N_gold
    , const uint64_t NmK_gold
    , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > mc_profile_gold_test_patterns
    , const Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > mc_profile_gold_miscorrections_binarized
    , std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > &Gs
    , std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > &Hs
    , std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > &Rs
    , const bool constrain_mc_exists
    , const bool constrain_mc_not_exists
    , const enum error_operator eo
    , const bool early_exit
)
{
    // prepare unknown code of GOLDEN LENGTH
    std::cout << "[INFO] preparing unknown Hamming code with " << K_gold << " data bits" << std::endl;
    z3::context z3_ctx;
    z3::solver s(z3_ctx);
    gf2::init_with_z3_context(z3_ctx);
    unk_hamming_code hc_unk(z3_ctx, s, K_gold);
    
    // solve for all possible GHR matrices
    bool is_rerun = false;
    int n_solutions = 0;
    std::cout << "[INFO] solving for GHR matrices using binarized miscorrection profile" << std::endl;
    while(1)
    {
        Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > G;
        Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > H;
        Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > R;
        std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > Hs_excluded;
        if(!Hs.empty())
        {
            const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &latest_base_solution = Hs.back();
            std::vector< int > permutations(NmK_gold);
            for(uint64_t i = 0; i < NmK_gold; i++)
                permutations[i] = i;

            // generate every 'equivalent code' based on permuting the parity bits and row-reducing
            do
            {
                // generate the next 'equivalent code'
                Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > this_sol;
                this_sol.resize(NmK_gold, N_gold);
                this_sol.block(0, 0, NmK_gold, K_gold) = latest_base_solution.block(0, 0, NmK_gold, K_gold); 
                for(uint64_t i = 0; i < NmK_gold; i++)
                    this_sol.col(K_gold + i) = latest_base_solution.col(K_gold + permutations[i]);

                // RREF the matrix
                Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > this_sol_rref = this_sol;
                gf2_rref(this_sol_rref);

                if(g_verbosity > 3)
                {
                    std::cout << "[DEBUG4] Pushing back permutation ";
                    for(const auto &i : permutations)
                        std::cout << i << " ";
                    std::cout << std::endl << this_sol << std::endl;
                    std::cout << std::endl << std::endl << this_sol_rref << std::endl;
                }
                
                // insert this into the discovered solutions
                Hs_excluded.push_back(this_sol_rref);
            } while(std::next_permutation(permutations.begin(), permutations.end()));
        }
        bool success = solve_for_GHR(hc_unk, mc_profile_gold_test_patterns, mc_profile_gold_miscorrections_binarized
            , G, H, R, Hs_excluded, is_rerun, constrain_mc_exists, constrain_mc_not_exists, eo);
        if(!success)
            break;

        // print out the G/H/R matrices that we found this iteration
        if(g_verbosity > 3)
        {
            std::cout << "[DEBUG4]     G:" << std::endl << G << std::endl;
            std::cout << "[DEBUG4]     H:" << std::endl << H << std::endl;
            std::cout << "[DEBUG4]     R:" << std::endl << R << std::endl;
        }

        // add the newly-discovered G/H/R to our running list of solutions
        assert((uint64_t)G.rows() == K_gold && (uint64_t)G.cols() == N_gold);
        assert((uint64_t)H.rows() == NmK_gold && (uint64_t)H.cols() == N_gold);
        assert((uint64_t)R.rows() == K_gold && (uint64_t)R.cols() == N_gold);
        Gs.push_back(G);
        Hs.push_back(H);
        Rs.push_back(R);
        n_solutions++;
        is_rerun = true;

        if(g_verbosity >= 1)
        {
            std::cout << "G:" << std::endl << G << std::endl;
            std::cout << "H:" << std::endl << H << std::endl;
            std::cout << "R:" << std::endl << R << std::endl;\
        }

        // output a JSON version of the discovered result!
        if(g_verbosity >= 2)
        {
            rapidjson::Document d;
            d.SetObject();
            rapidjson::Value G_json(rapidjson::kArrayType);
            rapidjson::Value H_json(rapidjson::kArrayType);
            rapidjson::Value R_json(rapidjson::kArrayType);
            rapidjson::Document::AllocatorType &allocator = d.GetAllocator();
            for(int r = 0; r < G.rows(); r++)
            {
                rapidjson::Value row(rapidjson::kArrayType);
                for(int c = 0; c < G.cols(); c++)
                    row.PushBack(int(G(r, c)), allocator);
                G_json.PushBack(row, allocator);
            }
            for(int r = 0; r < H.rows(); r++)
            {
                rapidjson::Value row(rapidjson::kArrayType);
                for(int c = 0; c < H.cols(); c++)
                    row.PushBack(int(H(r, c)), allocator);
                H_json.PushBack(row, allocator);
            }
            for(int r = 0; r < R.rows(); r++)
            {
                rapidjson::Value row(rapidjson::kArrayType);
                for(int c = 0; c < R.cols(); c++)
                    row.PushBack(int(R(r, c)), allocator);
                R_json.PushBack(row, allocator);
            }
            d.AddMember("GT", G_json, allocator);
            d.AddMember("H", H_json, allocator);
            d.AddMember("R", R_json, allocator);
            
            rapidjson::Value mcp_json(rapidjson::kArrayType);
            assert(mc_profile_gold_test_patterns.rows() == mc_profile_gold_miscorrections_binarized.rows());
            assert(mc_profile_gold_test_patterns.cols() == mc_profile_gold_miscorrections_binarized.cols());
            for(int r = 0; r < mc_profile_gold_test_patterns.rows(); r++)
            {
                rapidjson::Value row_tp(rapidjson::kArrayType);
                rapidjson::Value row_mc(rapidjson::kArrayType);
                for(int c = 0; c < mc_profile_gold_test_patterns.cols(); c++)
                {
                    row_tp.PushBack(int(mc_profile_gold_test_patterns(r, c)), allocator);
                    row_mc.PushBack(int(mc_profile_gold_miscorrections_binarized(r, c)), allocator);
                }
                rapidjson::Value row(rapidjson::kArrayType);
                row.PushBack(row_tp, allocator);
                row.PushBack(row_mc, allocator);
                mcp_json.PushBack(row, allocator);
            }
            d.AddMember("miscorrection_profile", mcp_json, allocator);

            rapidjson::Value k_json;
            k_json.SetUint64(K_gold);
            d.AddMember("k", k_json, allocator);

            rapidjson::Value s_json;
            s_json.SetString("HSC");
            d.AddMember("s", s_json, allocator);

            rapidjson::Value P_json;
            P_json.SetInt64(-1);
            d.AddMember("p", P_json, allocator);

            rapidjson::Value UID_json;
            UID_json.SetUint64(compute_uid(G, H, R));
            d.AddMember("uid", UID_json, allocator);

            rapidjson::StringBuffer strbuf;
            rapidjson::PrettyWriter< rapidjson::StringBuffer > writer(strbuf);
            writer.SetFormatOptions(rapidjson::kFormatSingleLineArray);
            d.Accept(writer);
            
            std::string json_str;
            json_str = strbuf.GetString();
            json_str.erase(std::remove_if(json_str.begin(), json_str.end(), ::isspace), json_str.end());
            std::cout << "[JSON] " << json_str << std::endl;
        }

        std::cout << "[INFO] " << n_solutions << " solutions found so far" << std::endl;

        // if we exit after the first found solution
        if(early_exit)
        {
            std::cout << "[INFO] ceasing search for further solutions due to early exit flag" << std::endl;
            break;
        }
    }
    return n_solutions;
}


/**
 * @brief Applies a binarization to a set of observed errors
 * 
 * the noise filtering algorithm is very simple:
 *     for each test pattern:
 *          apply threshold (based on all non-CHARGED bits) to all non-CHARGED bits
 * 
 * @param test_patterns list of test patterns in the miscorrection profile
 * @param miscorrections list of miscorrection counts in the miscorrection profile
 * @param miscorrections_binarized returned binarized miscorrections
 * @param true_anti_cell_layout true- and anti-cell layout of the codewords
 * @param threshold_ratio_of_the_max [0.0-1.0] ratio of the maximum error count to consider as a "true miscorrection"
 */
void binarize_miscorrections_using_threshold_filter(
      const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &test_patterns
    , const Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &miscorrections
    , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &miscorrections_binarized
    , const std::string &true_anti_cell_layout
    , float threshold_ratio_of_the_max
)
{
    miscorrections_binarized.resize(miscorrections.rows(), miscorrections.cols());

    // set the exact threshold value to work with if a nonzero threshold is requested
    for(Eigen::Index tp_idx = 0; tp_idx < miscorrections.rows(); tp_idx++)
    {
        float threshold_value = 0;
        float meanval = 0;
        float meancount = 0;
        float maxval = std::numeric_limits< float >::min();
        float minval = std::numeric_limits< float >::max();
        if(threshold_ratio_of_the_max > 0)
        {
            for(int c = 0; c < miscorrections.cols(); c++)
            {
                if(test_patterns(tp_idx, c) == false)
                {
                    meanval += miscorrections(tp_idx, c);
                    meancount++;
                    if(miscorrections(tp_idx, c) > maxval)
                        maxval = miscorrections(tp_idx, c);
                    if(miscorrections(tp_idx, c) < minval)
                        minval = miscorrections(tp_idx, c);
                }
            }
            assert(maxval > 0);
            assert(minval <= maxval);
            meanval /= meancount;
            threshold_value = minval + threshold_ratio_of_the_max * (maxval - minval);
        }

        if(g_verbosity >= 2)
            std::cout << "[DEBUG2] test_pattern: [" << test_patterns.row(tp_idx) << "] miscorrections: [" << miscorrections.row(tp_idx) << "] mean: " << meanval << " minval: " << minval
                << " maxval: " << maxval << " threshold_ratio_of_the_max: " << threshold_ratio_of_the_max << " threshold_value: " << threshold_value << std::endl;

        // apply the threshold filter to the miscorrections
        for(int c = 0; c < miscorrections.cols(); c++)
        {
            // bool is_true_cell;
            // if(true_anti_cell_layout == "ALL_TRUE")
            //     is_true_cell = true;
            // else if(true_anti_cell_layout == "ALL_ANTI")
            //     is_true_cell = false;
            // else if(true_anti_cell_layout == "COLSTRIPE_T")
            //     is_true_cell = (c % 2) == 0;
            // else if(true_anti_cell_layout == "COLSTRIPE_A")
            //     is_true_cell = (c % 2) == 1;
            // else
            // {
            //     std::cout << "[ERROR] no matching true_anti_cell_layout: " << true_anti_cell_layout << std::endl;
            //     exit(-1);
            // }
            // bool is_charged = is_true_cell ? test_patterns(tp_idx, c) == true : test_patterns(tp_idx, c) == false;

            if(miscorrections(tp_idx, c) < 0)
                miscorrections_binarized(tp_idx, c) = miscorrections(tp_idx, c); // negative indicates no constraint - pass the negative value along
            // else if(is_charged) // miscorrection unclear - could be retention error
            //     miscorrections_binarized(tp_idx, c) = is_true_cell ? miscorrections(tp_idx, c) > threshold_value : ;
            else // ONLY miscorrections here
                miscorrections_binarized(tp_idx, c) = miscorrections(tp_idx, c) > threshold_value;
        }

    }
}

/**
 * @brief Function to run the BEER analysis on the provided ECC code
 * 
 * @param einsim_binary path to EINSim executable binary
 * @param ecc_code_json_cfg_fname configuration file containing ECC code parameters to read
 * @param n_hot_patterns set of the 'n' values for the n-hot patterns to use
 * @param noise_values 
 * @param true_anti_cell_layout true- and anti-cell layout of the codewords
 * @param threshold_values list of values: [0.0-1.0] ratio of the maximum error count to consider as a "true miscorrection"
 * @param constrain_mc_exists whether to apply the constraint that a miscorrection must exist 
 * @param constrain_mc_not_exists whether to apply the constraint that a miscorrection does not exist
 * @param eo descriptor for the error model to consider when solving for the ECC function
 * @param early_exit whether to immediately exit after the first discovered solution
 */
void run_beer
(
      const std::string &einsim_binary
    , const std::string &ecc_code_json_cfg_fname
    , const std::set< int > &n_hot_patterns
    , const std::set< float > &noise_values
    , const std::string &true_anti_cell_layout
    , const std::set< float > &threshold_values
    , const bool constrain_mc_exists
    , const bool constrain_mc_not_exists
    , const enum error_operator eo
    , const bool early_exit)
{
    uint64_t K_gold, N_gold, NmK_gold;
    Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > mc_profile_gold_test_patterns;
    Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > mc_profile_gold_miscorrections;
    
    // extract code information and get golden MC profile
    bool have_mc_profile_gold = extract_k_n_nmk_and_mcp_from_ecc_code_cfg_file(einsim_binary, ecc_code_json_cfg_fname, true_anti_cell_layout
        , K_gold, N_gold, NmK_gold, n_hot_patterns, mc_profile_gold_test_patterns, mc_profile_gold_miscorrections);
    Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > mc_profile_gold_miscorrections_binarized;
    if(!have_mc_profile_gold)
    {
        // USING EINSim
        auto start_time = std::chrono::high_resolution_clock::now();
        einsim_compute_golden_observed_errors_using_einsim(einsim_binary, K_gold, N_gold, ecc_code_json_cfg_fname, true_anti_cell_layout
            , n_hot_patterns, mc_profile_gold_test_patterns, mc_profile_gold_miscorrections);
        auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
        std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to generate miscorrection profile using EINSim" << std::endl;
        // std::cout << "EINSIM:" << std::endl << mc_profile_gold_miscorrections << std::endl;

        // binarize the miscorrection profile
        binarize_miscorrections_using_threshold_filter(mc_profile_gold_test_patterns, mc_profile_gold_miscorrections
            , mc_profile_gold_miscorrections_binarized, true_anti_cell_layout, 0);

        // // using the SAT solver we get the binarized profile immediately
        // start_time = std::chrono::high_resolution_clock::now();
        // compute_golden_miscorrection_profile_using_SAT_solver(K_gold, N_gold, ecc_code_json_cfg_fname, n_hot_patterns, mc_profile_gold_test_patterns, mc_profile_gold_miscorrections_binarized);
        // time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
        // std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to generate miscorrection profile using SAT solver" << std::endl;
        // std::cout << "SAT:" << std::endl << mc_profile_gold_miscorrections << std::endl;
    }
    else
    {
        // determine if any test patterns are missing and generate them

        // combine them
        mc_profile_gold_miscorrections_binarized.resize(mc_profile_gold_miscorrections.rows(), mc_profile_gold_miscorrections.cols());
        for(Eigen::Index r = 0; r < mc_profile_gold_miscorrections.rows(); r++)
            for(Eigen::Index c = 0; c < mc_profile_gold_miscorrections.cols(); c++)
            {
                long val = mc_profile_gold_miscorrections(r, c);
                assert(mc_profile_gold_miscorrections(r, c) == val && "golden profile read from file should already be binarized");
                mc_profile_gold_miscorrections_binarized(r, c) = val;
            }
    }

    if(g_verbosity >= 1)
    {
        std::cout << "[INFO] binarized golden miscorrection profile:" << std::endl;
        for(int r = 0; r < mc_profile_gold_miscorrections.rows(); r++)
        {
            std::stringstream ss; 
            ss << "[[" << mc_profile_gold_test_patterns.row(r) << "],[" << mc_profile_gold_miscorrections_binarized.row(r) << "]]";
            std::string output = ss.str();
            std::replace(output.begin(), output.end(), ' ', ',');
            std::cout << output << std::endl;
        }
    }

    // emit the JSON for the ECC code
    std::string json_cfg_file = add_miscorrection_profile_to_json_cfg_file(ecc_code_json_cfg_fname
        , mc_profile_gold_test_patterns, mc_profile_gold_miscorrections_binarized);
    std::cout << "[ECC] json cfg string start" << std::endl;
    std::cout << json_cfg_file << std::endl;
    std::cout << "[ECC] json cfg string end" << std::endl;

    // run the BEER analysis
    for(float noise : noise_values)
    {
        // generate an MC profile with noise level specified
        Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > mc_profile_miscorrections;
        if(noise == 0)
            mc_profile_miscorrections = mc_profile_gold_miscorrections;
        else
            compute_experimental_observed_errors_using_einsim(einsim_binary, K_gold, N_gold, ecc_code_json_cfg_fname
                , noise, true_anti_cell_layout, mc_profile_gold_test_patterns, mc_profile_miscorrections);

        if(g_verbosity >= 2)
        {
            std::cout << "[DEBUG2] golden mc profile:" << std::endl << mc_profile_gold_miscorrections << std::endl;
            std::cout << "[DEBUG2] experimental mc profile (noise: " << noise << "):" << std::endl << mc_profile_miscorrections << std::endl;
        }

        // run the solver using a particular threshold filter
        for(float threshold : threshold_values)
        {
            // binarize the miscorrection profile according to the threshold
            Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > mc_profile_miscorrections_binarized;
            binarize_miscorrections_using_threshold_filter(mc_profile_gold_test_patterns, mc_profile_miscorrections
                , mc_profile_miscorrections_binarized, true_anti_cell_layout, threshold);
            
            if(g_verbosity >= 2)
                std::cout << "[DEBUG2] experimental binarized mc profile (threshold: " << threshold 
                    << "):" << std::endl << mc_profile_miscorrections_binarized << std::endl;
    
            // run the solver
            std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > Gs;
            std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > Hs;
            std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > Rs;
            int n_solutions = find_all_solutions_for_mc_profile(K_gold, N_gold, NmK_gold
                , mc_profile_gold_test_patterns, mc_profile_miscorrections_binarized, Gs, Hs, Rs
                , constrain_mc_exists, constrain_mc_not_exists, eo, early_exit);

            if(n_solutions == 1)
                std::cout << "[INFO] SUCCESS! only one solution found :)" << std::endl;
            else if(n_solutions == 0)
                std::cout << "[INFO] FAILURE! no solutions found! system overconstrained - need better filtering" << std::endl;
            else if(n_solutions > 1)
                std::cout << "[INFO] FAILURE! more than one solution found" << std::endl;
            std::cout << "[DATA] noise: " << noise << " threshold: " << threshold << " n_solutions: " << n_solutions << std::endl;
        }
    }
}

/**
 * @brief BEER entry point that parses and handles CLI options
 * 
 * @param argc number of command line arguments
 * @param argv array of command line arguments
 * 
 * @return application return value
 */
int main(int argc, char ** argv) // canNOT const argv due to cxxopts requirements
{
    // save the comand line as a string
    std::stringstream command_line;
    for(int i = 0; i < argc; i++)
        command_line << ' ' << argv[i];
    std::cout << "[CLI]" << command_line.str() << std::endl;

    // parse the CLI options
    cxxopts::Options option_parser("beer", "ECC Reverse-Engineering Tool Suite");
    option_parser.add_options("Common")
        ("einsim_binary", "executable binary for EINSim", cxxopts::value< std::string >())
        ("m, mode", "obtain code via (f)ile or (g)eneration", cxxopts::value< char >()->default_value("f"))
        ("f, code_filename", "ECC code configuration filename", cxxopts::value< std::string >())
        ("k, data_bits", "dataword length of the code to generate", cxxopts::value< uint64_t >())
        ("p, patterns", "n-hot patterns to use (can provide multiple times)", cxxopts::value< std::vector< int > >()->default_value("1"))
        ("r, rseed", "random seed for code generation", cxxopts::value< uint64_t >()->default_value("0"))
        ("l, error_model", "true-/anti-cell layout to simulate", cxxopts::value< std::string >())
        ("n, noise", "noise ratio values to simulate (can provide multiple times)", cxxopts::value< std::vector< float > >()->default_value("0.0"))
        ("t, threshold", "miscorrection profile threshold values to simulate (can provide multiple times)", cxxopts::value< std::vector< float > >()->default_value("0.0"))
        ("c, omit_mc_exists_constraint", "omit the constraint that observed miscorrections must occur")
        ("z, omit_mc_not_exists_constraint", "omit the constraint that unobserved miscorrections must not occur")
        ("e, earlyexit", "terminate after finding first solution and don't check if others exist")
        ("v, verbose", "Print non-essential messages")
        ("h, help", "Show help")
        ;
    option_parser.parse_positional({"einsim_binary"});
    option_parser.positional_help("<string : einsim_binary path>");

    bool needs_help = (argc == 1);
    auto options = option_parser.parse(argc, argv);
    if(needs_help || options.count("help"))
    {
        std::cout << option_parser.help({"", "Common", "positional parameters"}) << std::endl;
        return 0;
    }

    // set g_verbosity
    g_verbosity = options.count("verbose");

    // prepare Eigen library to use multithreading
    Eigen::initParallel();

    // get the einsim binary path
    if(options.count("einsim_binary") == 0)
    {
        std::cout << "[ERROR] must provide EINSim binary path" << std::endl; 
        std::cout << option_parser.help({"", "Common", "positional parameters"}) << std::endl;
        return 0;
    }

    std::string einsim_binary = options["einsim_binary"].as< std::string >();
    std::ifstream f(einsim_binary, std::ifstream::in);
    if(f.good())
    {
        std::cout << "[INFO] einsim binary provided as " << einsim_binary << std::endl;
        f.close();
    }
    else
    {
        std::cout << "[ERROR] invalid einsim binary path: " << einsim_binary << std::endl;
        std::cout << option_parser.help({"", "Simulation"}) << std::endl;
        return -1;            
    }

    std::FILE* ecc_code_json_file = nullptr;
    std::string ecc_code_json_cfg_fname;
    if(options.count("mode") == 0 || std::tolower(options["mode"].as< char >()) == 'f')
    {
        if(options.count("code_filename") == 0)
        {
            std::cout << "[ERROR] ECC code configuration filename required in file mode" << std::endl;
            std::cout << option_parser.help({"", "Common"}) << std::endl;
            return -1;
        }
        ecc_code_json_cfg_fname = options["code_filename"].as< std::string >();
        
        // open the file
        ecc_code_json_file = std::fopen(ecc_code_json_cfg_fname.c_str(), "r");
        if(ecc_code_json_file == nullptr)
        {
            std::cout << "[ERROR] invalid ECC code configuration filename: " << ecc_code_json_cfg_fname << std::endl;
            std::cout << option_parser.help({"", "Simulation"}) << std::endl;
            return -1;            
        }
    }
    else
    {
        if(options.count("data_bits") == 0)
        {
            std::cout << "[ERROR] dataword length required in code generation" << std::endl;
            std::cout << option_parser.help({"", "Common"}) << std::endl;
            return -1;
        }

        uint64_t k = options["data_bits"].as< uint64_t >();
        uint64_t rseed = options["rseed"].as< uint64_t >();
        std::cout << "[RSEED] " << rseed << std::endl; 

        // use EINSim to generate an ECC code and write it to a file
        std::string ecc_code_json_data = einsim_generate_ecc_code(einsim_binary, k, rseed);
        if(g_verbosity > 1)
        {
            std::cout << "[DEBUG] generated ECC code:" << std::endl;
            std::cout << "    " << ecc_code_json_data << std::endl;
        }

        // create a temporary file
        ecc_code_json_cfg_fname = std::string(P_tmpdir) + "/beer_cpp_test";
        ecc_code_json_file = std::fopen(ecc_code_json_cfg_fname.c_str(), "w");
        if(ecc_code_json_file == nullptr)
        {
            std::cout << "[ERROR] invalid ECC code configuration filename: " << ecc_code_json_cfg_fname << std::endl;
            std::cout << option_parser.help({"", "Simulation"}) << std::endl;
            return -1;            
        }

        std::fputs(ecc_code_json_data.c_str(), ecc_code_json_file);
        std::fclose(ecc_code_json_file);
    }

    // n-hot patterns
    const std::vector< int > &n_hot_patterns_list = options["patterns"].as< std::vector< int > >();
    std::set< int > n_hot_patterns(n_hot_patterns_list.begin(), n_hot_patterns_list.end());
    std::cout << "[INFO] simulating " << n_hot_patterns.size() << " n_hot_patterns: ";
    for(const auto& c : n_hot_patterns)
    {
        std::cout << c << " "; 
        assert(c > 0);
    }
    std::cout << std::endl;

    // error operator
    if(options.count("error_model") != 1)
    {
        std::cout << "[ERROR] must provide exactly one error_model" << std::endl;
        return -1;
    }
    const std::string &error_model = options["error_model"].as< std::string >();
    std::cout << "[INFO] simulating error_model: " << error_model << std::endl;
    if(string_to_error_operator.count(error_model) == 0)
    {
        std::cout << "[ERROR] invalid true-/anti-cell layout requested: " << error_model << std::endl;
        return -1;
    }
    enum error_operator eo = string_to_error_operator.at(error_model);
    std::string true_anti_cell_layout;
    if(eo == EO_ALL_T)
        true_anti_cell_layout = "ALL_TRUE";
    else if(eo == EO_ALL_A)
        true_anti_cell_layout = "ALL_ANTI";
    else if(eo == EO_ALT_T
        || eo == EO_DATA_ALT_T_PARITY_ALT_A
        || eo == EO_DATA_ALT_T_PARITY_T
        || eo == EO_DATA_ALT_T_PARITY_A
        || eo == EO_DATA_ALT_T_PARITY_0T
        || eo == EO_DATA_ALT_T_PARITY_1T
        || eo == EO_DATA_ALT_T_PARITY_2T
        || eo == EO_DATA_ALT_T_PARITY_3T
        || eo == EO_DATA_ALT_T_PARITY_4T
        || eo == EO_DATA_ALT_T_PARITY_5T
        || eo == EO_DATA_ALT_T_PARITY_6T
        || eo == EO_DATA_ALT_T_PARITY_7T)
        true_anti_cell_layout = "COLSTRIPE_T";
    else if(eo == EO_ALT_A
        || eo == EO_DATA_ALT_A_PARITY_ALT_T
        || eo == EO_DATA_ALT_A_PARITY_T
        || eo == EO_DATA_ALT_A_PARITY_A
        || eo == EO_DATA_ALT_A_PARITY_0T
        || eo == EO_DATA_ALT_A_PARITY_1T
        || eo == EO_DATA_ALT_A_PARITY_2T
        || eo == EO_DATA_ALT_A_PARITY_3T
        || eo == EO_DATA_ALT_A_PARITY_4T
        || eo == EO_DATA_ALT_A_PARITY_5T
        || eo == EO_DATA_ALT_A_PARITY_6T
        || eo == EO_DATA_ALT_A_PARITY_7T)
        true_anti_cell_layout = "COLSTRIPE_A";


    // noise values
    const std::vector< float > &noise_values_list = options["noise"].as< std::vector< float > >();
    std::set< float > noise_values(noise_values_list.begin(), noise_values_list.end());
    std::cout << "[INFO] simulating " << noise_values.size() << " noise values: ";
    for(const auto& c : noise_values) 
    {
        std::cout << c << " "; 
        assert(c >= 0.0 && c <= 1.0);
    }
    std::cout << std::endl;

    // threshold values 
    const std::vector< float > &threshold_values_list = options["threshold"].as< std::vector< float > >();
    std::set< float > threshold_values(threshold_values_list.begin(), threshold_values_list.end());
    std::cout << "[INFO] simulating " << threshold_values.size() << " threshold values: ";
    for(const auto& c : threshold_values) 
    {
        std::cout << c << " "; 
        assert(c >= 0.0 && c <= 1.0);
    }
    std::cout << std::endl;

    // check which constratins to set
    bool constrain_mc_exists = options.count("omit_mc_exists_constraint") == 0;
    bool constrain_mc_not_exists = options.count("omit_mc_not_exists_constraint") == 0;
    std::cout << "[INFO] using miscorrection constraints for exists/not exists: " 
        << constrain_mc_exists << " / " << constrain_mc_not_exists << std::endl;

    // check if we should exist early
    bool early_exit = options.count("earlyexit") != 0;
    if(early_exit)
        std::cout << "[INFO] will exit after first solution found and will not search for more" << std::endl;

    // we have the code- jump to main routine
    auto start_time = std::chrono::high_resolution_clock::now();
    run_beer(einsim_binary, ecc_code_json_cfg_fname, n_hot_patterns, noise_values, true_anti_cell_layout, threshold_values
        , constrain_mc_exists, constrain_mc_not_exists, eo, early_exit);
    auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to run entire script" << std::endl;
    return 0;
}
