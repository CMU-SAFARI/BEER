/**
 * @file beer_utils.cpp
 *
 * @brief General utility functions used by BEER.
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <algorithm>

/* libraries */
#include "Eigen/Eigen"
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/writer.h"
#include "rapidjson/error/en.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"
#include "crc64/crc64.h"

/* project includes */
#include "beer_utils.h"

/* forward declaratoins for global variables */
extern int g_verbosity;

/**
 * @brief Routine for extracting the contents of a JSON file as a rapidjson::Document
 * 
 * @param ecc_code_json_cfg_fname filename of the JSON configuration file
 * @return rapidjson::Document rapidjson document representation
 */
rapidjson::Document read_json_cfg_file(const std::string &ecc_code_json_cfg_fname)
{
    std::ifstream ifs(ecc_code_json_cfg_fname);
    rapidjson::IStreamWrapper isw(ifs);

    rapidjson::Document d;
    rapidjson::ParseResult result = d.ParseStream< rapidjson::kParseCommentsFlag >(isw);
    if(!result) 
    {
        std::cout << "[ERROR] JSON parse error while parsing ECC configuration file " <<
            ecc_code_json_cfg_fname << ": " << rapidjson::GetParseError_En(result.Code()) 
            << " (" << result.Offset() << ")" << std::endl;
        exit(-1);
    }
    ifs.close();

    return d;
}

/**
 * @brief hash fucntion used to compute the UID of an ECC code (copied from EINSim)
 * 
 * @tparam T matrix type
 * @param mat_list list of matrices to compute the UID of
 * @return uint64_t computed UID
 */
template< typename T >
uint64_t hash_matrix(std::initializer_list< T > mat_list) 
{
    std::stringstream ss;
    for(const auto &mat : mat_list)
        for(int r = 0; r < mat.rows(); r++)
            for(int c = 0; c < mat.cols(); c++)
                ss << mat(r, c);
    return crc64(-1, ss.str().c_str(), ss.str().length());
};

/**
 * @brief Compute the UID for an ECC code in the same manner that EINSim would
 * 
 * @param G ECC generator matrix
 * @param H ECC parity-check matrix
 * @param R ECC degenerator matrix
 * @return uint64_t computed UID
 */
uint64_t compute_uid(
      const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &G
    , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &H
    , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &R)
{
    Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > GT = G.transpose();
    return hash_matrix({GT, H, R});
}


/**
 * @brief Initializes the test pattern with the first permutation in the sequence
 * 
 * @param test_pattern Test patern to initialize
 * @param K Desired word length
 * @param k_ones Number of bits that should be set
 */
void test_pattern_generator_k_ones_in_n_bits_initialize(std::vector< bool > &test_pattern, uint64_t K, uint64_t k_ones)
{
    test_pattern.resize(K);
    for(uint64_t i = 0; i < K; i++)
        test_pattern[i] = (i >= K - k_ones) ? 1 : 0;
    std::sort(test_pattern.begin(), test_pattern.end());
}

/**
 * @brief Advance the test pattern to the next sequential permutation (using
 * C++'s std::next_permutation)
 *
 * @param test_pattern Test patern to initialize
 * @param K Desired word length
 * @param k_ones Number of bits that should be set
 *
 * @return true More permutations exist
 * @return false No more permutations exist 
 */
bool test_pattern_generator_k_ones_in_n_bits_advance(std::vector< bool > &test_pattern, uint64_t K, uint64_t k_ones)
{
    return std::next_permutation(test_pattern.begin(), test_pattern.end());
}

/**
 * @brief Simple test routine to sanity-check the test pattern generator
 */
void tpg_test(void)
{
    std::vector< bool > test_pattern;
    uint64_t K = 5;
    for(uint64_t k_ones = 1; k_ones < 4; k_ones++)
    {
        std::cout << "K_ones: " << k_ones << std::endl;
        test_pattern_generator_k_ones_in_n_bits_initialize(test_pattern, K, k_ones);
        do
        {
            std::cout << "    ";
            for(bool bit : test_pattern) std:: cout << bit; 
            std::cout << std::endl;
        } while(test_pattern_generator_k_ones_in_n_bits_advance(test_pattern, K, k_ones));
    }
}