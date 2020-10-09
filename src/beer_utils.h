/**
 * @file beer_utils.h
 *
 * @brief General utility functions used by BEER
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#ifndef BEER_UTILS_H
#define BEER_UTILS_H

#include <vector>
#include <cmath>
#include <string>

#include "Eigen/Eigen"
#include "rapidjson/document.h"
// #include "z3++.h"

/**
 * @brief apply an AND-reduction on a vector of types
 * 
 * @tparam T type of each element
 * @param vec vector of elements to AND-reduce
 * 
 * @return T AND-reduction of vec
 */
template< typename T > 
T z3_and_reduce(const std::vector< T > &vec)
{
	assert(vec.size() > 0);

	T ret = vec.at(0);
	for(std::size_t i = 1; i < vec.size(); i++)
		ret = ret && vec.at(i);
	return ret;
}


/**
 * @brief apply an OR-reduction on a vector of types
 * 
 * @tparam T type of each element
 * @param vec vector of elements to OR-reduce
 * 
 * @return T OR-reduction of vec
 */
template< typename T >
T z3_or_reduce(const std::vector< T > &vec)
{
	assert(vec.size() > 0);

	T ret = vec.at(0);
	for(std::size_t i = 1; i < vec.size(); i++)
		ret = ret || vec.at(i);
	return ret;
}

// compute ECC code UID in the exact same way as EINSim does
uint64_t compute_uid(
      const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &G
    , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &H
    , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &R);

// read an ECC code configuration out of the provided JSON configuration file
rapidjson::Document read_json_cfg_file(const std::string &ecc_code_json_cfg_fname);

// test pattern generator for creating permutations of a base n-hot test pattern
void test_pattern_generator_k_ones_in_n_bits_initialize(std::vector< bool > &test_pattern, uint64_t K, uint64_t k_ones);
bool test_pattern_generator_k_ones_in_n_bits_advance(std::vector< bool > &test_pattern, uint64_t K, uint64_t k_ones);
void tpg_test(void);

#endif /* BEER_UTILS_H */
