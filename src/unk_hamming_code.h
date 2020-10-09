/**
 * @file unk_hamming_code.h
 *
 * @brief Representation of an unknown Hamming code and associated Z3 processing routines
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#ifndef UNK_HAMMING_CODE
#define UNK_HAMMING_CODE

#include <cmath>
#include <chrono>

#include "Eigen/Eigen"
#include "z3++.h"

#include "gf2.h"
#include "beer_utils.h"

// basic operations on GF(2) types that are wrapped in Eigen objects
Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > z3_eval_Eigen_type(const z3::model &model, const Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > &m);
Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > gf2_to_int(const Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > &m);
Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > gf2_to_bool(const Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > &m);
void gf2_rref(Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &M, bool left_rref=false);

/**
 * @brief Enum representation of different supported error-injection models
 */
enum error_operator
{
      EO_ALL_T
    , EO_ALL_A
    , EO_ALT_T
    , EO_DATA_ALT_T_PARITY_ALT_A
    , EO_DATA_ALT_T_PARITY_T
    , EO_DATA_ALT_T_PARITY_A
    , EO_ALT_A
    , EO_DATA_ALT_A_PARITY_ALT_T
    , EO_DATA_ALT_A_PARITY_T
    , EO_DATA_ALT_A_PARITY_A

    // brute-force
    , EO_DATA_ALT_T_PARITY_0T
    , EO_DATA_ALT_T_PARITY_1T
    , EO_DATA_ALT_T_PARITY_2T
    , EO_DATA_ALT_T_PARITY_3T
    , EO_DATA_ALT_T_PARITY_4T
    , EO_DATA_ALT_T_PARITY_5T
    , EO_DATA_ALT_T_PARITY_6T
    , EO_DATA_ALT_T_PARITY_7T

    , EO_DATA_ALT_A_PARITY_0T
    , EO_DATA_ALT_A_PARITY_1T
    , EO_DATA_ALT_A_PARITY_2T
    , EO_DATA_ALT_A_PARITY_3T
    , EO_DATA_ALT_A_PARITY_4T
    , EO_DATA_ALT_A_PARITY_5T
    , EO_DATA_ALT_A_PARITY_6T
    , EO_DATA_ALT_A_PARITY_7T
};

// convertion between enum and string representations
extern std::map< enum error_operator, std::string > error_operator_to_string;
extern std::map< std::string, enum error_operator > string_to_error_operator;

/**
 * @brief representation of an unknown Hamming code in terms of z3 objects
 * 
 */
class unk_hamming_code
{
private:
    uint64_t N, K, NmK;
    z3::context &ctx; // z3 solver
    z3::solver &sol; // z3 solver

    Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > G;
    Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > H;
    Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > R;
    Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > GT;
    Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > HT;
    Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > RT;

    // internal initializers for the unknown Hamming code
    void build_ecc_code_systematic_and_standard_form(void);
    void extract_ecc_code_from_cfg_file(const std::string &json_ecc_code_cfg_file);

    // internal worker to determine the ECC function using a given miscorrection profile
    void determine_ecc_function_from_binarized_miscorrection_profile_impl(
          const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
        , const Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections_binarized
        , const bool constrain_mc_exists
        , const bool constrain_mc_not_exists
        , const enum error_operator eo
    );

public:
    /**
     * @brief Construct a new unk hamming code object given a JSON configuration file
     * 
     * @param z3_ctx a Z3 context
     * @param z3_solver a Z3 solver instance
     * @param json_ecc_code_cfg_file JSON configuration file that specifies an ECC code
     */
    unk_hamming_code(z3::context &z3_ctx, z3::solver &z3_solver, const std::string &json_ecc_code_cfg_file) 
        : ctx(z3_ctx), sol(z3_solver)
    { 
        auto start_time = std::chrono::high_resolution_clock::now();
        extract_ecc_code_from_cfg_file(json_ecc_code_cfg_file);
        auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
        std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to extract hamming_code from cfg file (k: " << K << " n: " << N << ")" << std::endl;
    }

    /**
     * @brief Construct a new unk hamming code object given a JSON configuration file
     * 
     * @param z3_ctx a Z3 context
     * @param z3_solver a Z3 solver instance
     * @param k the ECC dataword length for the ECC code to generate
     */
    unk_hamming_code(z3::context &z3_ctx, z3::solver &z3_solver, uint64_t k) 
        : ctx(z3_ctx), sol(z3_solver)
    {
        this->K = k;
        this->NmK = compute_hamming_code_n_parity_bits(K);
        this->N = NmK + K;

        auto start_time = std::chrono::high_resolution_clock::now();
        build_ecc_code_systematic_and_standard_form();
        auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
        std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to build unknown hamming_code (k: " << K << " n: " << N << ")" << std::endl;
    }

    uint64_t get_n_code_bits(void) const { return N; } /**< get the ECC codeword length */
    uint64_t get_n_data_bits(void) const { return K; } /**< get the ECC dataword legnth */
    uint64_t get_n_parity_bits(void) const { return NmK; } /**< get the number of ECC parity-check bits */

    void solver_push(void) { sol.push(); } /**< save the current state of the z3 solver */
    void solver_pop(void) { sol.pop(); } /**< restore the state of the z3 solver */

    void apply_error_operator_to_codeword(    
          Eigen::Matrix< gf2, 1, Eigen::Dynamic > &codeword_p
        , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &codeword
        , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
        , const enum error_operator eo);

    Eigen::Matrix< gf2, 1, Eigen::Dynamic > calculate_syndrome_from_codeword(
          const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &codeword
        , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
        , const enum error_operator eo
    );

    Eigen::Matrix< gf2, 1, Eigen::Dynamic > calculate_syndrome_from_dataword(
          const Eigen::Matrix< bool, 1, Eigen::Dynamic > &dataword
        , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
        , const enum error_operator eo);

    void calculate_everything(
          const Eigen::Matrix< bool, 1, Eigen::Dynamic > &dataword
        , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
        , const enum error_operator eo
        , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &codeword
        , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &codeword_p
        , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &corrected_codeword
        , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &dataword_p
        , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &syndrome
    );

    z3::expr enforce_error_mask_causes_syndrome_idx(
          const Eigen::Matrix< bool, 1, Eigen::Dynamic > &dataword
        , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
        , uint64_t syndrome_idx
        , const enum error_operator eo);

    void compute_miscorrection_profile(
          const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
        , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections_binarized
        , const enum error_operator eo);

    bool determine_ecc_function_from_binarized_miscorrection_profile(
          const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
        , const Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections_binarized
        , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &G_ret
        , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &H_ret
        , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &R_ret
        , std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > Hs
        , const bool is_rerun
        , const bool constrain_mc_exists
        , const bool constrain_mc_not_exists
        , const enum error_operator eo
        );

    /**
     * @brief utility function to compute the number of parity-check bits used
     * by a Hamming code with k parity-check bits
     *
     * @param k ECC dataword length
     * @return uint64_t number of ECC parity-check bits
     */
    static uint64_t compute_hamming_code_n_parity_bits(uint64_t k)
    {
        uint64_t np = 0;
        while((1u << np) < (np + k + 1u))
            np += 1;
        assert(np == std::ceil(std::log2(k + np + 1)) && "MISMATCH - INCORRECT LOGIC");
        return np;
    }
};

#endif /* UNK_HAMMING_CODE */
