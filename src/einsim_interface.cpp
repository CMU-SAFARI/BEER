/**
 * @file einsim_interface.cpp
 *
 * @brief Utility routines that interface with the EINSim simulator and aren't
 * part of the core BEER methodology
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <set>
#include <algorithm>
#include <array>

/* libraries */
#include "Eigen/Eigen"
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/writer.h"
#include "rapidjson/error/en.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"

/* project includes */
#include "beer_utils.h"
#include "einsim_interface.h"

/* forward declaratoins for global variables */
extern int g_verbosity;

/**
 * @brief Subroutine to launch an external binary and return stdout in an std::string
 * 
 * @param cmd The executable command to run
 * @param stdout std::string to return stdout in
 * 
 * @return int return code of the launched process
 */
static int spawn_and_capture_stdout(const char* cmd, std::string &stdout) 
{
    FILE* pipe = popen(cmd, "r");
    if(!pipe)
    {
        std::cerr << "[ERROR] failed to start command \"" << cmd << "\"" << std::endl;
        exit(-1);
    }

    const int BUF_SIZE = 1024; 
    std::string result;
    std::array< char, BUF_SIZE > buffer;
    while(fgets(buffer.data(), BUF_SIZE, pipe) != NULL) 
        stdout.append(buffer.data());
    auto return_code = pclose(pipe);

    return return_code;
}

/**
 * @brief Runs EINSim to generate a Hamming SEC ECC code JSON configuration
 * string corresponding to the desired dataword length and random seed
 *
 * @param einsim_binary Path to the EINSim executable binary
 * @param desired_k Desired data word length
 * @param desired_p Random seed
 * 
 * @return std::string JSON string representing the desired ECC code
 */
std::string einsim_generate_ecc_code(std::string &einsim_binary, uint64_t desired_k, uint64_t desired_p)
{
    std::cout << "[INFO] generating Hamming code of k: " << desired_k << " p: " << desired_p << std::endl;
    
    // prepare EINSim configuration
    // most fields are irrelevant since we are simply harvesting the code itself
    std::stringstream einsim_command;
    einsim_command << einsim_binary << " "
        <<  "-m " << "s"
        << " -n " << "1"
        << " -b " << std::to_string(desired_k)
        << " -k " << std::to_string(desired_k)
        << " -w " << "BLOCKS"
        << " -c " << "ALL_TRUE" // irrelevant
        << " -d " << "RANDOM" // irrelevant
        << " -e " << "DATA_RETENTION,0.5" // irrelevant
        << " -o " << "PER_BIT_ERROR_COUNT"
        << " -s " << "HSC"
        << " -p " << std::to_string(desired_p)
        << " 2>&1";
    if(g_verbosity >= 2)
        std::cout << "[DEBUG2] EINSim command: " << einsim_command.str() << std::endl;

    // run EINSim
    std::string stdout;
    int retcode = spawn_and_capture_stdout(einsim_command.str().c_str(), stdout);
    if(retcode != 0)
    {
        std::cout << "[ERROR] EINSim returned error code " << retcode << std::endl;
        std::cout << "[DEBUG] command:" << einsim_command.str() << std::endl;
        std::cout << "[DEBUG] stdout:" << std::endl;
        std::cout << stdout << std::endl;
        exit(-1);
    }

    // extract useful data from output lines
    std::size_t found = stdout.find("[ECC] ");
    if(found == std::string::npos)
    {
        std::cout << "[ERROR] unable to find ECC scheme in EINSim output" << std::endl;
        std::cout << "[DEBUG] command:" << einsim_command.str() << std::endl;
        std::cout << "[DEBUG] stdout:" << std::endl;
        std::cout << stdout << std::endl;
        exit(-1);
    }

    // find the end
    found = found + 6; // offset the [ECC] 
    std::size_t end = stdout.find('\n', found);
    if(end == std::string::npos)
    {
        std::cout << "[ERROR] unable to find EOL in EINSim output" << std::endl;
        std::cout << "[DEBUG] command:" << einsim_command.str() << std::endl;
        std::cout << "[DEBUG] stdout:" << std::endl;
        std::cout << stdout << std::endl;
        exit(-1);
    }

    std::string raw_json = stdout.substr(found, end - found);

    // insert the miscorrection profile!
    return raw_json;
}

/**
 * @brief Uses EINSim to compute which bit positions can experience errors for a
 * given n-charged test pattern
 *
 * @param einsim_binary Path to executable EINSim binary
 * @param K_gold Number of data bits in the ECC code
 * @param N_gold Number of total (data + parity-check) bits in the ECC code
 * @param ecc_code_json_cfg_fname Path to the JSON configuration file for the ECC code
 * @param noise_ratio Amount of random noise to use when simulating data-retention error injection
 * @param true_anti_cell_layout True-/anti-cell layout
 * @param n_charged_pattern Bit vector showing which bits are charged (1) and discharged (0)
 * @param observed_errors Error counts observed at each bit position
 */
void einsim_compute_observed_errors_for_test_pattern(
      const std::string &einsim_binary
    , const uint64_t K_gold, const uint64_t N_gold 
    , const std::string &ecc_code_json_cfg_fname
    , const float noise_ratio
    , const std::string &true_anti_cell_layout
    , const Eigen::Matrix< bool, 1, Eigen::Dynamic > &n_charged_pattern
    , Eigen::Matrix< long, 1, Eigen::Dynamic > &observed_errors)
{
    // convert the n-charged-pattern to an actual test pattern
    std::string test_pattern_str = "0b";
    for(int i = 0; i < n_charged_pattern.cols(); i++)
    {
        bool test_pattern_bit;
        if(true_anti_cell_layout == "ALL_TRUE")
            test_pattern_bit = n_charged_pattern(i);
        else if(true_anti_cell_layout == "ALL_ANTI")
            test_pattern_bit = !n_charged_pattern(i);
        else if(true_anti_cell_layout == "COLSTRIPE_T")
            test_pattern_bit = n_charged_pattern(i) ? (i % 2 == 0) : (i % 2 == 1);
        else if(true_anti_cell_layout == "COLSTRIPE_A")
            test_pattern_bit = n_charged_pattern(i) ? (i % 2 == 1) : (i % 2 == 0);
        else
        {
            std::cout << "[ERROR] unsupported true-/anti-cell layout: " << true_anti_cell_layout << std::endl;
            assert(0 && "unsupported true-/anti-cell layout");
        }
        test_pattern_str += std::to_string(test_pattern_bit);
    }

    // std::cout << "Test pattern: " << test_pattern_str << std::endl;

    // prepare EINSim configuration
    std::stringstream einsim_command;
    einsim_command << einsim_binary << " "
        <<  "-m " << "s"
        << " -n " << "10000"
        << " -b " << std::to_string(K_gold)
        << " -w " << "BLOCKS"
        << " -c " << true_anti_cell_layout
        << " -d " << test_pattern_str
        << " -e " << (noise_ratio == 0.0 ? "DATA_RETENTION,0.5" : "DATA_RETENTION_NOISY,0.5," + std::to_string(noise_ratio))
        << " -o " << "PER_BIT_ERROR_COUNT"
        << " -s " << ecc_code_json_cfg_fname
        << " 2>&1";
    if(g_verbosity >= 2)
        std::cout << "[DEBUG2] EINSim command: " << einsim_command.str() << std::endl;

    // run EINSim
    std::string stdout;
    int retcode = spawn_and_capture_stdout(einsim_command.str().c_str(), stdout);
    if(retcode != 0)
    {
        std::cout << "[ERROR] EINSim returned error code " << retcode << std::endl;
        std::cout << "[DEBUG] command: " << einsim_command.str() << std::endl;
        std::cout << "[DEBUG] stdout: " << std::endl;
        std::cout << stdout << std::endl;
        exit(-1);
    }

    // std::cout << stdout << std::endl;

    // extract useful data from output lines (output won't be huge)
    std::size_t found;
    std::size_t end = 0;
    std::vector< uint64_t > dataword_counts(K_gold);
    std::vector< uint64_t > codeword_counts(N_gold);
    while((found = stdout.find("[DATA] ", end)) != std::string::npos)
    {
        // find the EOL
        end = stdout.find('\n', found);
        if(end == std::string::npos)
        {
            std::cout << "[ERROR] unable to find EOL in EINSim output" << std::endl;
            std::cout << "[DEBUG] command:" << einsim_command.str() << std::endl;
            std::cout << "[DEBUG] stdout:" << std::endl;
            std::cout << stdout << std::endl;
            exit(-1);
        }

        // extract the error comments
        std::size_t bracket_start = stdout.find('[', found + 7);
        assert(bracket_start != std::string::npos && bracket_start > found && bracket_start < end);
        std::size_t bracket_end = stdout.find(']', bracket_start);
        std::size_t colon_pos = stdout.find(':', bracket_start);
        std::string dataword = stdout.substr(bracket_start + 2, colon_pos - 1 - bracket_start - 2);
        std::string codeword = stdout.substr(colon_pos + 2, bracket_end - 1 - colon_pos - 2);

        // parse out the digits
        uint64_t i;
        std::vector< uint64_t > these_dataword_counts;
        std::vector< uint64_t > these_codeword_counts;
        std::stringstream ssd(dataword);
        std::stringstream ssc(codeword);
        while(ssd >> i)
            these_dataword_counts.push_back(i);
        assert(these_dataword_counts.size() == K_gold);
        while(ssc >> i)
            these_codeword_counts.push_back(i);
        assert(these_codeword_counts.size() == N_gold);

        // add the local values to the global arrays
        for(uint64_t i = 0; i < K_gold; i++) 
            dataword_counts[i] += these_dataword_counts[i];
        for(uint64_t i = 0; i < N_gold; i++) 
            codeword_counts[i] += these_codeword_counts[i];

        // std::cout << stdout.substr(bracket_start, bracket_end - bracket_start + 1) << std::endl;
        // std::cout << "D: ";
        // for(const uint64_t &i : these_dataword_counts) std::cout << i << " ";
        // std::cout << " C: ";
        // for(const uint64_t &i : these_codeword_counts) std::cout << i << " ";
        // std::cout << std::endl;
    }

    // DEBUG: ensure that we captured some errors HOWEVER: note that it is
    // possible for a test pattern with a correctable number of CHARGED cells to
    // show NO errors (i.e., if the parity-check symbols just happen to all be
    // DISCHARGED)
    // uint64_t dataword_count_sum = 0;
    // uint64_t codeword_count_sum = 0;
    // for(uint64_t i = 0; i < K_gold; i++) 
    //     dataword_count_sum += dataword_counts[i];
    // for(uint64_t i = 0; i < N_gold; i++) 
    //     codeword_count_sum += codeword_counts[i];
    // if(dataword_count_sum == 0 || codeword_count_sum == 0)
    // {
    //     std::cout << "[ERROR] no errors generated in EINSim output" << std::endl;
    //     std::cout << "[DEBUG] command:" << einsim_command.str() << std::endl;
    //     std::cout << "[DEBUG] stdout:" << std::endl;
    //     std::cout << stdout << std::endl;
    //     exit(-1);
    // }

    // return the appropriate values
    for(uint64_t i = 0; i < K_gold; i++)
        observed_errors(i) = dataword_counts[i];
}

/**
 * @brief Uses EINSim to compute which bit positions can experience errors for a
 * set of n-charged test patterns
 * 
 * @param einsim_binary Path to executable EINSim binary
 * @param K_gold Number of data bits in the ECC code
 * @param N_gold Number of total (data + parity-check) bits in the ECC code
 * @param ecc_code_json_cfg_fname Path to the JSON configuration file for the ECC code
 * @param noise_ratio Amount of random noise to use when simulating data-retention error injection
 * @param true_anti_cell_layout True-/anti-cell layout
 * @param test_pattern_list List of N-charged test patterns showing which bits are charged (1) and discharged (0)
 * @param observed_errors_per_test_pattern Error counts observed at each bit position for each test pattern
 */
void compute_experimental_observed_errors_using_einsim(
      const std::string &einsim_binary
    , const uint64_t K_gold, const uint64_t N_gold 
    , const std::string &ecc_code_json_cfg_fname
    , const float noise_ratio
    , const std::string &true_anti_cell_layout
    , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &test_pattern_list
    , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &observed_errors_per_test_pattern)
{
    observed_errors_per_test_pattern.resize(test_pattern_list.rows(), K_gold);
    for(Eigen::Index r = 0; r < test_pattern_list.rows(); r++)
    {
        const Eigen::Matrix< bool, 1, Eigen::Dynamic > &tp = test_pattern_list.row(r);
        Eigen::Matrix< long, 1, Eigen::Dynamic > observed_errors_this_tp(K_gold);
        einsim_compute_observed_errors_for_test_pattern(einsim_binary, K_gold, N_gold, ecc_code_json_cfg_fname
            , noise_ratio, true_anti_cell_layout, tp, observed_errors_this_tp);
        observed_errors_per_test_pattern.row(r) = observed_errors_this_tp;
    }
}

/**
 * @brief Computes both (1) list of test patterns and (2) observed errors for
 * each test pattern using EINSim to compute which bit positions can experience
 * errors per test pattern
 *
 * @param einsim_binary Path to executable EINSim binary
 * @param K_gold Number of data bits in the ECC code
 * @param N_gold Number of total (data + parity-check) bits in the ECC code
 * @param ecc_code_json_cfg_fname Path to the JSON configuration file for the ECC code
 * @param true_anti_cell_layout True-/anti-cell layout
 * @param n_hot_patterns Set of n's for the n-hot patterns to generate
 * @param test_pattern_list List of N-charged test patterns showing which bits are charged (1) and discharged (0)
 * @param observed_errors_per_test_pattern Error counts observed at each bit position for each test pattern
 */
void einsim_compute_golden_observed_errors_using_einsim(
      const std::string &einsim_binary
    , const uint64_t K_gold, const uint64_t N_gold 
    , const std::string &ecc_code_json_cfg_fname
    , const std::string &true_anti_cell_layout
    , const std::set< int > &n_hot_patterns
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &test_pattern_list
    , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &observed_errors_per_test_pattern)
{
    std::cout << "[INFO] generating miscorrection profile for code: " << ecc_code_json_cfg_fname << std::endl;
    std::cout << "[INFO] using n_hot_patterns:";
    for(const auto &c : n_hot_patterns) 
        std::cout << " " << c;
    std::cout << std::endl;

    // reset the matrices just in case
    test_pattern_list.resize(0, K_gold);
    observed_errors_per_test_pattern.resize(0, K_gold);

    std::vector< bool > test_pattern;
    for(uint64_t k_ones : n_hot_patterns)
    {
        test_pattern_generator_k_ones_in_n_bits_initialize(test_pattern, K_gold, k_ones);
        do
        {
            Eigen::Matrix< bool, 1, Eigen::Dynamic > test_pattern_mat(K_gold);
            Eigen::Matrix< long, 1, Eigen::Dynamic > miscorrections_mat(K_gold);
            for(uint64_t i = 0; i < K_gold; i++) 
                test_pattern_mat(i) = test_pattern[i];
            einsim_compute_observed_errors_for_test_pattern(einsim_binary, K_gold, N_gold, ecc_code_json_cfg_fname, 0.0, true_anti_cell_layout, test_pattern_mat, miscorrections_mat);

            // add this contribution to the miscorrection profile
            assert((uint64_t)test_pattern_list.cols() == K_gold);
            test_pattern_list.conservativeResize(test_pattern_list.rows() + 1, K_gold);
            test_pattern_list.row(test_pattern_list.rows() - 1) = test_pattern_mat;

            assert((uint64_t)observed_errors_per_test_pattern.cols() == K_gold);
            observed_errors_per_test_pattern.conservativeResize(observed_errors_per_test_pattern.rows() + 1, K_gold);
            observed_errors_per_test_pattern.row(observed_errors_per_test_pattern.rows() - 1) = miscorrections_mat;

            // std::cout << test_pattern_mat << " : " << miscorrections_mat << std::endl;
        } while(test_pattern_generator_k_ones_in_n_bits_advance(test_pattern, K_gold, k_ones));
    }
}

/**
 * @brief Extracts a pre-existing miscorrection profile from a JSON
 * configuration file and also generates entries for missing test patterns based
 * on the requested set of n-hot test patterns
 *
 * @param einsim_binary Path to executable EINSim binary
 * @param json_obj Rapidjson object representing the ECC configuration file
 * @param K_gold Number of data bits in the ECC code
 * @param N_gold Number of total (data + parity-check) bits in the ECC code
 * @param ecc_code_json_cfg_fname Path to the JSON configuration file for the ECC code
 * @param true_anti_cell_layout True-/anti-cell layout
 * @param n_hot_patterns Set of n's for the n-hot patterns to generate
 * @param mc_profile_test_patterns List of N-charged test patterns showing which bits are charged (1) and discharged (0)
 * @param mc_profile_miscorrections Error counts observed at each bit position for each test pattern
 */
void einsim_extract_miscorrection_profile_from_file_and_generate_missing_entries(
      const std::string &einsim_binary
    , const rapidjson::Value &json_obj
    , const uint64_t K_gold
    , const uint64_t N_gold 
    , const std::string &ecc_code_json_cfg_fname
    , const std::string &true_anti_cell_layout
    , const std::set< int > &n_hot_patterns
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
    , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections)
{
    std::cout << "[INFO] extracting miscorrection profile from JSON" << std::endl;
    std::cout << "[INFO] using n_hot_patterns:";
    for(const auto &c : n_hot_patterns) 
        std::cout << " " << c;
    std::cout << std::endl;

    // reset the matrices just in case
    mc_profile_test_patterns.resize(0, K_gold);
    mc_profile_miscorrections.resize(0, K_gold);

    // extract the JSON data
    Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > extracted_mc_test_patterns;
    Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > extracted_mc_miscorrections;
    for(const auto &mc_profile_entry : json_obj.GetArray())
    {
        Eigen::Matrix< bool, 1, Eigen::Dynamic > test_pattern_extracted(K_gold);
        Eigen::Matrix< long, 1, Eigen::Dynamic > miscorrections_extracted(K_gold);
        const rapidjson::Value &json_test_pattern = mc_profile_entry[0];
        const rapidjson::Value &json_miscorrections = mc_profile_entry[1];
        for(rapidjson::SizeType i = 0; i < json_test_pattern.Size(); i++)
        {
            test_pattern_extracted(i) = json_test_pattern[i].GetInt();
            miscorrections_extracted(i) = json_miscorrections[i].GetInt();
        }

        extracted_mc_test_patterns.conservativeResize(extracted_mc_test_patterns.rows() + 1, K_gold);
        extracted_mc_test_patterns.row(extracted_mc_test_patterns.rows() - 1) = test_pattern_extracted;
        assert((uint64_t)extracted_mc_test_patterns.cols() == K_gold);
        
        extracted_mc_miscorrections.conservativeResize(extracted_mc_miscorrections.rows() + 1, K_gold);
        extracted_mc_miscorrections.row(extracted_mc_miscorrections.rows() - 1) = miscorrections_extracted;
        assert((uint64_t)extracted_mc_miscorrections.cols() == K_gold);
    }

    // filter out exactly which ones we need
    std::vector< bool > test_pattern;
    for(uint64_t k_ones : n_hot_patterns)
    {
        test_pattern_generator_k_ones_in_n_bits_initialize(test_pattern, K_gold, k_ones);
        do
        {
            Eigen::Matrix< bool, 1, Eigen::Dynamic > test_pattern_mat(K_gold);
            Eigen::Matrix< long, 1, Eigen::Dynamic > miscorrections_mat(K_gold);
            for(uint64_t i = 0; i < K_gold; i++) 
                test_pattern_mat(i) = test_pattern[i];

            // find this pattern in the JSON object
            Eigen::Index i;
            for(i = 0; i < extracted_mc_test_patterns.rows(); i++)
            {
                if(extracted_mc_test_patterns.row(i) == test_pattern_mat)
                {
                    miscorrections_mat = extracted_mc_miscorrections.row(i);
                    break;
                }
            }
            if(i == extracted_mc_test_patterns.rows()) // nothing found - generate it
                einsim_compute_observed_errors_for_test_pattern(einsim_binary, K_gold, N_gold, ecc_code_json_cfg_fname, 0.0, true_anti_cell_layout, test_pattern_mat, miscorrections_mat);

            assert((uint64_t)mc_profile_test_patterns.cols() == K_gold);
            mc_profile_test_patterns.conservativeResize(mc_profile_test_patterns.rows() + 1, K_gold);
            mc_profile_test_patterns.row(mc_profile_test_patterns.rows() - 1) = test_pattern_mat;

            assert((uint64_t)mc_profile_miscorrections.cols() == K_gold);
            mc_profile_miscorrections.conservativeResize(mc_profile_miscorrections.rows() + 1, K_gold);
            mc_profile_miscorrections.row(mc_profile_miscorrections.rows() - 1) = miscorrections_mat;

            // std::cout << test_pattern_mat << " : " << miscorrections_mat << std::endl;
        } while(test_pattern_generator_k_ones_in_n_bits_advance(test_pattern, K_gold, k_ones));
    }
}
