#include <cassert>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <ranges>
#include <string>
#include <vector>

int main(int argc, char ** argv)
{
    constexpr int kmer_size{32};
    constexpr int number_of_errors{2};
    constexpr int destroyed_kmers{kmer_size * number_of_errors};

    // $FILENAME_QUERY_NAMES
    // $FILENAME_USER_BIN_IDS
    // $FILENAME_MANTIS_QUERY_RESULT
    if (argc != 4)
        throw std::runtime_error{"Please provide query.names user_bin.ids and mantis.results"};

    std::string line_buffer; // Buffer for file I/O

    // Parse query names
    std::vector<std::string> query_names;
    {
        std::ifstream query_names_file{argv[1]};
        std::cout << "Reading " << argv[1] << "... ";
        // Contains lines: "query_name"
        while (std::getline(query_names_file, line_buffer))
            query_names.push_back(line_buffer);
        std::cout << "Done " << std::endl;
    }

    // Parse user bin name-to-id info
    std::unordered_map<std::string, uint64_t> ub_name_to_id; // [reference_name] = number
    {
        std::ifstream user_bin_ids_file{argv[2]};
        std::cout << "Reading in " << argv[2] << "... ";
        // Contains lines: "some_number <tab> reference_name"
        while (std::getline(user_bin_ids_file, line_buffer))
        {
            auto tab_it{line_buffer.begin() + line_buffer.find('\t')};
            std::string_view id_value{line_buffer.begin(), tab_it};
            std::string_view name_key{++tab_it, line_buffer.end()};
            ub_name_to_id.emplace(name_key, std::atoi(id_value.data()));
        }
        std::cout << "Done " << std::endl;
    }

    // Process mantis results
    std::ifstream mantis_result_in{argv[3]};
    std::ofstream mantis_result_out{"mantis.ready"};

    std::cout << "Processing " << argv[3] << "... ";
    // ## Threshold:
    // Let:
    //   * p patternsize
    //   * k k-mer size
    //   * e errors
    //   * c k-mer count
    //   * t threshold
    // Then:
    //   * c = p - k + 1 [Lemma A]
    //   * p = c + k - 1 [Eq 1]
    //   * t = p - (e + 1) * k + 1 [Lemma B]
    //   * t = c + k - 1 - (e + 1) * k + 1 [Eq 1 + Lemma B]
    //   * t = c + 1 - 1 + k - k - e * k
    //   * t = c - e * k
    // However, for mantis, we need - 1. Maybe this is some edge case in mantis that isn't covered.
    // 155 is the threshold, 154 gets reported.
    auto threshold = [destroyed_kmers] (int const kmer_count)
    {
        assert(kmer_count > 0);
        assert(kmer_count - 1 >= destroyed_kmers);
        return kmer_count - (destroyed_kmers) - 1;
    };

    int mantis_threshold{}; // Needs to be set for each query.

    // ## Mantis results ##
    // There is no query ID, instead they are enumerated: seq0 - seqX
    // For each read:
    //   * seqX <tab> kmers in query
    //   * For each user bin that has kmers, the count is reported:
    //     * absolute_path/to/ub.squeakr <tab> count
    // E.g.:
    // seq0    219
    // [...]/GCF_016028855.1_ASM1602885v1_genomic.squeakr       173
    // [...]/GCF_016128175.1_ASM1612817v1_genomic.squeakr       205
    // [...]/GCF_020162095.1_ASM2016209v1_genomic.squeakr       1
    // seq1    219
    std::string ub_name_buffer{};
    auto parse_user_bin_id = [&ub_name_buffer, &ub_name_to_id] (std::string const & line)
    {
        ub_name_buffer.assign(line.begin() + line.find_last_of('/') + 1, // Skip absolute path
                              line.begin() + line.find_last_of('.')); // Skip .squeakr extension
        return ub_name_to_id[ub_name_buffer];
    };

    auto parse_kmer_count = [] (std::string const & line)
    {
        std::string_view const sv{line.begin() + line.find('\t') + 1, // Skip seqX
                                  line.end()};
        return std::atoi(sv.data());
    };

    size_t current_query_number{};
    std::vector<uint64_t> results;

    // // process first line
    // std::getline(mantis_result_in, line);
    // mantis_result_out << query_names[i] << '\t';
    // std::string kmer_count_str{line.begin() + line.find('\t') + 1, line.end()};
    // int kmer_count = std::atoi(kmer_count_str.data());
    // mantis_threshold = threshold(kmer_count);
    // ++i;
    // ^^^ old
    // vvv TODO: Why is this a special case? Shouldn't it be covered by the while loop?

    std::string result_buffer{};

    while (std::getline(mantis_result_in, line_buffer))
    {
        if (line_buffer.starts_with("seq")) // new query
        {
            // // Process the results of the previous query.
            // std::sort(results.begin(), results.end());
            // if (!results.empty())
            //     mantis_result_out << results.front();
            // for (unsigned r = 1; r < results.size(); ++r)
            //     mantis_result_out << ',' << results[r];
            // mantis_result_out << '\n';
            // results.clear();
            // ^^^ old
            // vvv TODO: this could span everything regarding results; then the first line doesn't need special treatment?
            // Process the results of the previous query.
            if (!results.empty())
            {
                std::sort(results.begin(), results.end());
                for (size_t const ub : results)
                    result_buffer += std::to_string(ub) + ',';
                result_buffer.back() = '\n';
                mantis_result_out << result_buffer;
                result_buffer.clear();
                results.clear();
            }

            // Output new query name.
            assert(current_query_number < query_names.size());
            mantis_result_out << query_names[current_query_number] << '\t';

            // Compute threshold for current query.
            mantis_threshold = threshold(parse_kmer_count(line_buffer));

            ++current_query_number;
        }
        else
        {
            if (parse_kmer_count(line_buffer) >= mantis_threshold)
                results.push_back(parse_user_bin_id(line_buffer));
        }
    }
    // write out last result
    if (!results.empty())
    {
        std::sort(results.begin(), results.end());
        for (size_t const ub : results)
            result_buffer += std::to_string(ub) + ',';
        result_buffer.back() = '\n';
        mantis_result_out << result_buffer;
    }
    // std::sort(results.begin(), results.end());
    // if (!results.empty())
    //     mantis_result_out << results.front();
    // for (unsigned r = 1; r < results.size(); ++r)
    //     mantis_result_out << "," << results[r];
    // mantis_result_out << '\n';

    std::cout << "Done " << std::endl;
}
