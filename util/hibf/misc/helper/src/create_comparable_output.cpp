#include <cassert>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <ranges>
#include <string>
#include <vector>

int const kmer_size{32};
int const number_of_errors{2};

int main(int argc, char ** argv)
{
    if (argc != 4)
        throw std::runtime_error{"Please provide query.names user_bin.ids and mantis.results"};

    std::ifstream query_names_file{argv[1]};
    std::ifstream user_bin_ids_file{argv[2]};

    // read query names
    std::cout << "Reading in " << argv[1] << "... ";
    std::vector<std::string> query_names;
    std::string line;
    while (std::getline(query_names_file, line))
        query_names.push_back(line);
    std::cout << "Done " << std::endl;

    std::cout << "Reading in " << argv[2] << "... ";
    std::unordered_map<std::string, uint64_t> user_bin_ids;
    while (std::getline(user_bin_ids_file, line))
    {
        std::string value{line.begin(), line.begin() + line.find('\t')};
        std::string key{line.begin() + line.find('\t') + 1, line.end()};
        uint64_t value_as_number = std::atoi(value.data());
        user_bin_ids.emplace(key, value_as_number);
    }
    std::cout << "Done " << std::endl;

    std::ifstream mantis_result_in{argv[3]};
    std::ofstream mantis_result_out{"mantis.ready"};

    std::cout << "Processing " << argv[3] << "... ";
    size_t i{0};
    int mantis_threshold{0}; // set when reading seq
    std::vector<uint64_t> results;

    // process first line
    std::getline(mantis_result_in, line);
    mantis_result_out << query_names[i] << '\t';
    std::string kmer_count_str{line.begin() + line.find('\t') + 1, line.end()};
    int kmer_count = std::atoi(kmer_count_str.data());
    mantis_threshold = kmer_count - (kmer_size * number_of_errors);
    ++i;

    while (std::getline(mantis_result_in, line))
    {
        if (line.starts_with("seq"))
        { // new query
            std::sort(results.begin(), results.end());
            if (!results.empty())
                mantis_result_out << results.front();
            for (unsigned r = 1; r < results.size(); ++r)
                mantis_result_out << ',' << results[r];
            mantis_result_out << '\n';
            results.clear();

            assert(i < query_names.size());
            mantis_result_out << query_names[i] << '\t';

            std::string kmer_count_str{line.begin() + line.find('\t') + 1, line.end()};
            int kmer_count = std::atoi(kmer_count_str.data());
            mantis_threshold = kmer_count - (kmer_size * number_of_errors) - 1;

            ++i;
        }
        else
        {
            std::string name{line.begin() + line.find_last_of('/') + 1, line.begin() + line.find_last_of('.')};
            std::string kmer_count_str{line.begin() + line.find('\t') + 1, line.end()};

            int kmer_count = std::atoi(kmer_count_str.data());

            if (kmer_count >= mantis_threshold)
                results.push_back(user_bin_ids[name]);
        }
    }
    // write out last result
    std::sort(results.begin(), results.end());
    if (!results.empty())
        mantis_result_out << results.front();
    for (unsigned r = 1; r < results.size(); ++r)
        mantis_result_out << "," << results[r];
    mantis_result_out << '\n';

    std::cout << "Done " << std::endl;
}
