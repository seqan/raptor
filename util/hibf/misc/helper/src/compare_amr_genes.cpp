#include <cassert>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/std/algorithm>

int main(int argc, char ** argv)
{
    if (argc != 3)
        throw std::runtime_error{"Please provide a processed raptor result file and the truth file"};

    std::ifstream raptor_result{argv[1]};
    std::ifstream truth_file{argv[2]};

    if (!raptor_result.good())
        throw std::runtime_error{"Could not open file " + std::string{argv[1]}};
    if (!truth_file.good())
        throw std::runtime_error{"Could not open file " + std::string{argv[2]}};

    std::unordered_map<std::string, uint64_t> user_bin_ids;

    std::cout << "Reading in user bin ids from raptor header in " << argv[1] << "... ";
    std::string line;
    while (std::getline(raptor_result, line) && line[0] == '#' && line[1] != 'Q')
    {
        std::string const value{line.begin(), line.begin() + line.find('\t')};
        std::string const key{line.begin() + line.find_last_of('/') + 1, line.begin() + line.find_first_of('.')};
        uint64_t value_as_number = std::atoi(value.data());
        user_bin_ids.emplace(key, value_as_number);
    }
    std::cout << "Done " << std::endl;

    std::unordered_map<std::string, std::vector<uint64_t>> truth_set{};

    std::cout << "Reading in truth set from file " << argv[2] << "... ";
    // process header line
    std::getline(truth_file, line);
    std::vector<std::string> genes{};
    for (auto && gene : line | std::views::split('\t'))
    {
        std::string gene_str = [] (auto v) { std::string s; for (auto c : v) s.push_back(c); return s; }(gene);
        genes.push_back(gene_str);
    }
    // process rest of files
    while (std::getline(truth_file, line))
    {
        std::string sample_id{line.begin(), line.begin() + line.find('\t')};
        uint64_t sample_idx = user_bin_ids[sample_id];

        size_t current_pos{0};
        for (auto && occ : line | std::views::split('\t'))
        {
            if (std::ranges::equal(occ, std::string{"1"})) // not 0
                truth_set[genes[current_pos]].push_back(sample_idx);
            ++current_pos;
        }
    }
    std::cout << "Done - Truth set has size " << truth_set.size() << std::endl;

    std::cout << "Processing Results from raptor file " << argv[1] << "... ";

    std::ofstream false_positives_file{"raptor.fps"};
    std::ofstream false_negatives_file{"raptor.fns"};
    uint64_t true_positives{0};
    uint64_t false_positives{0};
    uint64_t false_negatives{0};
    uint64_t line_no{0};
    uint64_t all_raptor{0};
    size_t skipped_genes{};

    while (std::getline(raptor_result, line))
    {
        auto gv = line | std::views::split('|') | std::views::drop(5);
        std::string gene = [] (auto v) { std::string s; for (auto c : v) s.push_back(c); return s; }(*gv.begin());

        auto it = truth_set.find(gene);
        if (it == truth_set.end())
        {
            ++skipped_genes;
            std::cerr << "Warning: Could not find gene '" << gene << "' in truth set.\n";
            continue;
        }

        auto & truth_fields = it->second;

        std::string raptor_fields{line.begin() + line.find('\t') + 1, line.end()};
        auto raptor_fields_view = raptor_fields | std::views::split(',');

        auto truth_it = truth_fields.begin();
        auto raptor_it = raptor_fields_view.begin();

        while (truth_it != truth_fields.end() && raptor_it != raptor_fields_view.end())
        {
            std::string raptor_str = [] (auto v) { std::string s; for (auto c : v) s.push_back(c); return s; }(*raptor_it);
            uint64_t raptor_value = std::atoi(raptor_str.data());

            uint64_t truth_value = *truth_it;

            if (truth_value != raptor_value) // If mantis results are empty, then...?
            {
                if (truth_value < raptor_value)
                {
                    false_negatives_file << gene << ":" << truth_value << '\n';
                    ++false_negatives;
                    // ++all_mantis;
                    ++truth_it;
                }
                else
                {
                    false_positives_file << gene << ":" << raptor_value << '\n';
                    ++false_positives;
                    ++all_raptor;
                    ++raptor_it;
                }
            }
            else
            {
                ++true_positives;
                //  ++all_mantis;
                 ++all_raptor;
                 ++truth_it;
                 ++raptor_it;
            }
        }

        while (truth_it != truth_fields.end()) // process the rest of mantis
        {
            uint64_t truth_value = *truth_it;
            ++false_negatives;
            false_negatives_file << gene << ":" << truth_value << '\n';
            // ++all_mantis;
            ++truth_it;
        }

        while (raptor_it != raptor_fields_view.end()) // process the rest of raptor if any
        {
            std::string raptor_str = [] (auto v) { std::string s; for (auto c : v) s.push_back(c); return s; }(*raptor_it);
            uint64_t raptor_value = std::atoi(raptor_str.data());
            false_positives_file << gene << ":" << raptor_value << '\n';
            ++false_positives;
            ++all_raptor;
            ++raptor_it;
        }

        ++line_no;
    }

    std::cout << std::endl;
    // std::cout << "Mantis total #hits:" << all_mantis << std::endl;
    std::cout << "#Skipped genes: " << skipped_genes << std::endl;
    std::cout << "Raptor total #hits:" << all_raptor << std::endl;
    std::cout << "#True positives raptor: " << true_positives << std::endl;
    std::cout << "#False positives raptor: " << false_positives << std::endl;
    std::cout << "#False negatives raptor: " << false_negatives << std::endl;
}
