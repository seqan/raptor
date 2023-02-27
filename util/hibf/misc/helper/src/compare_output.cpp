#include <cassert>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <unordered_map>
#include <vector>

int main(int argc, char ** argv)
{
    if (argc != 4)
        throw std::runtime_error{
            "Please provide user_bin.id mantis.ready and raptor.ready"}; // $FILENAME_USER_BIN_IDS $FILENAME_MANTIS_READY_TO_COMPARE $FILENAME_RAPTOR_READY_TO_COMPARE

    std::ifstream user_bin_ids_file{argv[1]};

    std::cout << "Reading in " << argv[1] << "... ";
    std::string line;
    std::unordered_map<std::string, uint64_t> user_bin_ids;
    while (std::getline(user_bin_ids_file, line))
    {
        std::string value{line.begin(), line.begin() + line.find('\t')};
        std::string key{line.begin() + line.find('\t') + 1, line.end()};
        uint64_t value_as_number = std::atoi(value.data());
        user_bin_ids.emplace(key, value_as_number);
    }
    std::cout << "Done " << std::endl;

    std::ifstream mantis_result{argv[2]};
    std::ifstream raptor_result{argv[3]};

    if (!mantis_result.good())
        throw std::runtime_error{"Could not open file " + std::string{argv[2]}};
    if (!raptor_result.good())
        throw std::runtime_error{"Could not open file " + std::string{argv[3]}};

    std::cout << "Processing Results from " << argv[2] << " and " << argv[3] << "... ";

    std::string mantis_line, raptor_line;
    std::ofstream false_positives_file{"raptor.fps"};
    std::ofstream false_negatives_file{"raptor.fns"};
    uint64_t false_positives{0};
    uint64_t false_negatives{0};
    uint64_t line_no{0};
    uint64_t all_mantis{0};
    uint64_t all_raptor{0};

    while (std::getline(mantis_result, mantis_line) && std::getline(raptor_result, raptor_line))
    {
        std::string mantis_query_name{mantis_line.begin(), mantis_line.begin() + mantis_line.find('\t')};
        std::string raptor_query_name{raptor_line.begin(), raptor_line.begin() + raptor_line.find('\t')};
        if (mantis_query_name != raptor_query_name)
            throw std::runtime_error{"Query names do not match, something went wrong"};

        std::string query_name{mantis_query_name.begin(),
                               mantis_query_name.begin() + mantis_query_name.find("genomic") + 7};
        uint64_t query_id = user_bin_ids[query_name];
        bool found_query_id_in_mantis{false};
        bool found_query_id_in_raptor{false};

        std::string mantis_fields{mantis_line.begin() + mantis_line.find('\t') + 1, mantis_line.end()};
        std::string raptor_fields{raptor_line.begin() + raptor_line.find('\t') + 1, raptor_line.end()};
        auto mantis_fields_view = mantis_fields | std::views::split(',');
        auto raptor_fields_view = raptor_fields | std::views::split(',');
        auto mantis_it = mantis_fields_view.begin();
        auto raptor_it = raptor_fields_view.begin();

        while (mantis_it != mantis_fields_view.end() && raptor_it != raptor_fields_view.end())
        {
            std::string mantis_str = [](auto v)
            {
                std::string s;
                for (auto c : v)
                    s.push_back(c);
                return s;
            }(*mantis_it);
            std::string raptor_str = [](auto v)
            {
                std::string s;
                for (auto c : v)
                    s.push_back(c);
                return s;
            }(*raptor_it);
            // Should also work:
            // std::string_view mantis_str{*mantis_it};
            // std::string_view raptor_str{*raptor_it};
            uint64_t mantis_value = std::atoi(mantis_str.data());
            uint64_t raptor_value = std::atoi(raptor_str.data());

            found_query_id_in_mantis = found_query_id_in_mantis || mantis_value == query_id;
            found_query_id_in_raptor = found_query_id_in_raptor || raptor_value == query_id;
            // Was:
            // if (mantis_value == query_id)
            //     found_query_id_in_mantis = true;
            // if (raptor_value == query_id)
            //     found_query_id_in_raptor = true;

            if (mantis_value != raptor_value) // If mantis results are empty, then...?
            {
                if (mantis_value < raptor_value)
                {
                    if (raptor_value != query_id)
                    {
                        false_negatives_file << mantis_query_name << ":" << mantis_value << '\n';
                        ++false_negatives;
                    }
                    ++all_mantis;
                    ++mantis_it;
                }
                else
                {
                    if (raptor_value != query_id)
                    {
                        false_positives_file << raptor_query_name << ":" << raptor_value << '\n';
                        ++false_positives;
                    }
                    ++all_raptor;
                    ++raptor_it;
                }
            }
            else
            {
                ++all_mantis;
                ++all_raptor;
                ++mantis_it;
                ++raptor_it;
            }
        }

        while (mantis_it != mantis_fields_view.end()) // process the rest of mantis
        {
            std::string mantis_str = [](auto v)
            {
                std::string s;
                for (auto c : v)
                    s.push_back(c);
                return s;
            }(*mantis_it);
            std::string query_name{mantis_query_name.begin(),
                                   mantis_query_name.begin() + mantis_query_name.find("genomic") + 7};
            // uint64_t query_id = user_bin_ids[query_name];
            uint64_t mantis_value = std::atoi(mantis_str.data());
            found_query_id_in_mantis = found_query_id_in_mantis || mantis_value == query_id;
            ++false_negatives;
            false_negatives_file << query_name << ":" << mantis_value << '\n';
            ++all_mantis;
            ++mantis_it;
        }

        while (raptor_it != raptor_fields_view.end()) // process the rest of raptor if any
        {
            std::string raptor_str = [](auto v)
            {
                std::string s;
                for (auto c : v)
                    s.push_back(c);
                return s;
            }(*raptor_it);
            std::string query_name{mantis_query_name.begin(),
                                   mantis_query_name.begin() + mantis_query_name.find("genomic") + 7};
            // uint64_t query_id = user_bin_ids[query_name];
            uint64_t raptor_value = std::atoi(raptor_str.data());
            if (raptor_value == query_id)
            {
                found_query_id_in_raptor = true;
            }
            else
            {
                false_positives_file << query_name << ":" << raptor_value << '\n';
                ++false_positives;
            }
            ++all_raptor;
            ++raptor_it;
        }

        if (!found_query_id_in_mantis)
            std::cerr << "Warning in line " << line_no << ": Could not find query " << mantis_query_name << "("
                      << query_name << ":" << query_id << ") in its respective gemnome in mantis." << std::endl;
        if (!found_query_id_in_raptor)
            std::cerr << "Warning in line " << line_no << ": Could not find query " << raptor_query_name << "("
                      << query_name << ":" << query_id << ") in its respective gemnome in raptor." << std::endl;

        ++line_no;
    }

    while (std::getline(mantis_result, mantis_line))
        std::cerr << "WARNING: Missing line of mantis in comparison: " << mantis_line;
    while (std::getline(raptor_result, raptor_line))
        std::cerr << "WARNING: Missing line of raptor in comparison: " << raptor_line;

    std::cout << std::endl;
    std::cout << "Mantis total #hits:" << all_mantis << std::endl;
    std::cout << "Raptor total #hits:" << all_raptor << std::endl;
    std::cout << "#False positives raptor: " << false_positives << std::endl;
    std::cout << "#False negatives raptor: " << false_negatives << std::endl;
}
