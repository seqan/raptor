#include <shared.hpp>

#pragma once

void do_cerealisation_out(std::vector<size_t> const & vec, search_arguments const & arguments);
bool do_cerealisation_in(std::vector<size_t> & vec, search_arguments const & arguments);
std::vector<size_t> precompute_threshold(size_t const pattern_size,
                                         size_t const window_size,
                                         uint8_t const kmer_size,
                                         size_t const errors,
                                         double const tau);
