#pragma once

#include <cstdint>
#include <vector>

namespace raptor::detail
{

std::vector<size_t> pascal_row(size_t n);

std::tuple<double, std::vector<double>>
simple_model(size_t const kmer_size, std::vector<double> const & proba_x, std::vector<double> const & indirect_errors);

void impl(size_t const minimizers_left,
          std::vector<double> const & proba,
          std::vector<size_t> error_distribution,
          size_t const current_error_index,
          double & result);

double enumerate_all_errors(size_t const number_of_minimizers, size_t const errors, std::vector<double> const & proba);

std::vector<double> destroyed_indirectly_by_error(size_t const pattern_size,
                                                  size_t const window_size,
                                                  uint8_t const kmer_size);

} // namespace raptor::detail
