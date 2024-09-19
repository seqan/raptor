// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cstddef>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

namespace raptor::util::generate_reads
{

void print_weights(std::vector<size_t> const & vec, std::string_view const first_line)
{
    std::cout << '\n' << first_line << '\n';
    for (size_t i = 0; i < vec.size(); ++i)
        std::cout << vec[i] << ' ';
    std::cout << '\n';

    auto const sum_of_weights = std::reduce(vec.begin(), vec.end());
    std::cout << "  Sum: " << sum_of_weights << '\n';
}

} // namespace raptor::util::generate_reads
