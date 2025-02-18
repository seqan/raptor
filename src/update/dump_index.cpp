// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::dump_index.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/contrib/std/enumerate_view.hpp>

#include <raptor/update/dump_index.hpp>

namespace raptor
{

void dump_index(raptor_index<index_structure::hibf> const & index)
{
    dump_index(index.ibf());
}

void dump_index(seqan::hibf::hierarchical_interleaved_bloom_filter const & hibf)
{
    std::cerr << "\nDumping index\n";
    // std::cerr << "Window size: " << index.window_size() << '\n';
    // std::cerr << "Shape:       " << index.shape().to_string() << '\n';
    // std::cerr << "Parts:       " << index.parts() << '\n';
    // std::cerr << "FPR:         " << index.fpr() << '\n';
    // std::cerr << "Type:        " << (index.is_hibf() ? "HIBF" : "IBF") << '\n';
    // std::cerr << "Bin path:    " << index.bin_path().front().front() << '\n';
    for (auto const [i, to_user_bin_id] : seqan::stl::views::enumerate(hibf.ibf_bin_to_user_bin_id))
    {
        std::cerr << "User bin " << i << " [" << to_user_bin_id.size() << "]:   [";
        char sep{};

        for (auto const [j, val] : seqan::stl::views::enumerate(to_user_bin_id))
        {
            switch (val)
            {
            case seqan::hibf::bin_kind::deleted:
                std::cerr << sep << 'D';
                break;
            case seqan::hibf::bin_kind::merged:
                std::cerr << sep << 'M' << hibf.next_ibf_id[i][j];
                break;
            default:
                std::cerr << sep << val;
            }

            sep = ',';
        }
        std::cerr << "]\n";
    }

    for (auto const [i, ibf] : seqan::stl::views::enumerate(hibf.ibf_vector))
    {
        std::cerr << "IBF " << i << '\n';
        std::cerr << "  Num bins: " << ibf.bin_count() << '\n';
        std::cerr << "  Occupancy[" << ibf.occupancy.size() << "]:   [";
        char sep{};
        for (auto const val : ibf.occupancy)
        {
            std::cerr << sep << val;
            sep = ',';
        }
        std::cerr << "]\n";
    }
    std::cerr << std::flush;
}

} // namespace raptor
