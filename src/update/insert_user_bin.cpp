// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/contrib/robin_hood.hpp>

#include <raptor/argument_parsing/update_arguments.hpp>
#include <raptor/file_reader.hpp>
#include <raptor/index.hpp>
#include <raptor/update/dump_index.hpp> // DEBUG

#include "insert/is_fpr_exceeded.hpp"
#include "insert/strong_types.hpp"

namespace raptor::detail
{

insert_location get_location(std::vector<ibf_max> const & max_ibf_sizes,
                             size_t const kmer_count,
                             raptor_index<index_structure::hibf> & index);

rebuild_location insert_tb_and_parents(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                                       insert_location insert_location,
                                       raptor_index<index_structure::hibf> & index);

robin_hood::unordered_flat_set<uint64_t> compute_kmers(std::filesystem::path const & ub_file,
                                                       raptor_index<index_structure::hibf> const & index)
{
    robin_hood::unordered_flat_set<uint64_t> kmers{};
    raptor::file_reader<raptor::file_types::sequence> reader{index.shape(), static_cast<uint32_t>(index.window_size())};
    reader.hash_into(ub_file, std::inserter(kmers, kmers.begin()));
    return kmers;
}

// ceil(BITS / (-HASH / log(1 - exp(log(FPR) / HASH))))
size_t max_elements(max_elements_parameters const & params)
{
    assert(params.hash_count > 0);
    assert(params.fpr > 0.0);
    assert(params.fpr < 1.0);

    double const numerator{params.bin_size * std::log(1 - std::exp(std::log(params.fpr) / params.hash_count))};
    double const denominator{-static_cast<double>(params.hash_count)};

    double const result{std::ceil(numerator / denominator)};
    return result;
}

// Shouldn't this be the occupancy? However, bins might get cleared.
// Answer: Split bin correction :)
std::vector<ibf_max> max_ibf_sizes(raptor_index<index_structure::hibf> const & index)
{
    auto const & ibf_vector = index.ibf().ibf_vector;
    std::vector<ibf_max> max_sizes{};
    max_sizes.reserve(ibf_vector.size());

    for (size_t i = 0; i < ibf_vector.size(); ++i)
    {
        auto const & ibf = ibf_vector[i];
        size_t const max_kmers = max_elements({.fpr = index.fpr(), //
                                               .hash_count = ibf.hash_function_count(),
                                               .bin_size = ibf.bin_size()});
        max_sizes.push_back({.max_elements = max_kmers, .ibf_idx = i});
    }
    std::ranges::sort(max_sizes);
    return max_sizes;
}

} // namespace raptor::detail

namespace raptor
{

void get_ubs(robin_hood::unordered_flat_set<uint64_t> & ub_ids,
             raptor_index<index_structure::hibf> & index,
             size_t const ibf_idx)
{
    auto const & user_bin_ids = index.ibf().ibf_bin_to_user_bin_id[ibf_idx];

    for (size_t i = 0; i < user_bin_ids.size(); ++i)
    {
        size_t const ub = user_bin_ids[i];
        switch (ub)
        {
        case seqan::hibf::bin_kind::merged:
            get_ubs(ub_ids, index, index.ibf().next_ibf_id[ibf_idx][i]);
            break;
        case seqan::hibf::bin_kind::deleted:
            break;
        default:
            ub_ids.emplace(ub);
        }
    }
}

void get_overwrite_ibf(std::vector<size_t> & result, raptor_index<index_structure::hibf> & index, size_t const ibf_idx)
{
    result.push_back(ibf_idx);
    auto const & user_bin_ids = index.ibf().ibf_bin_to_user_bin_id[ibf_idx];
    for (size_t i = 0; i < user_bin_ids.size(); ++i)
    {
        if (user_bin_ids[i] == seqan::hibf::bin_kind::merged)
            get_overwrite_ibf(result, index, index.ibf().next_ibf_id[ibf_idx][i]);
    }
}

void partial_rebuild(update_arguments const & arguments,
                     detail::rebuild_location const & rebuild_location,
                     raptor_index<index_structure::hibf> & index)
{
    std::cout << "Partial Rebuild\n";
    assert(index.ibf().ibf_bin_to_user_bin_id[rebuild_location.ibf_idx][rebuild_location.bin_idx]
           == seqan::hibf::bin_kind::merged);
    size_t const child_ibf_id = index.ibf().next_ibf_id[rebuild_location.ibf_idx][rebuild_location.bin_idx];

    std::vector<size_t> const ub_ids = [&]()
    {
        robin_hood::unordered_flat_set<uint64_t> ub_ids{};
        get_ubs(ub_ids, index, child_ibf_id);
        return std::vector<size_t>{ub_ids.begin(), ub_ids.end()};
    }();

    std::vector<size_t> const overwrite_ibf_ids = [&]()
    {
        std::vector<size_t> result{};
        get_overwrite_ibf(result, index, child_ibf_id);
        return result;
    }();

    auto input_fn = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        raptor::file_reader<raptor::file_types::sequence> reader{index.shape(),
                                                                 static_cast<uint32_t>(index.window_size())};
        reader.hash_into(index.bin_path()[ub_ids[user_bin_id]], it);
    };

    seqan::hibf::config config{index.config()};
    // config.tmax = 0u;
    config.input_fn = input_fn;
    config.number_of_user_bins = ub_ids.size();
    config.threads = arguments.threads;
    // config.validate_and_set_defaults();

    seqan::hibf::hierarchical_interleaved_bloom_filter subindex{config};

    auto & original_hibf = index.ibf();
    size_t const offset = original_hibf.ibf_vector.size() - 1u;

    // "Delete" children:
    for (size_t i = 1u; i < overwrite_ibf_ids.size(); ++i)
    {
        size_t const ibf_id = overwrite_ibf_ids[i];

        original_hibf.ibf_vector[ibf_id] = seqan::hibf::interleaved_bloom_filter{};
        original_hibf.next_ibf_id[ibf_id].clear();
        original_hibf.prev_ibf_id[ibf_id] = {seqan::hibf::bin_kind::deleted, seqan::hibf::bin_kind::deleted};
        original_hibf.ibf_bin_to_user_bin_id[ibf_id].clear();
    }

    // Handle the first IBF
    original_hibf.ibf_vector[child_ibf_id] = std::move(subindex.ibf_vector[0]);

    auto & first_ibf_next_ibf_id = subindex.next_ibf_id[0];
    std::ranges::for_each(first_ibf_next_ibf_id,
                          [&](auto & id)
                          {
                              switch (id)
                              {
                              case 0:
                                  id = child_ibf_id;
                                  break;
                              default:
                                  id += offset;
                              }
                          });
    original_hibf.next_ibf_id[child_ibf_id] = std::move(first_ibf_next_ibf_id);

    auto & first_ibf_bin_to_user_bin_id = subindex.ibf_bin_to_user_bin_id[0];
    std::ranges::for_each(first_ibf_bin_to_user_bin_id,
                          [&](auto & id)
                          {
                              switch (id)
                              {
                              case seqan::hibf::bin_kind::deleted:
                                  break;
                              case seqan::hibf::bin_kind::merged:
                                  break;
                              default:
                                  id = ub_ids[id];
                              }
                          });
    original_hibf.ibf_bin_to_user_bin_id[child_ibf_id] = std::move(first_ibf_bin_to_user_bin_id);
    // Prev_ibf_id does not change for first IBF

    assert(subindex.ibf_vector[0].data() == nullptr);
    assert(subindex.next_ibf_id[0].empty());
    assert(subindex.ibf_bin_to_user_bin_id[0].empty());

    // Handle the rest of the IBFs
    for (size_t i = 1; i < subindex.ibf_vector.size(); ++i)
    {
        original_hibf.ibf_vector.push_back(std::move(subindex.ibf_vector[i]));

        auto & ibf_next_ibf_id = subindex.next_ibf_id[i];
        std::ranges::for_each(ibf_next_ibf_id,
                              [&](auto & id)
                              {
                                  id += offset;
                              });
        original_hibf.next_ibf_id.push_back(std::move(ibf_next_ibf_id));

        auto & ibf_bin_to_user_bin_id = subindex.ibf_bin_to_user_bin_id[i];
        std::ranges::for_each(ibf_bin_to_user_bin_id,
                              [&](auto & id)
                              {
                                  switch (id)
                                  {
                                  case seqan::hibf::bin_kind::deleted:
                                      break;
                                  case seqan::hibf::bin_kind::merged:
                                      break;
                                  default:
                                      id = ub_ids[id];
                                  }
                              });
        original_hibf.ibf_bin_to_user_bin_id.push_back(std::move(ibf_bin_to_user_bin_id));

        auto prev_idx = subindex.prev_ibf_id[i];
        if (prev_idx.ibf_idx == 0)
            prev_idx.ibf_idx = child_ibf_id;
        else
            prev_idx.ibf_idx += offset;
        original_hibf.prev_ibf_id.push_back(prev_idx);
    }
}

static constexpr bool consider_lower_level_tmax{false};

void full_rebuild(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    std::cout << "Full Rebuild\n";
    auto bin_path = index.bin_path();
    auto const shape = index.shape();
    auto const window_size = static_cast<uint32_t>(index.window_size());

    auto input_fn = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        raptor::file_reader<raptor::file_types::sequence> reader{shape, window_size};
        reader.hash_into(bin_path[user_bin_id], it);
    };

    seqan::hibf::config config{index.config()};
    // Force auto-tmax, overriding user tmax. It should always perform better, especially when inserting many bins.
    config.tmax = 0u;
    config.input_fn = std::move(input_fn);
    config.number_of_user_bins = bin_path.size();
    config.threads = arguments.threads;
    config.validate_and_set_defaults();

    index = {};
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};
    index = raptor_index<index_structure::hibf>{window{window_size},
                                                shape,
                                                1u,
                                                std::move(bin_path),
                                                config,
                                                std::move(hibf)};
}

enum class tmax_check : uint8_t
{
    no_rebuild = 0u,
    full_rebuild = 1u,
    partial_rebuild = 2u
};

tmax_check check_tmax_rebuild(update_arguments const & arguments,
                              raptor_index<index_structure::hibf> & index,
                              size_t const ibf_idx)
{
    if (index.ibf().ibf_vector[ibf_idx].bin_count() > seqan::hibf::next_multiple_of_64(index.config().tmax))
    {
        if (ibf_idx == 0u)
        {
            return tmax_check::full_rebuild;
        }
        else if constexpr (consider_lower_level_tmax)
        {
            auto const parent = index.ibf().prev_ibf_id[ibf_idx];
            partial_rebuild(arguments, detail::rebuild_location{parent.ibf_idx, parent.bin_idx}, index);
            return tmax_check::partial_rebuild;
        }

        return tmax_check::no_rebuild;
    }

    return tmax_check::no_rebuild;
}

void insert_user_bin(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    auto full_rebuild_bin_path = index.bin_path();
    full_rebuild_bin_path.insert(full_rebuild_bin_path.end(),
                                 arguments.user_bins_to_insert.begin(),
                                 arguments.user_bins_to_insert.end());

    for (auto const & ub : arguments.user_bins_to_insert)
    {
        if (ub.size() > 1u)
            throw std::runtime_error{"Currently not supporting multiple files per UB for insert."};

        for (auto const & path : ub)
        {
            auto const kmers = detail::compute_kmers(path, index);
            size_t const kmer_count = kmers.size();

            std::vector<detail::ibf_max> const max_kmers = detail::max_ibf_sizes(index);
            assert(std::ranges::is_sorted(max_kmers));

            auto const insert_location = detail::get_location(max_kmers, kmer_count, index);
            index.append_bin_path({path}); // TODO: update_bookkeeping, but it doesn't have the args
            auto const rebuild_location = detail::insert_tb_and_parents(kmers, insert_location, index);

            if (rebuild_location.ibf_idx != std::numeric_limits<size_t>::max())
            {
                if (check_tmax_rebuild(arguments, index, rebuild_location.ibf_idx) == tmax_check::no_rebuild)
                {
                    if (rebuild_location.ibf_idx == 0u && is_fpr_exceeded(index, rebuild_location))
                    {
                        index.replace_bin_path(std::move(full_rebuild_bin_path));
                        full_rebuild(arguments, index);
                        return;
                    }
                    else
                    {
                        // some downstream fpr too high
                        partial_rebuild(arguments, rebuild_location, index);
                    }
                }
            }
            else
            {
                // tmax too high
                if (check_tmax_rebuild(arguments, index, insert_location.ibf_idx) == tmax_check::full_rebuild)
                {
                    index.replace_bin_path(std::move(full_rebuild_bin_path));
                    full_rebuild(arguments, index);
                    return;
                }
            }
        }
    }

    // TODO: If possible, check whether a full rebuild is needed before doing the partial rebuild.
    // TODO: In original code there is a check in partial_rebuild. It shortcircuits if a full rebuild is needed.
}

} // namespace raptor
