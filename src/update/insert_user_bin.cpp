// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <iterator>
#include <optional>

#include <hibf/contrib/robin_hood.hpp>

#include <raptor/argument_parsing/update_arguments.hpp>
#include <raptor/file_reader.hpp>
#include <raptor/index.hpp>

#include "insert/is_fpr_exceeded.hpp"
#include "insert/strong_types.hpp"
#include "subtree.hpp"

namespace raptor::detail
{

std::optional<insert_location> get_location(std::vector<ibf_max> const & max_ibf_sizes,
                                            size_t const kmer_count,
                                            raptor_index<index_structure::hibf> & index);

std::optional<rebuild_location> insert_tb_and_parents(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                                                      insert_location insert_location,
                                                      raptor_index<index_structure::hibf> & index);

/*!\brief Reads the hashes of all files of a user bin into `target`.
 * \details
 * The file type is deduced from the first path; raptor::bin_validator rejects mixing file types within a user bin.
 * The type is deduced per user bin because the files an index was built from and the files being inserted may be of
 * different types.
 */
template <std::output_iterator<uint64_t> it_t>
void hash_paths_into(std::vector<std::string> const & paths,
                     it_t target,
                     seqan3::shape const shape,
                     uint32_t const window_size)
{
    assert(!paths.empty());
    assert(std::ranges::all_of(paths,
                               [](std::filesystem::path const & path)
                               {
                                   return path.extension() == ".minimiser";
                               })
           || std::ranges::none_of(paths,
                                   [](std::filesystem::path const & path)
                                   {
                                       return path.extension() == ".minimiser";
                                   }));

    if (std::filesystem::path{paths.front()}.extension() == ".minimiser")
        file_reader<file_types::minimiser>{}.hash_into(paths, target);
    else
        file_reader<file_types::sequence>{shape, window_size}.hash_into(paths, target);
}

//!\brief Computes the set of k-mers of a user bin, using the shape and window size the index was built with.
robin_hood::unordered_flat_set<uint64_t> compute_kmers(std::vector<std::string> const & user_bin,
                                                       raptor_index<index_structure::hibf> const & index)
{
    robin_hood::unordered_flat_set<uint64_t> kmers{};
    hash_paths_into(user_bin,
                    std::inserter(kmers, kmers.begin()),
                    index.shape(),
                    static_cast<uint32_t>(index.window_size()));
    return kmers;
}

/*!\brief The maximum number of elements a technical bin can hold without exceeding `params.fpr`.
 * \details Inverts the Bloom filter FPR formula: `ceil(BITS / (-HASH / log(1 - exp(log(FPR) / HASH))))`.
 */
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

//!\brief Whether the IBF at `ibf_idx` has been removed by a partial rebuild.
bool is_deleted(raptor_index<index_structure::hibf> const & index, size_t const ibf_idx)
{
    return index.ibf().prev_ibf_id[ibf_idx].ibf_idx == seqan::hibf::bin_kind::deleted;
}

/*!\brief The capacity of a technical bin of each IBF, sorted ascendingly. IBFs deleted by a partial rebuild are omitted.
 * \details
 * This is the capacity implied by the bin size, not the current occupancy: a user bin that does not fit into a
 * single technical bin is split over several of them, so the capacity is what decides how many bins it needs.
 */
std::vector<ibf_max> max_ibf_sizes(raptor_index<index_structure::hibf> const & index)
{
    auto const & ibf_vector = index.ibf().ibf_vector;
    std::vector<ibf_max> max_sizes{};
    max_sizes.reserve(ibf_vector.size());

    for (size_t i = 0; i < ibf_vector.size(); ++i)
    {
        if (is_deleted(index, i))
            continue;
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

//!\brief Collects the IDs of all user bins stored in the subtree rooted at `ibf_idx` into `ub_ids`.
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

/*!\brief Recomputes the layout of the subtree hanging in the merged bin given by `rebuild_location`.
 * \param[in] arguments The update arguments.
 * \param[in] rebuild_location The merged bin whose subtree is rebuilt.
 * \param[in,out] index The index to rebuild in.
 * \details
 * Builds a new HIBF from the user bins of that subtree and splices it into `index` in place of the old subtree. The
 * first IBF of the new subtree reuses the slot of the old one, so the merged bin in the parent stays valid and the
 * rest of the index is untouched. The remaining new IBFs are appended to `ibf_vector`; IBFs of the old subtree that
 * are not reused are retired, see raptor::detail::retire_ibf. Nothing outside the old subtree points into it, so no
 * `next_ibf_id` entry is left dangling.
 */
void partial_rebuild(update_arguments const & arguments,
                     detail::rebuild_location const & rebuild_location,
                     raptor_index<index_structure::hibf> & index)
{
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
        detail::collect_subtree(result, index, child_ibf_id);
        return result;
    }();

    auto input_fn = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        detail::hash_paths_into(index.bin_path()[ub_ids[user_bin_id]],
                                it,
                                index.shape(),
                                static_cast<uint32_t>(index.window_size()));
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

    // "Delete" children. The first IBF is kept, its slot is reused by the new subtree below.
    for (size_t i = 1u; i < overwrite_ibf_ids.size(); ++i)
        detail::retire_ibf(index, overwrite_ibf_ids[i]);

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

    index.mark_resized(child_ibf_id);
    index.sync_resized(true);
}

/*!\brief Rebuilds `index` from scratch from all user bins in its bin path.
 * \param[in] arguments The update arguments.
 * \param[in,out] index The index to rebuild.
 * \details
 * Computes a fresh layout with automatic tmax, overriding the tmax the index was built with: it should always
 * perform better, especially when many user bins were inserted since.
 */
void full_rebuild(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    auto bin_path = index.bin_path();
    auto const shape = index.shape();
    auto const window_size = static_cast<uint32_t>(index.window_size());

    auto input_fn = [&bin_path, shape, window_size](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        detail::hash_paths_into(bin_path[user_bin_id], it, shape, window_size);
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

//!\brief The rebuild an IBF needs because of its number of technical bins.
enum class tmax_check : uint8_t
{
    no_rebuild,
    full_rebuild,
    partial_rebuild
};

/*!\brief Checks whether the IBF at `ibf_idx` has grown too large.
 * \details
 * The top-level IBF must not exceed `tmax`; growing it beyond that requires a full rebuild. A lower-level IBF may
 * exceed `tmax` by raptor::detail::tmax_slack bins before it is rebuilt.
 */
tmax_check check_tmax_rebuild(raptor_index<index_structure::hibf> const & index, size_t const ibf_idx)
{
    size_t const bin_count = index.ibf().ibf_vector[ibf_idx].bin_count();
    size_t const tmax = index.config().tmax;

    if (ibf_idx == 0u)
        return bin_count > tmax ? tmax_check::full_rebuild : tmax_check::no_rebuild;

    return bin_count > tmax + detail::tmax_slack ? tmax_check::partial_rebuild : tmax_check::no_rebuild;
}

//!\brief The rebuild the index needs after a user bin was inserted.
enum class rebuild_kind : uint8_t
{
    none,
    partial,
    full
};

//!\brief What to rebuild after a user bin was inserted, and where.
struct rebuild_decision
{
    rebuild_kind kind{rebuild_kind::none};
    //!\brief The merged bin to rebuild. Only used for raptor::rebuild_kind::partial.
    detail::rebuild_location location{};
};

//!\brief The merged bin the IBF at `ibf_idx` hangs in, as a rebuild location.
detail::rebuild_location parent_location(raptor_index<index_structure::hibf> const & index, size_t const ibf_idx)
{
    assert(ibf_idx != 0u);
    auto const parent = index.ibf().prev_ibf_id[ibf_idx];
    return {.ibf_idx = parent.ibf_idx, .bin_idx = parent.bin_idx};
}

/*!\brief Decides whether the index needs to be rebuilt after a user bin has been inserted.
 * \param[in] index The index.
 * \param[in] insert_location Where the user bin was inserted.
 * \param[in] rebuild_location Where the FPR was exceeded, if it was exceeded at all.
 */
rebuild_decision determine_rebuild(raptor_index<index_structure::hibf> const & index,
                                   detail::insert_location const & insert_location,
                                   std::optional<detail::rebuild_location> const & rebuild_location)
{
    // Exceeding tmax takes precedence: it says that an IBF has too many bins, regardless of their FPR.
    size_t const ibf_idx = rebuild_location.has_value() ? rebuild_location->ibf_idx : insert_location.ibf_idx;

    switch (check_tmax_rebuild(index, ibf_idx))
    {
    case tmax_check::full_rebuild:
        return {.kind = rebuild_kind::full};
    case tmax_check::partial_rebuild:
        return {.kind = rebuild_kind::partial, .location = parent_location(index, ibf_idx)};
    case tmax_check::no_rebuild:
        break;
    }

    if (!rebuild_location.has_value())
        return {};

    // A rebuild location pointing to the top-level IBF has two meanings:
    //   1) The child IBF exceeds FPR, i.e. partial rebuild is necessary.
    //   2) A bin in the top-level IBF exceeds the FPR, i.e. a full rebuild is necessary.
    if (rebuild_location->ibf_idx == 0u && is_fpr_exceeded(index, *rebuild_location))
        return rebuild_decision{.kind = rebuild_kind::full};

    // A bin holding a user bin has no subtree to rebuild, and rebuilding below it would not lower its FPR anyway.
    // What has to be rebuilt is the subtree the containing IBF hangs in, so that the user bin can be laid out
    // anew. The top-level IBF hangs in no merged bin, hence the whole index has to be rebuilt.
    // This branch shouldn't be used because user bins should never exceed the FPR, but we guard anyway.
    if (index.ibf().ibf_bin_to_user_bin_id[rebuild_location->ibf_idx][rebuild_location->bin_idx]
        != seqan::hibf::bin_kind::merged)
    {
        if (rebuild_location->ibf_idx == 0u)
            return {.kind = rebuild_kind::full};

        return {.kind = rebuild_kind::partial, .location = parent_location(index, rebuild_location->ibf_idx)};
    }

    return {.kind = rebuild_kind::partial, .location = *rebuild_location};
}

void insert_user_bin(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    auto full_rebuild_bin_path = index.bin_path();
    full_rebuild_bin_path.insert(full_rebuild_bin_path.end(),
                                 arguments.user_bins_to_insert.begin(),
                                 arguments.user_bins_to_insert.end());

    for (std::vector<std::string> const & user_bin : arguments.user_bins_to_insert)
    {
        if (user_bin.size() > 1u)
            throw std::runtime_error{"Currently not supporting multiple files per UB for insert."};

        auto const kmers = detail::compute_kmers(user_bin, index);

        std::vector<detail::ibf_max> const max_kmers = detail::max_ibf_sizes(index);
        assert(std::ranges::is_sorted(max_kmers));

        auto const insert_location = detail::get_location(max_kmers, kmers.size(), index);

        // A valueless location means that no IBF can take the user bin without exceeding tmax. The index has to be
        // rebuilt in full, hence the user bin is not inserted into the to-be-discarded index first.
        rebuild_decision decision{.kind = rebuild_kind::full};

        if (insert_location.has_value())
        {
            index.append_bin_path(user_bin); // TODO: update_bookkeeping, but it doesn't have the args
            auto const rebuild_location = detail::insert_tb_and_parents(kmers, *insert_location, index);
            decision = determine_rebuild(index, *insert_location, rebuild_location);
        }

        // A full rebuild is always a valid replacement for a partial one.
        if (decision.kind == rebuild_kind::partial && arguments.disable_partial_rebuild)
            decision.kind = rebuild_kind::full;

        switch (decision.kind)
        {
        case rebuild_kind::none:
            break;
        case rebuild_kind::partial:
            partial_rebuild(arguments, decision.location, index);
            break;
        case rebuild_kind::full:
            index.replace_bin_path(std::move(full_rebuild_bin_path));
            full_rebuild(arguments, index);
            return;
        }
    }

    // TODO: If possible, check whether a full rebuild is needed before doing the partial rebuild.
    // TODO: In original code there is a check in partial_rebuild. It shortcircuits if a full rebuild is needed.
}

} // namespace raptor
