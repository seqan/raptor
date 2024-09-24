// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::raptor_index.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cereal/types/string.hpp>

#include <sharg/exceptions.hpp>

#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/strong_types.hpp>

namespace raptor
{

namespace index_structure
{

using ibf = seqan::hibf::interleaved_bloom_filter;
using hibf = seqan::hibf::hierarchical_interleaved_bloom_filter;

template <typename index_t>
concept is_ibf = std::same_as<index_t, index_structure::ibf>;

template <typename index_t>
concept is_hibf = std::same_as<index_t, index_structure::hibf>;

template <typename index_t>
concept is_valid = is_ibf<index_t> || is_hibf<index_t>;

} // namespace index_structure

class index_upgrader;

template <index_structure::is_valid data_t = index_structure::ibf>
class raptor_index
{
private:
    template <index_structure::is_valid friend_data_t>
    friend class raptor_index;

    uint64_t window_size_{};
    seqan3::shape shape_{};
    uint8_t parts_{};
    std::vector<std::vector<std::string>> bin_path_{};
    bool is_hibf_{index_structure::is_hibf<data_t>};
    double fpr_{};
    seqan::hibf::config config_{};
    data_t ibf_{};

public:
    static constexpr uint32_t version{3u};

    raptor_index() = default;
    raptor_index(raptor_index const &) = default;
    raptor_index(raptor_index &&) = default;
    raptor_index & operator=(raptor_index const &) = default;
    raptor_index & operator=(raptor_index &&) = default;
    ~raptor_index() = default;

    explicit raptor_index(window const window_size,
                          seqan3::shape const shape,
                          uint8_t const parts,
                          std::vector<std::vector<std::string>> const & bin_path,
                          seqan::hibf::config const & config,
                          data_t && ibf)
        requires index_structure::is_hibf<data_t>
        :
        window_size_{window_size.v},
        shape_{shape},
        parts_{parts},
        bin_path_{bin_path},
        fpr_{config.maximum_fpr},
        config_{config},
        ibf_{std::move(ibf)}
    {}

    explicit raptor_index(build_arguments const & arguments)
        requires index_structure::is_ibf<data_t>
        :
        window_size_{arguments.window_size},
        shape_{arguments.shape},
        parts_{arguments.parts},
        bin_path_{arguments.bin_path},
        fpr_{arguments.fpr},
        ibf_{seqan::hibf::bin_count{arguments.bins},
             seqan::hibf::bin_size{arguments.bits / arguments.parts},
             seqan::hibf::hash_function_count{arguments.hash}}
    {}

    uint64_t window_size() const
    {
        return window_size_;
    }

    seqan3::shape shape() const
    {
        return shape_;
    }

    uint8_t parts() const
    {
        return parts_;
    }

    std::vector<std::vector<std::string>> const & bin_path() const
    {
        return bin_path_;
    }

    double fpr() const
    {
        return fpr_;
    }

    auto const & config() const
        requires index_structure::is_hibf<data_t>
    {
        return config_;
    }

    bool is_hibf() const
    {
        return is_hibf_;
    }

    data_t & ibf()
    {
        return ibf_;
    }

    data_t const & ibf() const
    {
        return ibf_;
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly.
     * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
     */
    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        uint32_t parsed_version{raptor_index<>::version};
        archive(parsed_version);
        if (parsed_version == raptor_index<>::version)
        {
            try
            {
                archive(window_size_);
                archive(shape_);
                archive(parts_);
                archive(bin_path_);
                archive(fpr_);
                archive(is_hibf_);
                archive(config_);
                archive(ibf_);
            }
            // GCOVR_EXCL_START
            catch (std::exception const & e)
            {
                throw sharg::parser_error{"Cannot read index: " + std::string{e.what()}};
            }
            // GCOVR_EXCL_STOP
        }
        else
        {
            throw sharg::parser_error{"Unsupported index version. Check raptor upgrade."}; // GCOVR_EXCL_LINE
        }
    }

    /* \brief Serialisation support function. Do not load the actual data.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_input_archive.
     * \param[in] archive The archive being serialised from/to.
     * \param[in] version Index version.
     *
     * \attention These functions are never called directly.
     * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
     */
    template <seqan3::cereal_input_archive archive_t>
    void load_parameters(archive_t & archive)
    {
        uint32_t parsed_version{};
        archive(parsed_version);
        if (parsed_version == version)
        {
            try
            {
                archive(window_size_);
                archive(shape_);
                archive(parts_);
                archive(bin_path_);
                archive(fpr_);
                archive(is_hibf_);
                archive(config_);
            }
            // GCOVR_EXCL_START
            catch (std::exception const & e)
            {
                throw sharg::parser_error{"Cannot read index: " + std::string{e.what()}};
            }
            // GCOVR_EXCL_STOP
        }
        else
        {
            throw sharg::parser_error{"Unsupported index version. Check raptor upgrade."}; // GCOVR_EXCL_LINE
        }
    }
    //!\endcond
};

} // namespace raptor
