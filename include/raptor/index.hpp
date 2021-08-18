// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

#include <raptor/shared.hpp>

namespace raptor
{

template <seqan3::data_layout data_layout_mode_ = seqan3::data_layout::uncompressed>
class raptor_index
{
private:
    template <seqan3::data_layout data_layout_mode>
    friend class raptor_index;

    using ibf_t = seqan3::interleaved_bloom_filter<data_layout_mode_>;

    uint64_t window_size_{};
    uint8_t kmer_size_{};
    uint8_t parts_{};
    bool compressed_{};
    std::vector<std::vector<std::string>> bin_path_{};
    ibf_t ibf_{};

public:
    static constexpr seqan3::data_layout data_layout_mode = data_layout_mode_;

    static constexpr uint32_t version{1u};

    raptor_index() = default;
    raptor_index(raptor_index const &) = default;
    raptor_index(raptor_index &&) = default;
    raptor_index & operator=(raptor_index const &) = default;
    raptor_index & operator=(raptor_index &&) = default;
    ~raptor_index() = default;

    explicit raptor_index(window const window_size,
                          kmer const kmer_size,
                          uint8_t const parts,
                          bool const compressed,
                          std::vector<std::vector<std::string>> const & bin_path,
                          ibf_t && ibf)
    :
        window_size_{window_size.v},
        kmer_size_{kmer_size.v},
        parts_{parts},
        compressed_{compressed},
        bin_path_{bin_path},
        ibf_{std::move(ibf)}
    {}

    explicit raptor_index(build_arguments const & arguments)
        requires (data_layout_mode == seqan3::data_layout::uncompressed)
    :
        window_size_{arguments.window_size},
        kmer_size_{arguments.kmer_size},
        parts_{arguments.parts},
        compressed_{arguments.compressed},
        bin_path_{arguments.bin_path},
        ibf_{seqan3::bin_count{arguments.bins},
             seqan3::bin_size{arguments.bits / arguments.parts},
             seqan3::hash_function_count{arguments.hash}}
    {}

    explicit raptor_index(raptor_index<seqan3::data_layout::uncompressed> const & other)
        requires (data_layout_mode == seqan3::data_layout::compressed)
    {
        window_size_ = other.window_size_;
        kmer_size_ = other.kmer_size_;
        parts_ = other.parts_;
        compressed_ = true;
        bin_path_ = other.bin_path_;
        ibf_ = ibf_t{other.ibf_};
    }

    explicit raptor_index(raptor_index<seqan3::data_layout::uncompressed> && other)
        requires (data_layout_mode == seqan3::data_layout::compressed)
    {
        window_size_ = std::move(other.window_size_);
        kmer_size_ = std::move(other.kmer_size_);
        parts_ = std::move(other.parts_);
        compressed_ = true;
        bin_path_ = std::move(other.bin_path_);
        ibf_ = std::move(ibf_t{std::move(other.ibf_)});
    }

    uint64_t window_size() const
    {
        return window_size_;
    }

    uint8_t kmer_size() const
    {
        return kmer_size_;
    }

    uint8_t parts() const
    {
        return parts_;
    }

    bool compressed() const
    {
        return compressed_;
    }

    std::vector<std::vector<std::string>> const & bin_path() const
    {
        return bin_path_;
    }

    ibf_t & ibf()
    {
        return ibf_;
    }

    ibf_t const & ibf() const
    {
        return ibf_;
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     * \param[in] version Index version.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive, uint32_t const version)
    {
        if (version == 1u)
        {
            try
            {
                archive(window_size_);
                seqan3::shape shape = kmer_size_ ? seqan3::shape{seqan3::ungapped{kmer_size_}} : seqan3::shape{};
                archive(shape);
                kmer_size_ = shape.size();
                archive(parts_);
                archive(compressed_);
                if ((data_layout_mode == seqan3::data_layout::compressed && !compressed_) ||
                    (data_layout_mode == seqan3::data_layout::uncompressed && compressed_))
                {
                    throw seqan3::argument_parser_error{"Data layouts of serialised and specified index differ."};
                }
                archive(bin_path_);
                archive(ibf_);
            }
            catch (std::exception const & e)
            {
                throw seqan3::argument_parser_error{"Cannot read index: " + std::string{e.what()}};
            }
        }
        else
        {
            throw seqan3::argument_parser_error{"Unsupported index version. Check raptor upgrade."};
        }
    }

    /* \brief Serialisation support function. Do not load the actual data.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_input_archive.
     * \param[in] archive The archive being serialised from/to.
     * \param[in] version Index version.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <seqan3::cereal_input_archive archive_t>
    void load_parameters(archive_t & archive)
    {
        uint32_t version{};
        archive(version);
        if (version == 1u)
        {
            try
            {
                archive(window_size_);
                seqan3::shape shape{};
                archive(shape);
                kmer_size_ = shape.size();
                archive(parts_);
                archive(compressed_);
                archive(bin_path_);
            }
            catch (std::exception const & e)
            {
                throw seqan3::argument_parser_error{"Cannot read index: " + std::string{e.what()}};
            }
        }
        else
        {
            throw seqan3::argument_parser_error{"Unsupported index version. Check raptor upgrade."};
        }
    }
    //!\endcond

};

} // namespace raptor

CEREAL_CLASS_VERSION(raptor::raptor_index<seqan3::data_layout::uncompressed>, raptor::raptor_index<>::version);
CEREAL_CLASS_VERSION(raptor::raptor_index<seqan3::data_layout::compressed>, raptor::raptor_index<>::version);
