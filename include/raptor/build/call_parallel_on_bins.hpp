#pragma once

#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/utility/views/chunk.hpp>

#include <raptor/shared.hpp>

namespace raptor
{

template <typename algorithm_t>
void call_parallel_on_bins(algorithm_t && worker, build_arguments const & arguments)
{
    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(arguments.bins / arguments.threads),
                                                 8u,
                                                 64u);
    auto chunked_view = seqan3::views::zip(arguments.bin_path, std::views::iota(0u)) |
                        seqan3::views::chunk(chunk_size);
    seqan3::detail::execution_handler_parallel executioner{arguments.threads};
    executioner.bulk_execute(std::move(worker), std::move(chunked_view), [](){});
}

} // namespace raptor
