// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-FileCopyrightText: 2020-2025 Thomas Steinke & Zuse Institute Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#if defined(FPGA_EMULATOR) == defined(FPGA_HARDWARE)
#    error "Either FPGA_EMULATOR or FPGA_HARDWARE have to be defined."
#endif

#include <dlfcn.h> // dlopen, dlsym
#include <filesystem>
#include <format>
#include <string>

#include <cereal/archives/binary.hpp>

#include <raptor/index.hpp>

#include <min_ibf_fpga/backend_sycl/exception_handler.hpp>
#include <min_ibf_fpga/backend_sycl/shared.hpp>

#include <sycl/ext/intel/fpga_extensions.hpp>
#include <sycl/sycl.hpp>

namespace raptor
{

template <size_t chunk_bits, bool profile = false>
class min_ibf_fpga_oneapi
{
    using HostSizeType = min_ibf_fpga::backend_sycl::HostSizeType;
    using Chunk = ac_int<chunk_bits, false>;

    static constexpr bool is_emulator =
#ifdef FPGA_EMULATOR
        true;
#else
        false;
#endif

private:
    raptor::raptor_index<raptor::index_structure::ibf> index;
    size_t technical_bins{}; // TODO: accessor?

    size_t const minimalNumberOfMinimizers{};
    size_t const maximalNumberOfMinimizers{};
    sycl::queue utilitiesQueue;

    std::function<void(void *)> utilities_QueueDeleter = [this](void * ptr)
    {
        sycl::free(ptr, utilitiesQueue);
    };

    template <typename T>
    using unique_utilities_ptr = std::unique_ptr<T, decltype(utilities_QueueDeleter)>;

    template <typename T>
    unique_utilities_ptr<T> malloc_device_utilities(size_t const count)
    {
        return {sycl::malloc_device<T>(count, utilitiesQueue), utilities_QueueDeleter};
    }

    template <typename T>
    unique_utilities_ptr<T> malloc_host_utilities(size_t const count)
    {
        return {sycl::malloc_host<T>(count, utilitiesQueue), utilities_QueueDeleter};
    }

    sycl::queue kernelQueue;
    std::vector<size_t> thresholds_size_t;
    size_t const bufferSizeBytes;
    size_t const numberOfKernelCopys;
    using timeUnit = std::chrono::duration<double, std::milli>;
    std::chrono::nanoseconds constructorDuration{};
    std::array<sycl::event, 2> setupEvents; // 0: transferThresholds, 1: transferIBF
    struct buffer_data
    {
        size_t numberOfQueries;
        std::vector<std::string> ids;

        unique_utilities_ptr<char> queries_host;
        unique_utilities_ptr<HostSizeType> querySizes_host;
        unique_utilities_ptr<Chunk> results_host;

        std::vector<sycl::event> kernelEvents;
    };
    std::array<buffer_data, 2> double_buffer;

    void setup_fpga()
    {
        auto device_selector = []()
        {
            if constexpr (is_emulator)
                return sycl::ext::intel::fpga_emulator_selector_v;
            else
                return sycl::ext::intel::fpga_selector_v;
        }();

        auto create_queue = [&device_selector]()
        {
            if constexpr (profile)
                return sycl::queue(device_selector,
                                   fpga_tools::exception_handler,
                                   sycl::property::queue::enable_profiling());
            else
                return sycl::queue(device_selector, fpga_tools::exception_handler);
        };

        utilitiesQueue = create_queue();

        sycl::device dev = utilitiesQueue.get_device();

        if (!dev.has(sycl::aspect::usm_device_allocations) /*dev.get_info<info::device::usm_device_allocations>()*/)
        {
            std::cerr << "ERROR: The selected device does not support USM device allocations\n";
            std::terminate();
        }
        if (!dev.has(sycl::aspect::usm_host_allocations) /*dev.get_info<info::device::usm_host_allocations>()*/)
        {
            std::cerr << "ERROR: The selected device does not support USM host allocations\n";
            std::terminate();
        }

        kernelQueue = create_queue();
    }

public:
    min_ibf_fpga_oneapi(std::filesystem::path const & index_path,
                        size_t minimalNumberOfMinimizers,
                        size_t maximalNumberOfMinimizers,
                        std::vector<size_t> thresholds,
                        uint8_t const bufferSizeMiB,
                        uint8_t const numberOfKernelCopys) :
        minimalNumberOfMinimizers{minimalNumberOfMinimizers},
        maximalNumberOfMinimizers{maximalNumberOfMinimizers},
        thresholds_size_t{std::move(thresholds)},
        bufferSizeBytes(bufferSizeMiB * 1'048'576),
        numberOfKernelCopys{numberOfKernelCopys}
    {

        setup_fpga();
        std::chrono::steady_clock::time_point start;

        if constexpr (profile)
            start = std::chrono::steady_clock::now();

        {
            std::ifstream archiveStream{index_path, std::ios::binary};
            cereal::BinaryInputArchive archive{archiveStream};
            archive(index);
        }

        technical_bins = seqan::hibf::next_multiple_of_64(index.ibf().bin_count());

        if constexpr (profile)
        {
            std::chrono::steady_clock::time_point const end = std::chrono::steady_clock::now();
            constructorDuration = end - start;
        }
    }

    size_t hash_function_count() const
    {
        return index.ibf().hash_function_count();
    }

    std::string get_library_name() const
    {
        return std::format("libraptor_search_fpga_kernel_w{}_k{}_b{}_kernels{}.so",
                           index.window_size(),
                           index.shape().size(),
                           technical_bins,
                           numberOfKernelCopys);
    }

    void count(std::filesystem::path const & query_path, std::filesystem::path const & output_path)
    {
        std::chrono::steady_clock::time_point countStart;
        std::chrono::steady_clock::time_point hostStart;

        if constexpr (profile)
            countStart = std::chrono::steady_clock::now();

        // Was (data.bit_size() + 63) / 64 * 64 / 8;
        size_t const data_size_bytes = seqan::hibf::divide_and_ceil(index.ibf().bit_size(), 64u) * 8;
        size_t const number_of_chunks = seqan::hibf::divide_and_ceil(data_size_bytes, sizeof(Chunk));
        Chunk const * ibfData_host = reinterpret_cast<Chunk const *>(index.ibf().data());

        auto ibfData_device = malloc_device_utilities<Chunk>(number_of_chunks);
        setupEvents[1] = utilitiesQueue.copy(ibfData_host, ibfData_device.get(), number_of_chunks);

        if constexpr (profile)
        {
            printDuration("IBF I/O:\t", countStart, constructorDuration);
            hostStart = std::chrono::steady_clock::now();
        }

        // RunKernelType is a function pointer to a function that takes multiple arguments and returns void.
        using RunKernelType = void (*)(sycl::queue &,
                                       char const *,
                                       HostSizeType const *,
                                       HostSizeType const,
                                       Chunk const *,
                                       HostSizeType const,
                                       HostSizeType const,
                                       HostSizeType const,
                                       HostSizeType const,
                                       HostSizeType const *,
                                       Chunk *,
                                       std::vector<sycl::event> &);

        std::string const library_name = get_library_name();
        void * handle = dlopen(library_name.c_str(), RTLD_NOW);

        if (!handle)
        {
            throw std::runtime_error{dlerror()};
        }

        dlerror(); // Clear any existing error

        std::string const symbol_name = "RunKernel";
        RunKernelType RunKernel = reinterpret_cast<RunKernelType>(dlsym(handle, symbol_name.c_str()));

        char * error_message = dlerror();

        if (error_message != NULL)
        {
            throw std::runtime_error{error_message};
        }

        static_assert(sizeof(HostSizeType) == sizeof(size_t));

        auto thresholds_device = malloc_device_utilities<HostSizeType>(thresholds_size_t.size());
        HostSizeType const * thresholds_host = reinterpret_cast<HostSizeType const *>(thresholds_size_t.data());
        setupEvents[0] = utilitiesQueue.copy(thresholds_host, thresholds_device.get(), thresholds_size_t.size());

        std::ofstream outputStream(output_path, std::ios::app);
        std::string result_string{};

        const size_t max_numberOfQueries = bufferSizeBytes / MIN_QUERY_LENGTH;
        const size_t chunks_per_result = technical_bins / chunk_bits;

        for (buffer_data & state : double_buffer)
        {
            state.queries_host = malloc_host_utilities<char>(bufferSizeBytes);
            state.querySizes_host = malloc_host_utilities<HostSizeType>(max_numberOfQueries);
            state.results_host = malloc_host_utilities<Chunk>(max_numberOfQueries * chunks_per_result);

            state.kernelEvents.reserve(numberOfKernelCopys * 2 + 2); // +2: Distributor, Collector
            state.numberOfQueries = 0;
        }

        auto const queueToFPGA = [&](buffer_data & state, size_t computeIteration)
        {
            if (computeIteration == 0)
            {
                utilitiesQueue.wait();
                if constexpr (profile)
                {
                    printDurationFromEvent("Thresh. Trans.:\t", setupEvents[0]);
                    printDurationFromEvent("IBF Transfer:\t", setupEvents[1]);
                }
            }

            std::chrono::steady_clock::time_point start;

            if constexpr (profile)
                start = std::chrono::steady_clock::now();

            RunKernel(kernelQueue,
                      state.queries_host.get(),
                      state.querySizes_host.get(),
                      state.numberOfQueries,
                      ibfData_device.get(),
                      index.ibf().bin_size(),
                      std::countl_zero(index.ibf().bin_size()), // TODO: accessor index.ibf().hash_shift()
                      minimalNumberOfMinimizers,
                      maximalNumberOfMinimizers,
                      thresholds_device.get(),
                      state.results_host.get(),
                      state.kernelEvents);

            if constexpr (profile)
                printDuration("Queue:\t\t", start);
        };

        auto const waitOnFPGA = [&](buffer_data & state)
        {
            std::chrono::steady_clock::time_point start;

            if constexpr (profile)
                start = std::chrono::steady_clock::now();

            for (sycl::event e : state.kernelEvents)
                e.wait();

            if constexpr (profile)
                printDuration("Wait:\t\t", start);
        };

        auto const outputResults = [&](buffer_data & state)
        {
            std::chrono::steady_clock::time_point start;

            if constexpr (profile)
                start = std::chrono::steady_clock::now();

            result_string.clear();

            size_t const elements_per_result = technical_bins / 64;
            size_t const elements_per_chunk = chunk_bits / 64;

            if (elements_per_chunk * chunks_per_result != elements_per_result)
                throw std::runtime_error("outputResults: elements_per_result/chunks_per_result mismatch");

            for (size_t queryIndex = 0; queryIndex < state.ids.size(); queryIndex++)
            {
                // TODO: Use string view?
                result_string += state.ids.at(queryIndex).substr(1, std::string::npos) + '\t';

                for (size_t chunkOffset = 0; chunkOffset < chunks_per_result; ++chunkOffset)
                {
                    Chunk & chunk = state.results_host.get()[queryIndex * chunks_per_result + chunkOffset];

                    for (size_t elementOffset = 0; elementOffset < elements_per_chunk; ++elementOffset)
                    {
                        uint64_t const & element = ((uint64_t *)&chunk)[elementOffset];

                        if (!element)
                        {
                            continue;
                        }

                        for (size_t byteOffset = 0; byteOffset < 8; ++byteOffset)
                        {
                            uint8_t const & value = ((uint8_t *)&element)[byteOffset];

                            if (!value)
                            {
                                continue;
                            }

                            uint8_t mask = 1;

                            for (size_t bitOffset = 0; bitOffset < 8; ++bitOffset)
                            {
                                if (value & mask)
                                {
                                    result_string += std::to_string(chunkOffset * chunk_bits + elementOffset * 64
                                                                    + byteOffset * 8 + bitOffset)
                                                   + ",";
                                }
                                mask <<= 1;
                            }
                        }
                    }
                }
                if (auto & last_char = result_string.back(); last_char == ',')
                    last_char = '\n';
                else
                    result_string += '\n';
            }

            outputStream << result_string;

            if constexpr (profile)
                printDuration("Output:\t\t", start);

            if constexpr (profile)
            {
                std::stringstream profilingOutput;

                profilingOutput << "\n"
                                << "Kernels:\t" << getRuntime(state.kernelEvents) << " ms\n";

                for (size_t kernelId = 0; kernelId < state.kernelEvents.size(); ++kernelId)
                    profilingOutput << "Kernel " << kernelId << ":\t" << getRuntime(state.kernelEvents[kernelId])
                                    << " ms\n";

                profilingOutput << "\n";

                std::clog << profilingOutput.str();
            }

            state.kernelEvents.clear();
        };

        std::ifstream inputStream(query_path, std::ios::in | std::ios::binary);

        std::string id;
        std::string query;

        size_t computeIteration = 0;
        size_t currentBufferSize = 0;

        buffer_data * currentBufferData = double_buffer.data();

        while (std::getline(inputStream, id))
        //for (auto && [id, query] : fin)
        {
            std::getline(inputStream, query);

            if (currentBufferSize + query.size() > bufferSizeBytes)
            {
                if (currentBufferData->numberOfQueries == 0)
                    throw std::runtime_error("Buffer size to small");

                if constexpr (profile)
                    printDuration("Host:\t\t", hostStart);

                queueToFPGA(*currentBufferData, computeIteration);
                ++computeIteration;
                currentBufferData = std::addressof(double_buffer[computeIteration % 2]);

                if (computeIteration >= 2)
                {
                    waitOnFPGA(*currentBufferData);
                    outputResults(*currentBufferData); // could run async up to the next kernel start
                }

                if constexpr (profile)
                    hostStart = std::chrono::steady_clock::now();

                currentBufferData->ids.clear();
                currentBufferData->numberOfQueries = 0;

                currentBufferSize = 0;
            }

            currentBufferData->ids.emplace_back(std::move(id));

            currentBufferData->querySizes_host.get()[currentBufferData->numberOfQueries] = query.size();

            // Copy query to queries
            std::memcpy(currentBufferData->queries_host.get() + currentBufferSize, query.data(), query.size());

            currentBufferSize += query.size();
            currentBufferData->numberOfQueries++;

            // Discard two lines from input stream (delimiter and quality)
            for (size_t i = 0; i < 2; ++i)
                inputStream.ignore(std::numeric_limits<std::streamsize>::max(), inputStream.widen('\n'));
        }

        // Handle left over queries
        if (currentBufferData->numberOfQueries)
        {
            if constexpr (profile)
                printDuration("Host:\t\t", hostStart);

            queueToFPGA(*currentBufferData, computeIteration);
            ++computeIteration;
        }

        // Wait for the remaining ones to finish
        for (size_t i = std::min<size_t>(computeIteration, 2ul); i > 0; --i)
        {
            currentBufferData = std::addressof(double_buffer[(computeIteration + i) % 2]);
            waitOnFPGA(*currentBufferData);
            outputResults(*currentBufferData);
        }

        if constexpr (profile)
            printDuration("Count total:\t", countStart);
    }

private:
    static void printDuration(std::string const & label,
                              std::chrono::steady_clock::time_point const & start,
                              std::chrono::steady_clock::duration duration = std::chrono::steady_clock::duration())
    {
        auto const end = std::chrono::steady_clock::now();

        duration += end - start;
        auto const runtime = std::chrono::duration_cast<timeUnit>(duration).count();

        std::stringstream profilingOutput;

        profilingOutput << label << runtime << " ms\n";

        std::clog << profilingOutput.str();
    }

    static void printDurationFromEvent(std::string const & label, sycl::event const event)
    {
        uint64_t const start = event.get_profiling_info<sycl::info::event_profiling::command_start>();
        uint64_t const end = event.get_profiling_info<sycl::info::event_profiling::command_end>();

        uint64_t const duration_ns = end - start;
        float const duration_ms = static_cast<float>(duration_ns) / 1000000;

        std::stringstream profilingOutput;

        profilingOutput << label << duration_ms << " ms\n";

        std::clog << profilingOutput.str();
    }

    // Helper for SYCL profiling info

    static double getRuntime(std::vector<sycl::event> const events)
    {
        auto const comparator_start = [](sycl::event const & lhs, sycl::event const & rhs)
        {
            uint64_t lhsValue, rhsValue;

            lhsValue = lhs.get_profiling_info<sycl::info::event_profiling::command_start>();
            rhsValue = rhs.get_profiling_info<sycl::info::event_profiling::command_start>();

            return lhsValue < rhsValue;
        };

        auto const comparator_end = [](sycl::event const & lhs, sycl::event const & rhs)
        {
            uint64_t lhsValue, rhsValue;

            lhsValue = lhs.get_profiling_info<sycl::info::event_profiling::command_end>();
            rhsValue = rhs.get_profiling_info<sycl::info::event_profiling::command_end>();

            return lhsValue < rhsValue;
        };

        auto const startEvent = *std::min_element(events.begin(),
                                                  events.end(),
                                                  [&](auto const & lhs, auto const & rhs)
                                                  {
                                                      return comparator_start(lhs, rhs);
                                                  });
        auto const endEvent = *std::max_element(events.begin(),
                                                events.end(),
                                                [&](auto const & lhs, auto const & rhs)
                                                {
                                                    return comparator_end(lhs, rhs);
                                                });

        return getRuntime(startEvent, endEvent);
    };

    static double getRuntime(sycl::event const & event)
    {
        return getRuntime(event, event);
    }

    static double getRuntime(sycl::event const & startEvent, sycl::event const & endEvent)
    {
        uint64_t start, end;

        start = startEvent.get_profiling_info<sycl::info::event_profiling::command_start>();
        end = endEvent.get_profiling_info<sycl::info::event_profiling::command_end>();

        std::chrono::nanoseconds const duration{end - start};

        return std::chrono::duration_cast<timeUnit>(duration).count();
    }
};

} // namespace raptor
