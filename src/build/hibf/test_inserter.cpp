#include <raptor/build/hibf/build_data.hpp>
#include <raptor/file_reader.hpp>

namespace other
{

void hash_into(raptor::hibf::build_data const & data, size_t const user_bin_number, std::vector<uint64_t> & to)
{
    if (data.arguments.input_is_minimiser)
    {
        raptor::file_reader<raptor::file_types::minimiser> const reader{};
        reader.hash_into(data.filenames[user_bin_number], std::back_inserter(to));
    }
    else
    {
        raptor::file_reader<raptor::file_types::sequence> const reader{data.arguments.shape,
                                                                       data.arguments.window_size};
        reader.hash_into(data.filenames[user_bin_number], std::back_inserter(to));
    }
}

} // namespace other

namespace raptor::hibf
{
void test_inserter(raptor::hibf::build_data const & data, size_t const user_bin_number, std::vector<uint64_t> & to)
{
    other::hash_into(data, user_bin_number, to);
}
} // namespace raptor::hibf
