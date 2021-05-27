#include <raptor/argument_parsing/top_level.hpp>

namespace raptor
{

void init_top_level_parser(seqan3::argument_parser & parser)
{
    init_shared_meta(parser);
    parser.info.description.emplace_back("Binning Directories are a datastruture that can be used in various ways. "
                                         "What's a bin, how can it be used, etc.");

    parser.info.examples = {"./raptor build --help", "./raptor search --help"};
};

} // namespace raptor
