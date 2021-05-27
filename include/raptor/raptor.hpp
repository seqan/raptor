#pragma once

#include <raptor/shared.hpp>

namespace raptor
{

void try_parsing(seqan3::argument_parser & parser);
void init_top_level_parser(seqan3::argument_parser & parser);
void run_build(seqan3::argument_parser & parser);
void run_search(seqan3::argument_parser & parser);

} // namespace raptor
