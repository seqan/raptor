#pragma once

#include <raptor/argument_parsing/shared.hpp>

namespace raptor
{

void init_build_parser(seqan3::argument_parser & parser, build_arguments & arguments);
void run_build(seqan3::argument_parser & parser);

} // namespace raptor
