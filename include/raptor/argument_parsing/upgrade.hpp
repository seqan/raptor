#pragma once

#include <raptor/argument_parsing/shared.hpp>

namespace raptor
{

void init_upgrade_parser(seqan3::argument_parser & parser, upgrade_arguments & arguments);
void run_upgrade(seqan3::argument_parser & parser);

} // namespace raptor
