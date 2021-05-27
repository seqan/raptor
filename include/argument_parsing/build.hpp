#pragma once

#include <argument_parsing/shared.hpp>

void init_build_parser(seqan3::argument_parser & parser, build_arguments & arguments);
void run_build(seqan3::argument_parser & parser);

