#include <seqan3/io/all.hpp>

int main(int argc, char ** argv)
{
    if (argc != 2)
        throw std::runtime_error{"provide a fastq file!"};

    seqan3::sequence_file_input fin{argv[1]};

    for (auto & rec : fin)
        std::cout << "parsed record:" << rec.id() << std::endl;
}
