#include <raptor/argument_parsing/build.hpp>
#include <raptor/build.hpp>

void init_build_parser(seqan3::argument_parser & parser, build_arguments & arguments)
{
    init_shared_meta(parser);
    init_shared_options(parser, arguments);
    parser.add_positional_option(arguments.bin_path,
                                 "Provide a list of input files (one file per bin). Alternatively, provide a text file "
                                 "containing the paths to the bins (one line per path to a bin). ",
                                 bin_validator{});
    parser.add_option(arguments.out_path,
                      '\0',
                      "output",
                      "Provide an output filepath or an output directory if --compute-minimiser is used.",
                      seqan3::option_spec::required);
    parser.add_option(arguments.size,
                      '\0',
                      "size",
                      "Choose the size of the resulting IBF.",
                      seqan3::option_spec::required,
                      size_validator{"\\d+\\s{0,1}[k,m,g,t,K,M,G,T]"});
    parser.add_option(arguments.hash,
                      '\0',
                      "hash",
                      "Choose the number of hashes.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{1, 5});
    parser.add_flag(arguments.compressed,
                    '\0',
                    "compressed",
                    "Build a compressed IBF.");
    parser.add_flag(arguments.compute_minimiser,
                    '\0',
                    "compute-minimiser",
                    "Computes minimisers using cutoffs from Mantis (Pandey et al.). Does not create the index.");
}

void run_build(seqan3::argument_parser & parser)
{
    build_arguments arguments{};
    init_build_parser(parser, arguments);
    try_parsing(parser);

    // ==========================================
    // Process bin_path
    // ==========================================
    if (arguments.bin_path.size() == 1u) // Either only one bin or a file containing bin paths
    {
        auto & file = arguments.bin_path[0];

        if (file.extension() != ".minimiser")
        {
            try
            {
                bin_validator{}.sequence_file_validator(file);
            }
            catch (seqan3::validation_error const & exception)
            {
                decltype(arguments.bin_path) new_values;
                std::ifstream list_of_files{file};
                std::string line;
                while (std::getline(list_of_files, line))
                    new_values.emplace_back(line);
                while (new_values.back().empty())
                    new_values.pop_back();
                arguments.bin_path = std::move(new_values);
            }
        }
    }

    // ==========================================
    // Various checks.
    // ==========================================

    arguments.bins = arguments.bin_path.size();

    if (parser.is_option_set("window"))
    {
        if (arguments.kmer_size > arguments.window_size)
            throw seqan3::argument_parser_error{"The k-mer size cannot be bigger than the window size."};
    }
    else
        arguments.window_size = arguments.kmer_size;

    if (parser.is_option_set("compute-minimiser"))
    {
        try
        {
            seqan3::output_directory_validator{}(arguments.out_path);
        }
        catch (seqan3::argument_parser_error const & ext)
        {
            std::cerr << "[Error] " << ext.what() << '\n';
            std::exit(-1);
        }
    }
    else
    {
        try
        {
            seqan3::output_file_validator{}(arguments.out_path);
        }
        catch (seqan3::argument_parser_error const & ext)
        {
            std::cerr << "[Error] " << ext.what() << '\n';
            std::exit(-1);
        }
    }

    // ==========================================
    // Process --size.
    // ==========================================
    arguments.size.erase(std::remove(arguments.size.begin(), arguments.size.end(), ' '), arguments.size.end());

    size_t multiplier{};

    switch (std::tolower(arguments.size.back()))
    {
        case 't':
            multiplier = 8ull * 1024ull * 1024ull * 1024ull * 1024ull;
            break;
        case 'g':
            multiplier = 8ull * 1024ull * 1024ull * 1024ull;
            break;
        case 'm':
            multiplier = 8ull * 1024ull * 1024ull;
            break;
        case 'k':
            multiplier = 8ull * 1024ull;
            break;
        default:
            throw seqan3::argument_parser_error{"Use {k, m, g, t} to pass size. E.g., --size 8g."};
    }

    size_t size{};
    std::from_chars(arguments.size.data(), arguments.size.data() + arguments.size.size() - 1, size);
    size *= multiplier;
    arguments.bits = size / (((arguments.bins + 63) >> 6) << 6);

    // ==========================================
    // Dispatch
    // ==========================================
    raptor_build(arguments);
};
