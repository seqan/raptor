#include <raptor/build/build.hpp>
#include <raptor/build/build_from_minimiser.hpp>
#include <raptor/build/compute_minimiser.hpp>
#include <raptor/build/run_program.hpp>

namespace raptor
{

void raptor_build(build_arguments const & arguments)
{
    if (arguments.compute_minimiser)
    {
        compute_minimiser(arguments);
        return;
    }

    std::filesystem::path const first_file_path{arguments.bin_path[0][0]};
    if (first_file_path.extension() == ".minimiser")
    {
        if (arguments.compressed)
            build_from_minimiser<true>(arguments);
        else
            build_from_minimiser<false>(arguments);
        return;
    }

    if (arguments.compressed)
        run_program<true>(arguments);
    else
        run_program<false>(arguments);
    return;
}

} // namespace raptor
