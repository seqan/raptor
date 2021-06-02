#include <raptor/search/run_program_single.hpp>
#include <raptor/search/run_program_single_socks.hpp>
#include <raptor/search/run_program_multiple.hpp>

namespace raptor
{

void raptor_search(search_arguments const & arguments)
{
    if (arguments.parts == 1)
    {
        if (arguments.is_socks)
        {
            if (arguments.compressed)
                run_program_single_socks<true>(arguments);
            else
                run_program_single_socks<false>(arguments);
        }
        else
        {
            if (arguments.compressed)
                run_program_single<true>(arguments);
            else
                run_program_single<false>(arguments);
        }
    }
    else
    {
        if (arguments.compressed)
            run_program_multiple<true>(arguments);
        else
            run_program_multiple<false>(arguments);
    }

    return;
}

} // namespace raptor
