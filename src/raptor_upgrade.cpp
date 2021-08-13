#include <raptor/upgrade/upgrade.hpp>
#include <raptor/upgrade/upgrade_index.hpp>

namespace raptor
{

void raptor_upgrade(upgrade_arguments const & arguments)
{
    if (arguments.compressed)
        upgrade_index<true>(arguments);
    else
        upgrade_index<false>(arguments);
    return;
}

} // namespace raptor
