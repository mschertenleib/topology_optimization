#include "utility.hpp"

#include <cassert>

namespace
{

std::vector<Profile_entry> profile_entries {};
std::uint32_t current_level {};

} // namespace

void profile_begin_frame()
{
    const auto old_capacity = profile_entries.capacity();
    profile_entries.clear();
    // Clearing the vector should leave its capacity unchanged, but this is not
    // guaranteed. If the capacity did change, just reserve again. Having one
    // allocation per frame before starting the profiling is not a big deal,
    // what we absolutely don't want is to reallocate for every profiled
    // section.
    profile_entries.reserve(old_capacity);

    current_level = 0;
}

std::size_t profile_begin(const char *label)
{
    profile_entries.emplace_back();
    profile_entries.back().label = label;
    profile_entries.back().level = current_level;
    ++current_level;
    profile_entries.back().t_start = std::chrono::steady_clock::now();
    return profile_entries.size() - 1;
}

void profile_end(std::size_t index)
{
    const auto t_end = std::chrono::steady_clock::now();
    profile_entries[index].t_end = t_end;
    --current_level;
}

const std::vector<Profile_entry> &profile_end_frame()
{
    assert((current_level == 0) && "Unmatched profile_begin/profile_end");
    return profile_entries;
}
