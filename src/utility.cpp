#include "utility.hpp"

#include <fstream>
#include <filesystem>

namespace
{

std::vector<Profile_entry> g_profile_entries{};
std::uint32_t g_current_level{};

} // namespace

Scope_profiler::Scope_profiler(const char *label)
{
    m_index = g_profile_entries.size();
    auto &entry = g_profile_entries.emplace_back();
    entry.label = label;
    entry.level = g_current_level;
    ++g_current_level;
    entry.t_start = std::chrono::steady_clock::now();
}

Scope_profiler::~Scope_profiler() noexcept
{
    stop();
}

void Scope_profiler::stop() const noexcept
{
    const auto t_end = std::chrono::steady_clock::now();
    if (g_profile_entries[m_index].t_end ==
        std::chrono::steady_clock::time_point{})
    {
        g_profile_entries[m_index].t_end = t_end;
        --g_current_level;
    }
}

void reset_profile_entries()
{
    const auto old_capacity = g_profile_entries.capacity();
    g_profile_entries.clear();
    // Clearing the vector should leave its capacity unchanged, but this is not
    // guaranteed. If the capacity did change, just reserve again. Having one
    // allocation per frame before starting the profiling is not a big deal,
    // what we absolutely don't want is to reallocate for every profiled
    // section.
    g_profile_entries.reserve(old_capacity);
    g_current_level = 0;
}

const std::vector<Profile_entry> &get_profile_entries()
{
    return g_profile_entries;
}

Eigen::ArrayXXf read_matrix_file(const char *file_name,
                                 int rows,
                                 int cols)
{
    if (!std::filesystem::exists(file_name))
    {
        std::ostringstream oss;
        oss << "File \"" << file_name << "\" does not exist\n";
        throw std::runtime_error(oss.str());
    }

    const auto expected_size = static_cast<std::uintmax_t>(rows * cols) * sizeof
                               (float);
    const auto file_size = std::filesystem::file_size(file_name);
    if (file_size != expected_size)
    {
        std::ostringstream oss;
        oss << "Wrong file size for \"" << file_name << "\": expected " << rows
            << "*" << cols << "*" << sizeof(float) << "=" << expected_size <<
            ", got " << file_size << "\n";
        throw std::runtime_error(oss.str());
    }

    std::ifstream file(file_name, std::ios::binary);
    if (!file.is_open())
    {
        std::ostringstream oss;
        oss << "Failed to open file \"" << file_name << "\"\n";
        throw std::runtime_error(oss.str());
    }

    Eigen::ArrayXXf array;
    array.resize(rows, cols);
    file.read(reinterpret_cast<char *>(array.data()),
              static_cast<std::streamsize>(file_size));

    return array;
}
