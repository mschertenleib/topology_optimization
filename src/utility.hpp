#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/SparseExtra>

#include <chrono>
#include <stdexcept>
#include <utility>

template <typename F>
class Scope_exit
{
public:
    explicit Scope_exit(F &&f) : m_f(std::move(f))
    {
    }

    explicit Scope_exit(F &f) : m_f(f)
    {
    }

    ~Scope_exit() noexcept
    {
        m_f();
    }

    Scope_exit(const Scope_exit &) = delete;
    Scope_exit(Scope_exit &&) = delete;
    Scope_exit &operator=(const Scope_exit &) = delete;
    Scope_exit &operator=(Scope_exit &&) = delete;

private:
    F m_f;
};

template <typename F>
class Scope_fail
{
public:
    explicit Scope_fail(F &&f)
        : m_exception_count(std::uncaught_exceptions()), m_f(std::move(f))
    {
    }

    explicit Scope_fail(F &f)
        : m_exception_count(std::uncaught_exceptions()), m_f(f)
    {
    }

    ~Scope_fail() noexcept
    {
        if (std::uncaught_exceptions() > m_exception_count)
        {
            m_f();
        }
    }

    Scope_fail(const Scope_fail &) = delete;
    Scope_fail(Scope_fail &&) = delete;
    Scope_fail &operator=(const Scope_fail &) = delete;
    Scope_fail &operator=(Scope_fail &&) = delete;

private:
    int m_exception_count;
    F m_f;
};

#define CONCATENATE_IMPL(s1, s2) s1##s2
#define CONCATENATE(s1, s2)      CONCATENATE_IMPL(s1, s2)

#define SCOPE_EXIT(f) const Scope_exit CONCATENATE(scope_exit_, __LINE__)(f)
#define SCOPE_FAIL(f) const Scope_fail CONCATENATE(scope_fail_, __LINE__)(f)

struct Profile_entry
{
    const char *label {};
    std::uint32_t level {};
    std::chrono::steady_clock::time_point t_start {};
    std::chrono::steady_clock::time_point t_end {};
};

class Scope_profiler
{
public:
    explicit Scope_profiler(const char *label);
    ~Scope_profiler() noexcept;
    void stop() const noexcept;

    Scope_profiler(const Scope_profiler &) = delete;
    Scope_profiler(Scope_profiler &&) = delete;
    Scope_profiler &operator=(const Scope_profiler &) = delete;
    Scope_profiler &operator=(Scope_profiler &&) = delete;

private:
    std::size_t m_index;
};

void reset_profile_entries();
[[nodiscard]] const std::vector<Profile_entry> &get_profile_entries();

template <typename Matrix>
void load_market(Matrix &mat, const char *file_name)
{
    int sym;
    bool is_complex;
    bool is_dense;
    if (!Eigen::getMarketHeader(file_name, sym, is_complex, is_dense))
    {
        std::ostringstream oss;
        oss << "Failed to open file \"" << file_name << '\"';
        throw std::runtime_error(oss.str());
    }

    if (sym != 0 || is_complex)
    {
        std::ostringstream oss;
        oss << "Unsupported matrix type for file \"" << file_name << '\"';
        throw std::runtime_error(oss.str());
    }

    constexpr auto is_sparse_type = requires { Matrix::isCompressed; };
    if (is_sparse_type == is_dense)
    {
        std::ostringstream oss;
        if constexpr (is_sparse_type)
        {
            oss << "Reading \"" << file_name << "\" as sparse, but it is dense";
        }
        else
        {
            oss << "Reading \"" << file_name << "\" as dense, but it is sparse";
        }
        throw std::runtime_error(oss.str());
    }

    bool result;
    if constexpr (is_sparse_type)
    {
        result = Eigen::loadMarket(mat, file_name);
    }
    else
    {
        result = Eigen::loadMarketDense(mat, file_name);
    }
    if (!result)
    {
        std::ostringstream oss;
        oss << "Failed to parse file \"" << file_name << '\"';
        throw std::runtime_error(oss.str());
    }
}

template <typename T>
[[nodiscard]] Eigen::SparseMatrix<T> load_market_sparse(const char *file_name)
{
    Eigen::SparseMatrix<T> mat;
    load_market(mat, file_name);
    return mat;
}

template <typename T>
[[nodiscard]] Eigen::MatrixX<T> load_market_dense(const char *file_name)
{
    Eigen::MatrixX<T> mat;
    load_market(mat, file_name);
    return mat;
}

#endif // UTILITY_HPP
