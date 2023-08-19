#ifndef CONTEXT_HPP
#define CONTEXT_HPP

#include <concepts>
#include <stdexcept>

class Context
{
public:
    constexpr Context() = default;
    Context(int window_width, int window_height);
    ~Context() noexcept;

    Context(Context &&other) noexcept;
    Context &operator=(Context &&other) noexcept;

    Context(const Context &) = delete;
    Context &operator=(const Context &) = delete;

    void init(int window_width, int window_height);

    void run(std::invocable auto &&update) const
    {
        if (m_window == nullptr)
        {
            throw std::runtime_error(
                "Context needs to be initialized before calling run()");
        }

        while (!window_should_close())
        {
            begin_frame();
            update();
            end_frame();
        }
    }

private:
    [[nodiscard]] bool window_should_close() const;
    static void begin_frame();
    void end_frame() const;

    static constinit bool s_initialized;
    struct GLFWwindow *m_window {};
};

#endif // CONTEXT_HPP
