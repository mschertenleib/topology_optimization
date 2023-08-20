
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
#include "imgui.h"
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "implot.h"

#include <GLFW/glfw3.h>

#include <Eigen/Dense>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace
{

template <typename F>
struct Scope_exit
{
    explicit Scope_exit(F &&f) : fn {std::forward<F>(f)}
    {
    }
    ~Scope_exit() noexcept
    {
        fn();
    }
    F fn;
};

template <typename F>
Scope_exit(F &&) -> Scope_exit<F>;

#define CONCAT_IMPL(s1, s2) s1##s2
#define CONCAT(s1, s2)      CONCAT_IMPL(s1, s2)

#define SCOPE_EXIT(f) const Scope_exit CONCAT(scope_exit_, __LINE__)(f)

void plot_matrix(const char *title_id, const Eigen::MatrixXf &matrix)
{
    const auto rows {static_cast<int>(matrix.rows())};
    const auto cols {static_cast<int>(matrix.cols())};

    const auto region_size {ImGui::GetContentRegionAvail()};
    if (region_size.x <= 0.0f || region_size.y <= 0.0f)
    {
        return;
    }
    const auto matrix_aspect_ratio {static_cast<float>(cols) /
                                    static_cast<float>(rows)};
    const auto region_aspect_ratio {region_size.x / region_size.y};
    const auto cursor_pos {ImGui::GetCursorPos()};
    auto width {region_size.x};
    auto height {region_size.y};
    auto x {cursor_pos.x};
    auto y {cursor_pos.y};
    if (matrix_aspect_ratio >= region_aspect_ratio)
    {
        height = width / matrix_aspect_ratio;
        y += (region_size.y - height) * 0.5f;
    }
    else
    {
        width = height * matrix_aspect_ratio;
        x += (region_size.x - width) * 0.5f;
    }
    ImGui::SetCursorPos({x, y});

    ImPlot::PushColormap(ImPlotColormap_Jet);
    if (ImPlot::BeginPlot(title_id, {width, height}))
    {
        ImPlot::SetupAxes(nullptr,
                          nullptr,
                          ImPlotAxisFlags_NoDecorations,
                          ImPlotAxisFlags_NoDecorations);
        ImPlot::SetupAxesLimits(0, 1, 0, 1, ImPlotCond_Always);
        ImPlot::PlotHeatmap(
            "##heatmap", matrix.data(), rows, cols, 0, 0, nullptr);
        ImPlot::EndPlot();
    }
    ImPlot::PopColormap();
}

void update()
{
    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

    ImPlot::ShowDemoWindow();

    constexpr auto rows {5};
    constexpr auto cols {10};
    Eigen::MatrixXf matrix {Eigen::VectorXf::LinSpaced(rows * cols, 0.0f, 1.0f)
                                .reshaped(rows, cols)};

    if (ImGui::Begin("Test window"))
    {
        plot_matrix("##matrix_plot", matrix);
    }
    ImGui::End();
}

} // namespace

int main()
{
    try
    {
        glfwSetErrorCallback(
            [](int error, const char *description) {
                std::cerr << "GLFW error " << error << ": " << description
                          << std::endl;
            });

        if (!glfwInit())
        {
            return EXIT_FAILURE;
        }
        SCOPE_EXIT(glfwTerminate);

#if defined(__APPLE__)
        constexpr auto glsl_version {"#version 150"};
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#else
        constexpr auto glsl_version {"#version 130"};
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
#endif

        const auto window {
            glfwCreateWindow(1280, 720, "mechanisms", nullptr, nullptr)};
        if (window == nullptr)
        {
            return EXIT_FAILURE;
        }
        SCOPE_EXIT([window] { glfwDestroyWindow(window); });

        glfwMakeContextCurrent(window);
        glfwSwapInterval(1);

        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        SCOPE_EXIT([] { ImGui::DestroyContext(); });

        auto &io {ImGui::GetIO()};
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;
        io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

        if (!ImGui_ImplGlfw_InitForOpenGL(window, true))
        {
            return EXIT_FAILURE;
        }
        SCOPE_EXIT(ImGui_ImplGlfw_Shutdown);

        if (!ImGui_ImplOpenGL3_Init(glsl_version))
        {
            return EXIT_FAILURE;
        }
        SCOPE_EXIT(ImGui_ImplOpenGL3_Shutdown);

        ImPlot::CreateContext();
        SCOPE_EXIT([] { ImPlot::DestroyContext(); });

        while (!glfwWindowShouldClose(window))
        {
            glfwPollEvents();

            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            update();

            ImGui::Render();
            int framebuffer_width {};
            int framebuffer_height {};
            glfwGetFramebufferSize(
                window, &framebuffer_width, &framebuffer_height);
            glViewport(0, 0, framebuffer_width, framebuffer_height);
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

            glfwSwapBuffers(window);
        }

        return EXIT_SUCCESS;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception thrown: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cerr << "Unknown exception thrown\n";
        return EXIT_FAILURE;
    }
}
