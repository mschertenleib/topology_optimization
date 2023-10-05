#include "application.hpp"
#include "fea.hpp"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "implot.h"

#include <GLFW/glfw3.h>

#include <cstdlib>
#include <iostream>
#include <semaphore>
#include <thread>
#include <utility>

namespace
{

std::binary_semaphore g_fea_computing_finished {0};
std::binary_semaphore g_fea_read_finished {0};

template <typename F>
class Scope_guard
{
public:
    explicit Scope_guard(F &&f) : m_f {std::forward<F>(f)}
    {
    }

    ~Scope_guard() noexcept
    {
        m_f();
    }

    Scope_guard(const Scope_guard &) = delete;
    Scope_guard(Scope_guard &&) noexcept = delete;
    Scope_guard &operator=(const Scope_guard &) = delete;
    Scope_guard &operator=(Scope_guard &&) noexcept = delete;

private:
    F m_f;
};

void plot_matrix(
    const char *title_id,
    const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        &matrix)
{
    const auto scale_min = static_cast<double>(matrix.minCoeff());
    const auto scale_max = static_cast<double>(matrix.maxCoeff());
    ImPlot::PushColormap(ImPlotColormap_Jet);
    if (ImPlot::BeginPlot(title_id, ImVec2(800, 400)))
    {
        ImPlot::SetupAxes(nullptr,
                          nullptr,
                          ImPlotAxisFlags_NoDecorations,
                          ImPlotAxisFlags_NoDecorations);
        ImPlot::SetupAxesLimits(0, 1, 0, 1, ImPlotCond_Always);
        ImPlot::PlotHeatmap("##heat_map",
                            matrix.data(),
                            static_cast<int>(matrix.rows()),
                            static_cast<int>(matrix.cols()),
                            scale_min,
                            scale_max,
                            nullptr);
        ImPlot::EndPlot();
    }
    ImGui::SameLine();
    ImPlot::ColormapScale(
        "##heat_scale", scale_min, scale_max, ImVec2(100, 400));
    ImPlot::PopColormap();
}

void update_ui(
    FEA_state &fea,
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        &displacements_x,
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        &displacements_y,
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        &densities)
{
    if (g_fea_computing_finished.try_acquire())
    {
        displacements_x = fea.displacements(Eigen::seq(0, Eigen::last, 2))
                              .eval()
                              .reshaped(fea.num_nodes_y, fea.num_nodes_x);
        displacements_y = fea.displacements(Eigen::seq(1, Eigen::last, 2))
                              .eval()
                              .reshaped(fea.num_nodes_y, fea.num_nodes_x);
        densities = fea.design_variables_physical.reshaped(fea.num_elements_y,
                                                           fea.num_elements_x);
        g_fea_read_finished.release();
    }

    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

    if (ImGui::Begin("Displacements"))
    {
        plot_matrix("Displacements X", displacements_x);
        plot_matrix("Displacements Y", displacements_y);
    }
    ImGui::End();
    if (ImGui::Begin("Densities"))
    {
        plot_matrix("Densities", densities);
    }
    ImGui::End();
}

} // namespace

int application_main(FEA_state &fea)
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
    const Scope_guard glfw_guard([] { glfwTerminate(); });

    constexpr auto glsl_version {"#version 150"};
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    const auto window {
        glfwCreateWindow(1000, 850, "mechanisms", nullptr, nullptr)};
    if (window == nullptr)
    {
        return EXIT_FAILURE;
    }
    const Scope_guard window_guard([window] { glfwDestroyWindow(window); });

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    const Scope_guard imgui_guard([] { ImGui::DestroyContext(); });

    auto &io {ImGui::GetIO()};
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

    float y_scale {1.0f};
    glfwGetWindowContentScale(window, nullptr, &y_scale);
    ImFontConfig font_config {};
    font_config.SizePixels = 13.0f * y_scale;
    io.Fonts->AddFontDefault(&font_config);

    if (!ImGui_ImplGlfw_InitForOpenGL(window, true))
    {
        return EXIT_FAILURE;
    }
    const Scope_guard imgui_platform_guard([] { ImGui_ImplGlfw_Shutdown(); });

    if (!ImGui_ImplOpenGL3_Init(glsl_version))
    {
        return EXIT_FAILURE;
    }
    const Scope_guard imgui_renderer_guard([]
                                           { ImGui_ImplOpenGL3_Shutdown(); });

    ImPlot::CreateContext();
    const Scope_guard implot_guard([] { ImPlot::DestroyContext(); });

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        displacements_x(fea.num_nodes_y, fea.num_nodes_x);
    displacements_x.setZero();
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        displacements_y(fea.num_nodes_y, fea.num_nodes_x);
    displacements_y.setZero();
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        densities(fea.num_elements_y, fea.num_elements_x);
    densities.setZero();

    std::jthread fea_thread(
        [&fea](const std::stop_token &stop_token)
        {
            while (!stop_token.stop_requested())
            {
                g_fea_read_finished.acquire();

                fea_optimization_step(fea);

                g_fea_computing_finished.release();
            }
        });

    g_fea_read_finished.release();

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        update_ui(fea, displacements_x, displacements_y, densities);

        ImGui::Render();
        int framebuffer_width {};
        int framebuffer_height {};
        glfwGetFramebufferSize(window, &framebuffer_width, &framebuffer_height);
        glViewport(0, 0, framebuffer_width, framebuffer_height);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    return EXIT_SUCCESS;
}
