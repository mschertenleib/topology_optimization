#include "application.hpp"
#include "fea.hpp"
#include "utility.hpp"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "implot.h"

#include <GLFW/glfw3.h>

#include <cstdlib>
#include <iostream>

namespace
{

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
        "##heat_scale",
        scale_min,
        scale_max,
        ImVec2(100, 400));
    ImPlot::PopColormap();
}

void update_ui(
    FEA_state &fea,
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    &displacements_x,
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    &displacements_y,
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    &densities,
    const std::vector<Profile_entry> &profile_entries)
{
    displacements_x = fea.displacements(Eigen::seq(0, Eigen::last, 2))
                         .eval()
                         .reshaped(fea.num_nodes_y, fea.num_nodes_x);
    displacements_y = fea.displacements(Eigen::seq(1, Eigen::last, 2))
                         .eval()
                         .reshaped(fea.num_nodes_y, fea.num_nodes_x);
    densities = fea.design_variables_physical.reshaped(fea.num_elements_y,
        fea.num_elements_x);

    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

    if (ImGui::Begin("Stats"))
    {
        ImGui::Text("%.3f ms/frame (%.1f fps)",
                    1000.0 / static_cast<double>(ImGui::GetIO().Framerate),
                    static_cast<double>(ImGui::GetIO().Framerate));

        double parent_duration_ms{};
        for (std::size_t i{0}; i < profile_entries.size(); ++i)
        {
            const auto &entry = profile_entries[i];

            std::ostringstream indent;
            for (std::uint32_t level{0}; level < entry.level; ++level)
            {
                indent << "|   ";
            }

            const auto duration_ms =
                std::chrono::duration<double>(entry.t_end - entry.t_start)
                .count() *
                1000.0;

            double duration_percentage;
            if (entry.level == 0)
            {
                parent_duration_ms = duration_ms;
                duration_percentage = 100.0;
            }
            else
            {
                duration_percentage =
                    duration_ms / parent_duration_ms * 100.0;
            }

            ImGui::Text("%s%s: %.3f ms (%.2f%%)",
                        indent.str().c_str(),
                        entry.label,
                        duration_ms,
                        duration_percentage);
        }
    }
    ImGui::End();

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
        [](int error, const char *description)
        {
            std::cerr << "GLFW error " << error << ": " << description
                << std::endl;
        });

    if (!glfwInit())
    {
        return EXIT_FAILURE;
    }
    SCOPE_EXIT([] { glfwTerminate(); });

    constexpr auto glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    const auto window =
        glfwCreateWindow(1000, 850, "Topology Optimization", nullptr, nullptr);
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

    auto &io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

    float y_scale{1.0f};
    glfwGetWindowContentScale(window, nullptr, &y_scale);
    ImFontConfig font_config{};
    font_config.SizePixels = 13.0f * y_scale;
    io.Fonts->AddFontDefault(&font_config);

    if (!ImGui_ImplGlfw_InitForOpenGL(window, true))
    {
        return EXIT_FAILURE;
    }
    SCOPE_EXIT([] { ImGui_ImplGlfw_Shutdown(); });

    if (!ImGui_ImplOpenGL3_Init(glsl_version))
    {
        return EXIT_FAILURE;
    }
    SCOPE_EXIT([] { ImGui_ImplOpenGL3_Shutdown(); });

    ImPlot::CreateContext();
    SCOPE_EXIT([] { ImPlot::DestroyContext(); });

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        displacements_x(fea.num_nodes_y, fea.num_nodes_x);
    displacements_x.setZero();
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        displacements_y(fea.num_nodes_y, fea.num_nodes_x);
    displacements_y.setZero();
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        densities(fea.num_elements_y, fea.num_elements_x);
    densities.setZero();

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        reset_profile_entries();

        fea_optimization_step(fea);

        const auto &profile_entries = get_profile_entries();

        update_ui(
            fea,
            displacements_x,
            displacements_y,
            densities,
            profile_entries);

        ImGui::Render();
        int framebuffer_width{};
        int framebuffer_height{};
        glfwGetFramebufferSize(window, &framebuffer_width, &framebuffer_height);
        glViewport(0, 0, framebuffer_width, framebuffer_height);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    return EXIT_SUCCESS;
}
