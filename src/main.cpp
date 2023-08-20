
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

void update()
{
    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

    ImGui::ShowDemoWindow();

    ImPlot::ShowDemoWindow();
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
