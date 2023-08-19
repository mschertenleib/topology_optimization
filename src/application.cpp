#include "application.hpp"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
#include "imgui.h"
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#include "implot.h"

Application::Application() : m_context(1280, 720)
{
    ImGui::StyleColorsDark();
}

void Application::run()
{
    m_context.run([this] { update(); });
}

void Application::update()
{
    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

    ImGui::ShowDemoWindow();

    ImPlot::ShowDemoWindow();
}
