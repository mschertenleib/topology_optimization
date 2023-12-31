cmake_minimum_required(VERSION 3.19)

project(topology_optimization LANGUAGES CXX)


set(CPM_DOWNLOAD_VERSION 0.38.2)

if (CPM_SOURCE_CACHE)
    set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
elseif (DEFINED ENV{CPM_SOURCE_CACHE})
    set(CPM_DOWNLOAD_LOCATION "$ENV{CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
else ()
    set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
endif ()

if (NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
    message(STATUS "Downloading CPM.cmake to ${CPM_DOWNLOAD_LOCATION}")
    file(DOWNLOAD
            https://github.com/TheLartians/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake
            ${CPM_DOWNLOAD_LOCATION}
    )
endif ()

include(${CPM_DOWNLOAD_LOCATION})


find_package(OpenGL REQUIRED)

CPMAddPackage(
        GITLAB_REPOSITORY libeigen/eigen
        GIT_TAG fee5d60b
        DOWNLOAD_ONLY YES
)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
CPMAddPackage(
        GITHUB_REPOSITORY glfw/glfw
        GIT_TAG 3.3.8
)

CPMAddPackage(
        GITHUB_REPOSITORY ocornut/imgui
        GIT_TAG v1.89.8-docking
        DOWNLOAD_ONLY YES
)

CPMAddPackage(
        GITHUB_REPOSITORY epezent/implot
        GIT_TAG v0.15
        DOWNLOAD_ONLY YES
)


set(is_gnu_or_clang "$<CXX_COMPILER_ID:GNU,Clang,AppleClang>")
add_compile_options(
        $<${is_gnu_or_clang}:-march=native>
        $<${is_gnu_or_clang}:-ffast-math>
)
add_compile_definitions(
        #EIGEN_NO_AUTOMATIC_RESIZING
        EIGEN_RUNTIME_NO_MALLOC
        EIGEN_INITIALIZE_MATRICES_BY_NAN
)


#@formatter:off
add_compile_definitions(IMGUI_USER_CONFIG="${CMAKE_SOURCE_DIR}/src/imgui_config.h")
#@formatter:on
add_library(imgui STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui.h
        ${imgui_SOURCE_DIR}/imgui_internal.h
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/imgui_demo.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui SYSTEM
        PRIVATE
        src
        PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
)
target_link_libraries(imgui PRIVATE OpenGL::GL glfw)


add_library(implot STATIC
        ${implot_SOURCE_DIR}/implot.cpp
        ${implot_SOURCE_DIR}/implot.h
        ${implot_SOURCE_DIR}/implot_internal.h
        ${implot_SOURCE_DIR}/implot_items.cpp
        ${implot_SOURCE_DIR}/implot_demo.cpp
)
target_include_directories(implot SYSTEM
        PUBLIC
        ${implot_SOURCE_DIR}
        PRIVATE
        ${imgui_SOURCE_DIR}
)
target_link_libraries(implot PRIVATE imgui)


add_executable(topology_optimization
        src/main.cpp
        src/application.cpp src/application.hpp
        src/fea.cpp src/fea.hpp
        src/utility.cpp src/utility.hpp
)
target_include_directories(topology_optimization SYSTEM PRIVATE ${eigen_SOURCE_DIR})
target_compile_features(topology_optimization PRIVATE cxx_std_20)
target_link_libraries(topology_optimization PRIVATE
        OpenGL::GL
        glfw
        imgui
        implot
)


set(clang_warnings
        -Wall
        -Wextra
        -Wshadow
        -Wnon-virtual-dtor
        -Wold-style-cast
        -Wcast-align
        -Wunused
        -Woverloaded-virtual
        -Wpedantic
        -Wconversion
        -Wsign-conversion
        -Wnull-dereference
        -Wdouble-promotion
)
set(gcc_warnings
        ${clang_warnings}
        -Wmisleading-indentation
        -Wduplicated-cond
        -Wduplicated-branches
        -Wlogical-op
)
if (CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    target_compile_options(topology_optimization PRIVATE ${clang_warnings})
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(topology_optimization PRIVATE ${gcc_warnings})
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(topology_optimization PRIVATE /W4)
else ()
    message(WARNING "No warnings set for compiler ${CMAKE_CXX_COMPILER_ID}")
endif ()


message(CHECK_START "Checking for IPO support")
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported_result OUTPUT ipo_supported_output)
if (ipo_supported_result)
    message(CHECK_PASS "supported")
    set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION_MIN_SIZE_REL TRUE)
else ()
    message(CHECK_FAIL "not supported")
endif ()
