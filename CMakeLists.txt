cmake_minimum_required(VERSION 3.16)
project(Practical1)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_OPENGL		"Build target igl::opengl"		ON)
option(LIBIGL_GLFW		"Build target igl::glfw"		ON)
option(LIBIGL_IMGUI            "Build target igl::imgui"            ON)


include(libigl)

# Add your project files
file(GLOB SRCFILES *.cpp *.h *.c)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin PUBLIC igl::core igl::glfw igl::opengl igl::imgui)
