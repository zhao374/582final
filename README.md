# 582final
TOCOMPILE:

Note: you need to first download eigen3, libigl, glad, and glfw
Then, remember to change line 231 to your absolute path of vertices.txt on your local

clang++ -std=c++17 -I /Users/zhaoyanpeng/eigen3/ -I /Users/zhaoyanpeng/libigl/include/ -I /Users/zhaoyanpeng/glad/include -I /usr/local/include -o my_program main.cpp /Users/zhaoyanpeng/glad/src/glad.c -I /Users/zhaoyanpeng/glfw-3.3.8.bin.MACOS/include -L /opt/homebrew/lib -lglfw -framework OpenGL -lpthread -ldl

Main.cpp
this file generate an animation of twisting a rod. I didn't consider collision and gravity and fixed position of starting or end vertices. So the animation are just illustration aim.
Main_static.cpp
this file shows how bishop frame are parallel transporting from start point to the end. 
Overview
The code reads a list of vertex positions from a file, calculates the edges and edge lengths between consecutive vertices, and simulates the motion of the rod by applying forces and updating the rod's configuration.

The following steps are executed in the simulation loop:

Compute forces acting on the rod
Update the rod configuration based on the forces
Compute bishop frames and material frames
Visualize the updated rod configuration
Dependencies
Eigen: A C++ template library for linear algebra
libigl: A C++ library for geometry processing research and development
Code Structure
The program consists of the following functions:

read_vertex_positions: Reads the vertex positions from a file
create_edges_and_lengths: Creates edges between consecutive vertices and calculates edge lengths
compute_bishop_frames: Computes the bishop frames for the rod
compute_material_frames: Computes the material frames for the rod
compute_forces: Computes bending and twisting forces acting on the rod
update_rod_configuration: Updates the rod configuration based on the forces
draw_arrow: Draws an arrow in the viewer for visualization
main: The main function that runs the simulation
