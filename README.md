# 582final
TOCOMPILE:

Note: you need to first download eigen3, libigl, glad, and glfw
Then, remember to change line 231 to your absolute path of vertices.txt on your local

clang++ -std=c++17 -I /Users/zhaoyanpeng/eigen3/ -I /Users/zhaoyanpeng/libigl/include/ -I /Users/zhaoyanpeng/glad/include -I /usr/local/include -o my_program main.cpp /Users/zhaoyanpeng/glad/src/glad.c -I /Users/zhaoyanpeng/glfw-3.3.8.bin.MACOS/include -L /opt/homebrew/lib -lglfw -framework OpenGL -lpthread -ldl

Main.cpp
this file generate an animation of twisting a rod. I didn't consider collision and gravity and fixed position of starting or end vertices. So the animation are just illustration aim.
Main_static.cpp
this file shows how bishop frame are parallel transporting from start point to the end. 
