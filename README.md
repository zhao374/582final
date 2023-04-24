# 582final
TOCOMPILE:

Note: you need to first download eigen3, libigl, glad, and glfw

clang++ -std=c++17 -I /Users/zhaoyanpeng/eigen3/ -I /Users/zhaoyanpeng/libigl/include/ -I /Users/zhaoyanpeng/glad/include -I /usr/local/include -o my_program main.cpp /Users/zhaoyanpeng/glad/src/glad.c -I /Users/zhaoyanpeng/glfw-3.3.8.bin.MACOS/include -L /opt/homebrew/lib -lglfw -framework OpenGL -lpthread -ldl
