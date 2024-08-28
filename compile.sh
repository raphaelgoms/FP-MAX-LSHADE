#! /bin/sh

g++ -std=c++17 main.cpp lshade.cpp _fp_max_lshade.cpp search_algorithm.cpp nl_shade_lbc.cpp \
 k_means_lshade.cpp cec14_test_func.cpp ./external/lib/*.cpp -I./external/include -lmlpack -fopenmp -o solver

#  nl_shade_rsp_mid.cpp