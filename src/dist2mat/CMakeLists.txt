cmake_minimum_required(VERSION 2.6)
add_library(dist2mat STATIC dist2mat_cuda_utils.h dist2mat_helper_math.h dist2mat.h dist2mat.cu)
