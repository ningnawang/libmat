#pragma once

#include <assert.h>
#include <cuda_runtime.h>
#include <nvtx3/nvToolsExt.h>

#include <stdexcept>
#include <vector>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#define kWarpSize 32

#define CUDA_CHECK(ans) \
  { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char* file, int line,
                      bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPU assert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    if (abort) exit(code);
  }
}

template <typename T>
class GpuBuffer {
 public:
  GpuBuffer() {
    h_.num_ = 0;
    h_.ptr_ = nullptr;
    d_.num_ = 0;
    d_.ptr_ = nullptr;
  }

  GpuBuffer(size_t num) {
    d_.num_ = num;
    h_.num_ = num;

    const size_t bytes = num * sizeof(T);
    // device
    CUDA_CHECK(cudaMalloc(&d_.ptr_, bytes));
    // host
    h_.ptr_ = new T[num];
  }

  ~GpuBuffer() { Clear(); }

  void Clear() {
    if (d_.ptr_ != nullptr) {
      CUDA_CHECK(cudaFree(d_.ptr_));
    }

    if (h_.ptr_ != nullptr) {
      delete[] h_.ptr_;
    }

    d_.ptr_ = nullptr;
    h_.ptr_ = nullptr;
    d_.num_ = 0;
    h_.num_ = 0;
  }

  cudaError_t DResize(size_t num, float scale = 1.f) {
    if (d_.num_ >= num) {
      return cudaSuccess;
    }

    // when we do the resize, allocate more memory than needed
    d_.num_ = num * scale;
    const size_t bytes = d_.num_ * sizeof(T);
    // free first
    if (d_.ptr_ != nullptr) {
      cudaError_t tmp = cudaFree(d_.ptr_);
      if (tmp != cudaSuccess) {
        printf("cuda free error!\n");
        return tmp;
      }
    }
    // then resize
    return cudaMalloc(reinterpret_cast<void**>(&(d_.ptr_)), bytes);
  }

  void HResize(size_t num) {
    if (h_.num_ >= num) {
      return;
    }

    h_.num_ = num;
    // const size_t bytes = num * sizeof(T);
    // free first
    if (h_.ptr_ != nullptr) {
      delete[] h_.ptr_;
    }
    // then resize
    h_.ptr_ = new T[num];
  }

  cudaError_t DFillZero(size_t num) {
    if (d_.num_ < num) {
      printf("Alert: not enough memory for filling zero, need to resize!\n");
      cudaError_t tmp = DResize(num);
      if (tmp != cudaSuccess) {
        printf("DResize in DFillZero fail!\n");
        return tmp;
      }
    }
    return cudaMemset((void*)d_.ptr_, 0, sizeof(T) * num);
  }

  void D2H(bool async = false) {
    // check whether host has enough memory to hold data
    if (h_.num_ < d_.num_) {
      HResize(d_.num_);
    }
    if (async) {  // need to call cuda device/event synchronize afterwards
      CUDA_CHECK(cudaMemcpyAsync(h_.ptr_, d_.ptr_, d_.num_ * sizeof(T),
                                 cudaMemcpyDeviceToHost));
    } else {
      CUDA_CHECK(cudaMemcpy(h_.ptr_, d_.ptr_, d_.num_ * sizeof(T),
                            cudaMemcpyDeviceToHost));
    }
  }

  void H2D() {
    // check whether device has enough memory to hold data
    if (d_.num_ < h_.num_) {
      DResize(h_.num_);
    }
    CUDA_CHECK(cudaMemcpy(d_.ptr_, h_.ptr_, d_.num_ * sizeof(T),
                          cudaMemcpyHostToDevice));
  }

  size_t DSize() const { return d_.num_; }
  size_t HSize() const { return h_.num_; }

  T* DPtr() { return d_.ptr_; }
  T* HPtr() { return h_.ptr_; }

  const T* DPtr() const { return d_.ptr_; }
  const T* HPtr() const { return h_.ptr_; }

 protected:
  struct PtrPair {
    size_t num_;
    T* ptr_;

    PtrPair() : num_(0), ptr_(nullptr) {}
  };

  PtrPair h_;
  PtrPair d_;
};