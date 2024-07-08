#pragma once

#include <cuda_runtime.h>

#include <nlohmann/json.hpp>
#include <set>
#include <tuple>
#include <unordered_map>
#include <variant>
#include <vector>

namespace IO {

enum class AttributeType : int { FLOAT = 0, INT = 1, VECTOR = 2 };

/*
 * Refer to Partio (Big-endian vs Little-endian)
 */
struct BIGEND {
  template <class T>
  static void endianSwap(T& v) {
    T temp = v;
    char* src = (char*)&temp;
    char* dest = (char*)&v;
    for (unsigned int i = 0; i < sizeof(T); i++) {
      dest[i] = src[sizeof(T) - i - 1];
    }
  }
};

template <class E, class T>
void read(std::istream& input, T& d) {
  input.read((char*)&d, sizeof(T));
  E::endianSwap(d);
}

template <class E, class T1, class T2>
void read(std::istream& input, T1& d1, T2& d2) {
  read<E>(input, d1);
  read<E>(input, d2);
}

template <class E, class T1, class T2, class T3>
void read(std::istream& input, T1& d1, T2& d2, T3& d3) {
  read<E>(input, d1);
  read<E>(input, d2, d3);
}

template <class E, class T1, class T2, class T3, class T4>
void read(std::istream& input, T1& d1, T2& d2, T3& d3, T4& d4) {
  read<E>(input, d1);
  read<E>(input, d2, d3, d4);
}

template <class E, class T1, class T2, class T3, class T4, class T5>
void read(std::istream& input, T1& d1, T2& d2, T3& d3, T4& d4, T5& d5) {
  read<E>(input, d1);
  read<E>(input, d2, d3, d4, d5);
}

template <class E, class T1, class T2, class T3, class T4, class T5, class T6>
void read(std::istream& input, T1& d1, T2& d2, T3& d3, T4& d4, T5& d5, T6& d6) {
  read<E>(input, d1);
  read<E>(input, d2, d3, d4, d5, d6);
}

template <class E, class T>
void write(std::ostream& output, const T& d) {
  T copy = d;
  E::endianSwap(copy);
  output.write((char*)&copy, sizeof(T));
}

template <class E, class T1, class T2>
void write(std::ostream& output, const T1& d1, const T2& d2) {
  write<E>(output, d1);
  write<E>(output, d2);
}

template <class E, class T1, class T2, class T3>
void write(std::ostream& output, const T1& d1, const T2& d2, const T3& d3) {
  write<E>(output, d1);
  write<E>(output, d2, d3);
}

template <class E, class T1, class T2, class T3, class T4>
void write(std::ostream& output, const T1& d1, const T2& d2, const T3& d3,
           const T4& d4) {
  write<E>(output, d1);
  write<E>(output, d2, d3, d4);
}

template <class E, class T1, class T2, class T3, class T4, class T5>
void write(std::ostream& output, const T1& d1, const T2& d2, const T3& d3,
           const T4& d4, const T5& d5) {
  write<E>(output, d1);
  write<E>(output, d2, d3, d4, d5);
}

template <class E, class T1, class T2, class T3, class T4, class T5, class T6>
void write(std::ostream& output, const T1& d1, const T2& d2, const T3& d3,
           const T4& d4, const T5& d5, const T6& d6) {
  write<E>(output, d1);
  write<E>(output, d2, d3, d4, d5, d6);
}

// For writing attribute infos
void WriteHoudiniStr(std::ostream& ostream, const std::string& s);

using P_DATA_TYPE =
    std::variant<std::vector<unsigned int>, std::vector<int>,
                 std::vector<float>, std::vector<int2>, std::vector<int3>,
                 std::vector<int4>, std::vector<float2>, std::vector<float3>,
                 std::vector<float4>>;

using DEPRECATE_POLY_DATA_TYPE =
    std::variant<std::vector<uint2>, std::vector<uint3>, std::vector<uint4>,
                 std::vector<std::vector<unsigned>>>;
using POLY_DATA_TYPE = std::variant<uint2, uint3, uint4, std::vector<unsigned>>;

/**
 * \brief Contains all geometry information
 * Must has point attribute named as "P"
 * Must make sure point/primitive attribute size matches p_num_/poly_num_
 */
class Geometry final {
 public:
  size_t p_num_ = 0;
  size_t poly_num_ = 0;

  Geometry() = default;
  ~Geometry() = default;

  Geometry(const Geometry&) = delete;
  Geometry& operator=(const Geometry&) = delete;
  Geometry(Geometry&&) = delete;
  Geometry& operator=(Geometry&&) = delete;

  // particle attributes type
  std::vector<std::string> point_attr_name_;
  std::vector<P_DATA_TYPE> point_attr_data_;

  std::vector<POLY_DATA_TYPE> poly_indices_;
  std::vector<std::string> prim_attr_name_;
  std::vector<P_DATA_TYPE> prim_attr_data_;

  template <typename T>
  void AddParticleAttribute(const std::string& name,
                            const std::vector<T>& data);

  template <typename T>
  void AddPrimitiveAttribute(const std::string& name,
                             const std::vector<T>& data);

  template <typename T>
  void AddPolygon(const std::vector<T>& index);
};

class GeometryWriter final {
 public:
  std::string path_;

  GeometryWriter() = default;
  explicit GeometryWriter(const std::string& folder_path);
  ~GeometryWriter() = default;
  void OutputGeometry(std::string name, const Geometry& data) const;

 private:
  // Define point/polygon attribute before actually writing
  int DefineAttribute(std::ostream& output, int attrib_size,
                      std::vector<int>& num_per_attrib_inc,
                      const P_DATA_TYPE& attr_data) const;
  // Write point/polygon attribute
  void WriteAttribute(int* buffer, int idx, int offset,
                      const P_DATA_TYPE& attr_data) const;
  // Deprecated
  template <typename T>
  void ReinterpretPolyBuffer(T* buffer, std::ostream& output,
                             DEPRECATE_POLY_DATA_TYPE indices) const;

  // Cast each polygon indices in correct format & output in big-endian style
  template <typename T>
  void ReinterpretPolyBuffer(T* buffer, std::ostream& output,
                             const POLY_DATA_TYPE& indices) const;
};

}  // namespace IO
