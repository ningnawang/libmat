#include "io_utils.hpp"

#include <inttypes.h>

#include <fstream>

namespace IO {

void WriteHoudiniStr(std::ostream& ostream, const std::string& s) {
  write<BIGEND>(ostream, (short)s.size());
  ostream.write(s.c_str(), s.size());
}

GeometryWriter::GeometryWriter(const std::string& folder_path)
    : path_(folder_path) {}

#define IS_UNSIGNED(X) holds_alternative<vector<unsigned int>>(X)
#define IS_INTEGER(X) holds_alternative<vector<int>>(X)
#define IS_FLOAT(X) holds_alternative<vector<float>>(X)
#define IS_INT2(X) holds_alternative<vector<int2>>(X)
#define IS_INT3(X) holds_alternative<vector<int3>>(X)
#define IS_INT4(X) holds_alternative<vector<int4>>(X)
#define IS_FLOAT2(X) holds_alternative<vector<float2>>(X)
#define IS_FLOAT3(X) holds_alternative<vector<float3>>(X)
#define IS_FLOAT4(X) holds_alternative<vector<float4>>(X)

template void Geometry::AddParticleAttribute<unsigned int>(
    const std::string& name, const std::vector<unsigned int>& data);
template void Geometry::AddParticleAttribute<int>(const std::string& name,
                                                  const std::vector<int>& data);
template void Geometry::AddParticleAttribute<float>(
    const std::string& name, const std::vector<float>& data);
template void Geometry::AddParticleAttribute<int2>(
    const std::string& name, const std::vector<int2>& data);
template void Geometry::AddParticleAttribute<int3>(
    const std::string& name, const std::vector<int3>& data);
template void Geometry::AddParticleAttribute<int4>(
    const std::string& name, const std::vector<int4>& data);
template void Geometry::AddParticleAttribute<float2>(
    const std::string& name, const std::vector<float2>& data);
template void Geometry::AddParticleAttribute<float3>(
    const std::string& name, const std::vector<float3>& data);
template void Geometry::AddParticleAttribute<float4>(
    const std::string& name, const std::vector<float4>& data);

template void Geometry::AddPrimitiveAttribute<unsigned int>(
    const std::string& name, const std::vector<unsigned int>& data);
template void Geometry::AddPrimitiveAttribute<int>(
    const std::string& name, const std::vector<int>& data);
template void Geometry::AddPrimitiveAttribute<float>(
    const std::string& name, const std::vector<float>& data);
template void Geometry::AddPrimitiveAttribute<int2>(
    const std::string& name, const std::vector<int2>& data);
template void Geometry::AddPrimitiveAttribute<int3>(
    const std::string& name, const std::vector<int3>& data);
template void Geometry::AddPrimitiveAttribute<int4>(
    const std::string& name, const std::vector<int4>& data);
template void Geometry::AddPrimitiveAttribute<float2>(
    const std::string& name, const std::vector<float2>& data);
template void Geometry::AddPrimitiveAttribute<float3>(
    const std::string& name, const std::vector<float3>& data);
template void Geometry::AddPrimitiveAttribute<float4>(
    const std::string& name, const std::vector<float4>& data);

template void Geometry::AddPolygon<uint2>(const std::vector<uint2>& data);
template void Geometry::AddPolygon<uint3>(const std::vector<uint3>& data);
template void Geometry::AddPolygon<uint4>(const std::vector<uint4>& data);
template void Geometry::AddPolygon<std::vector<unsigned>>(
    const std::vector<std::vector<unsigned>>& data);

template <typename T>
void Geometry::AddParticleAttribute(const std::string& name,
                                    const std::vector<T>& data) {
  if (name == "P") {
    p_num_ = data.size();
    point_attr_data_.insert(point_attr_data_.begin(), data);

  } else {
    point_attr_data_.push_back(data);
    if (data.size() != p_num_ && p_num_ != 0) {
      throw std::runtime_error("data " + name +
                               " size does not match particle number");
    }
  }
  point_attr_name_.push_back(name);
}

template <typename T>
void Geometry::AddPrimitiveAttribute(const std::string& name,
                                     const std::vector<T>& data) {
  if (data.size() != poly_num_ && poly_num_ != 0) {
    throw std::runtime_error("data " + name +
                             " size does not match primitive size");
  }
  prim_attr_data_.push_back(data);
  prim_attr_name_.push_back(name);
}

template <typename T>
void Geometry::AddPolygon(const std::vector<T>& index) {
  poly_num_ += index.size();
  // poly_indices_.push_back(index);
  poly_indices_.insert(poly_indices_.end(), index.begin(), index.end());
}

void GeometryWriter::OutputGeometry(std::string name,
                                    const Geometry& data) const {
  using namespace std;
  if (data.p_num_ == 0) {
    throw runtime_error("Geometry data does not contain any particles");
  }

  const string filename = path_ + "/" + name + ".bgeo";
  const unique_ptr<ostream> fout(
      new ofstream(filename, ios::out | ios::binary));
  if (!*fout) {
    throw runtime_error("[IO::OutputGeometry] Failed to open file");
  } else {
    fout->imbue(locale::classic());
  }

  // Header part
  constexpr int magic = ((((('B' << 8) | 'g') << 8) | 'e') << 8) | 'o';
  constexpr char version_char = 'V';
  constexpr int version = 5;
  const int n_points = data.p_num_;
  const int n_prims = data.poly_num_;
  constexpr int n_point_groups = 0;
  constexpr int n_prim_groups = 0;
  const int n_point_attrib =
      data.point_attr_name_.size() - 1;  // do not count position
  constexpr int n_vertex_attrib = 0;
  const int n_prim_attrib = data.prim_attr_data_.size();  // primitive attribute
  constexpr int n_attrib = 0;
  write<BIGEND>(*fout, magic, version_char, version, n_points, n_prims,
                n_point_groups);
  write<BIGEND>(*fout, n_prim_groups, n_point_attrib, n_vertex_attrib,
                n_prim_attrib, n_attrib);

  int p_attrib_size = 4;  // at least contains position info
  vector<int> num_per_attrib_inc;
  num_per_attrib_inc.push_back(0);

  for (size_t i = 1; i < data.point_attr_name_.size(); ++i) {
    WriteHoudiniStr(*fout, data.point_attr_name_[i]);
    p_attrib_size = DefineAttribute(*fout, p_attrib_size, num_per_attrib_inc,
                                    data.point_attr_data_[i]);
  }

  // write particle data
  int* buffer = new int[p_attrib_size];
  for (int p_idx = 0; p_idx < n_points; p_idx++) {
    float3 pos = get<vector<float3>>(data.point_attr_data_[0])[p_idx];
    buffer[0] = reinterpret_cast<int&>(pos.x);
    buffer[1] = reinterpret_cast<int&>(pos.y);
    buffer[2] = reinterpret_cast<int&>(pos.z);
    float* w = (float*)&buffer[3];
    *w = 1.0f;
    for (size_t a_idx = 1; a_idx < data.point_attr_name_.size(); ++a_idx) {
      const int offset = num_per_attrib_inc[a_idx];
      WriteAttribute(buffer, p_idx, offset, data.point_attr_data_[a_idx]);
    }
    // big-endian
    for (int j = 0; j < p_attrib_size; j++) {
      BIGEND::endianSwap(buffer[j]);
    }
    fout->write((char*)buffer, p_attrib_size * sizeof(int));
  }
  delete[] buffer;

  // try to find maximum number of vertices inside a primitive
  int prim_size = 4;
  for (const auto& indices : data.poly_indices_) {
    if (holds_alternative<vector<unsigned>>(indices)) {
      const vector<unsigned>& prim = get<vector<unsigned>>(indices);
      prim_size = prim.size() > prim_size ? prim.size() : prim_size;
    }
  }

  // write primitive attributes
  int prim_attrib_size = 0;
  num_per_attrib_inc.clear();
  for (size_t i = 0; i < data.prim_attr_data_.size(); ++i) {
    WriteHoudiniStr(*fout, data.prim_attr_name_[i]);
    prim_attrib_size = DefineAttribute(
        *fout, prim_attrib_size, num_per_attrib_inc, data.prim_attr_data_[i]);
  }

  // write primitive data
  // NOTE: for detailed primitive definition, pls refer to HoudiniGPDlib
  // if the total number of the point is less than or equal to 65535,
  // 16 bit unsigned integers are used to store the point index.
  // otherwise 32 bit unsigned integers are used
  int* prim_attr_buffer = new int[prim_attrib_size];
  if (n_points > 1 << 16) {
    uint32_t* prim_buffer = new uint32_t[prim_size];
    for (int p_id = 0; p_id < data.poly_num_; ++p_id) {
      ReinterpretPolyBuffer<uint32_t>(prim_buffer, *fout,
                                      data.poly_indices_[p_id]);
      for (int a_id = 0; a_id < data.prim_attr_data_.size(); ++a_id) {
        WriteAttribute(prim_attr_buffer, p_id, num_per_attrib_inc[a_id],
                       data.prim_attr_data_[a_id]);
      }
      for (int i = 0; i < prim_attrib_size; ++i) {
        BIGEND::endianSwap(prim_attr_buffer[i]);
      }
      fout->write((char*)prim_attr_buffer, prim_attrib_size * sizeof(int));
    }
    delete[] prim_buffer;
  } else {
    uint16_t* prim_buffer = new uint16_t[prim_size];
    for (int p_id = 0; p_id < data.poly_num_; ++p_id) {
      ReinterpretPolyBuffer<uint16_t>(prim_buffer, *fout,
                                      data.poly_indices_[p_id]);
      for (int a_id = 0; a_id < data.prim_attr_data_.size(); ++a_id) {
        WriteAttribute(prim_attr_buffer, p_id, num_per_attrib_inc[a_id],
                       data.prim_attr_data_[a_id]);
      }
      for (int i = 0; i < prim_attrib_size; ++i) {
        BIGEND::endianSwap(prim_attr_buffer[i]);
      }
      fout->write((char*)prim_attr_buffer, prim_attrib_size * sizeof(int));
    }
    delete[] prim_buffer;
  }
  delete[] prim_attr_buffer;

  // Write extra
  write<BIGEND>(*fout, static_cast<char>(0x00));
  write<BIGEND>(*fout, static_cast<char>(0xff));
}

int GeometryWriter::DefineAttribute(std::ostream& output, int attrib_size,
                                    std::vector<int>& num_per_attrib_inc,
                                    const P_DATA_TYPE& attr_data) const {
  using namespace std;
  if (IS_UNSIGNED(attr_data)) {
    write<BIGEND>(output, static_cast<unsigned short>(1), AttributeType::INT,
                  0);
    num_per_attrib_inc.push_back(attrib_size);
    return attrib_size + 1;
  } else if (IS_INTEGER(attr_data)) {
    write<BIGEND>(output, static_cast<unsigned short>(1), AttributeType::INT,
                  0);
    num_per_attrib_inc.push_back(attrib_size);
    return attrib_size + 1;
  } else if (IS_FLOAT(attr_data)) {
    write<BIGEND>(output, static_cast<unsigned short>(1), AttributeType::FLOAT,
                  0.0f);
    num_per_attrib_inc.push_back(attrib_size);
    return attrib_size + 1;
  } else if (IS_INT2(attr_data)) {
    write<BIGEND>(output, static_cast<unsigned short>(2), AttributeType::INT, 0,
                  0);
    num_per_attrib_inc.push_back(attrib_size);
    return attrib_size + 2;
  } else if (IS_INT3(attr_data)) {
    write<BIGEND>(output, static_cast<unsigned short>(3), AttributeType::INT, 0,
                  0, 0);
    num_per_attrib_inc.push_back(attrib_size);
    return attrib_size + 3;
  } else if (IS_INT4(attr_data)) {
    write<BIGEND>(output, static_cast<unsigned short>(4), AttributeType::INT, 0,
                  0, 0, 0);
    num_per_attrib_inc.push_back(attrib_size);
    return attrib_size + 4;
  } else if (IS_FLOAT2(attr_data)) {
    write<BIGEND>(output, static_cast<unsigned short>(2), AttributeType::FLOAT,
                  0.0f, 0.0f);
    num_per_attrib_inc.push_back(attrib_size);
    return attrib_size + 2;
  } else if (IS_FLOAT3(attr_data)) {
    write<BIGEND>(output, static_cast<unsigned short>(3), AttributeType::FLOAT,
                  0.0f, 0.0f, 0.0f);
    num_per_attrib_inc.push_back(attrib_size);
    return attrib_size + 3;
  } else if (IS_FLOAT4(attr_data)) {
    write<BIGEND>(output, static_cast<unsigned short>(4), AttributeType::FLOAT,
                  0.0f, 0.0f, 0.0f, 0.0f);
    num_per_attrib_inc.push_back(attrib_size);
    return attrib_size + 4;
  } else {
    throw runtime_error("Unsupported attribute type!");
  }
}

void GeometryWriter::WriteAttribute(int* buffer, int idx, int offset,
                                    const P_DATA_TYPE& attr_data) const {
  using namespace std;
  if (IS_UNSIGNED(attr_data)) {
    unsigned int tmp = get<vector<unsigned int>>(attr_data)[idx];
    buffer[offset] = reinterpret_cast<int&>(tmp);
  } else if (IS_INTEGER(attr_data)) {
    int tmp = get<vector<int>>(attr_data)[idx];
    buffer[offset] = reinterpret_cast<int&>(tmp);
  } else if (IS_FLOAT(attr_data)) {
    float tmp = get<vector<float>>(attr_data)[idx];
    buffer[offset] = reinterpret_cast<int&>(tmp);
  } else if (IS_INT2(attr_data)) {
    int2 tmp = get<vector<int2>>(attr_data)[idx];
    buffer[offset] = reinterpret_cast<int&>(tmp.x);
    buffer[offset + 1] = reinterpret_cast<int&>(tmp.y);
  } else if (IS_INT3(attr_data)) {
    int3 tmp = get<vector<int3>>(attr_data)[idx];
    buffer[offset] = reinterpret_cast<int&>(tmp.x);
    buffer[offset + 1] = reinterpret_cast<int&>(tmp.y);
    buffer[offset + 2] = reinterpret_cast<int&>(tmp.z);
  } else if (IS_INT4(attr_data)) {
    int4 tmp = get<vector<int4>>(attr_data)[idx];
    buffer[offset] = reinterpret_cast<int&>(tmp.x);
    buffer[offset + 1] = reinterpret_cast<int&>(tmp.y);
    buffer[offset + 2] = reinterpret_cast<int&>(tmp.z);
    buffer[offset + 3] = reinterpret_cast<int&>(tmp.w);
  } else if (IS_FLOAT2(attr_data)) {
    float2 tmp = get<vector<float2>>(attr_data)[idx];
    buffer[offset] = reinterpret_cast<int&>(tmp.x);
    buffer[offset + 1] = reinterpret_cast<int&>(tmp.y);
  } else if (IS_FLOAT3(attr_data)) {
    float3 tmp = get<vector<float3>>(attr_data)[idx];
    buffer[offset] = reinterpret_cast<int&>(tmp.x);
    buffer[offset + 1] = reinterpret_cast<int&>(tmp.y);
    buffer[offset + 2] = reinterpret_cast<int&>(tmp.z);
  } else if (IS_FLOAT4(attr_data)) {
    float4 tmp = get<vector<float4>>(attr_data)[idx];
    buffer[offset] = reinterpret_cast<int&>(tmp.x);
    buffer[offset + 1] = reinterpret_cast<int&>(tmp.y);
    buffer[offset + 2] = reinterpret_cast<int&>(tmp.z);
    buffer[offset + 3] = reinterpret_cast<int&>(tmp.w);
  }
}

template <typename T>
void GeometryWriter::ReinterpretPolyBuffer(
    T* buffer, std::ostream& output, DEPRECATE_POLY_DATA_TYPE indices) const {
  using namespace std;
  if (holds_alternative<vector<uint2>>(indices)) {
    const vector<uint2>& line = get<vector<uint2>>(indices);
    for (auto idx : line) {
      write<BIGEND>(output, 0x00000001, 2, '<');
      buffer[0] = reinterpret_cast<T&>(idx.x);
      buffer[1] = reinterpret_cast<T&>(idx.y);
      for (int j = 0; j < 2; j++) {
        BIGEND::endianSwap(buffer[j]);
      }
      output.write((char*)buffer, 2 * sizeof(T));
    }
  } else if (holds_alternative<vector<uint3>>(indices)) {
    const vector<uint3>& tri = get<vector<uint3>>(indices);
    for (auto idx : tri) {
      write<BIGEND>(output, 0x00000001, 3, '<');
      buffer[0] = reinterpret_cast<T&>(idx.x);
      buffer[1] = reinterpret_cast<T&>(idx.y);
      buffer[2] = reinterpret_cast<T&>(idx.z);

      for (int j = 0; j < 3; j++) {
        BIGEND::endianSwap(buffer[j]);
      }
      output.write((char*)buffer, 3 * sizeof(T));
    }
  } else if (holds_alternative<vector<uint4>>(indices)) {
    const vector<uint4>& quad = get<vector<uint4>>(indices);
    for (auto idx : quad) {
      write<BIGEND>(output, 0x00000001, 4, '<');
      buffer[0] = reinterpret_cast<T&>(idx.x);
      buffer[1] = reinterpret_cast<T&>(idx.y);
      buffer[2] = reinterpret_cast<T&>(idx.z);
      buffer[3] = reinterpret_cast<T&>(idx.w);

      for (int j = 0; j < 4; j++) {
        BIGEND::endianSwap(buffer[j]);
      }
      output.write((char*)buffer, 4 * sizeof(T));
    }
  } else if (holds_alternative<vector<vector<unsigned>>>(indices)) {
    const vector<vector<unsigned>>& prim =
        get<vector<vector<unsigned>>>(indices);
    for (auto& p : prim) {
      const int num = p.size();
      write<BIGEND>(output, 0x00000001, num, '<');
      for (int i = 0; i < num; ++i) {
        unsigned idx = p[i];
        buffer[i] = reinterpret_cast<T&>(idx);
      }

      for (int j = 0; j < num; j++) {
        BIGEND::endianSwap(buffer[j]);
      }
      output.write((char*)buffer, num * sizeof(T));
    }
  } else {
    throw runtime_error("Unsupported data type for polygon indices");
  }
}

template <typename T>
void GeometryWriter::ReinterpretPolyBuffer(
    T* buffer, std::ostream& output, const POLY_DATA_TYPE& indices) const {
  using namespace std;
  if (holds_alternative<uint2>(indices)) {
    auto idx = get<uint2>(indices);
    write<BIGEND>(output, 0x00000001, 2, '<');
    buffer[0] = reinterpret_cast<T&>(idx.x);
    buffer[1] = reinterpret_cast<T&>(idx.y);
    for (int j = 0; j < 2; j++) {
      BIGEND::endianSwap(buffer[j]);
    }
    output.write((char*)buffer, 2 * sizeof(T));
  } else if (holds_alternative<uint3>(indices)) {
    auto idx = get<uint3>(indices);
    write<BIGEND>(output, 0x00000001, 3, '<');
    buffer[0] = reinterpret_cast<T&>(idx.x);
    buffer[1] = reinterpret_cast<T&>(idx.y);
    buffer[2] = reinterpret_cast<T&>(idx.z);

    for (int j = 0; j < 3; j++) {
      BIGEND::endianSwap(buffer[j]);
    }
    output.write((char*)buffer, 3 * sizeof(T));
  } else if (holds_alternative<uint4>(indices)) {
    auto idx = get<uint4>(indices);
    write<BIGEND>(output, 0x00000001, 4, '<');
    buffer[0] = reinterpret_cast<T&>(idx.x);
    buffer[1] = reinterpret_cast<T&>(idx.y);
    buffer[2] = reinterpret_cast<T&>(idx.z);
    buffer[3] = reinterpret_cast<T&>(idx.w);

    for (int j = 0; j < 4; j++) {
      BIGEND::endianSwap(buffer[j]);
    }
    output.write((char*)buffer, 4 * sizeof(T));
  } else if (holds_alternative<vector<unsigned>>(indices)) {
    const vector<unsigned> prim = get<vector<unsigned>>(indices);
    const int num = prim.size();
    write<BIGEND>(output, 0x00000001, num, '<');
    for (int i = 0; i < num; ++i) {
      unsigned idx = prim[i];
      buffer[i] = reinterpret_cast<T&>(idx);
    }

    for (int j = 0; j < num; j++) {
      BIGEND::endianSwap(buffer[j]);
    }
    output.write((char*)buffer, num * sizeof(T));
  } else {
    throw runtime_error("Unsupported data type for polygon indices");
  }
}
}  // namespace IO