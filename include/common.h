#pragma once

#include <limits.h>

#include <algorithm>
#include <array>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

#include "assert.h"

#define RANDOM_01() ((double)std::rand() / RAND_MAX)  // random double in [0, 1]
#define RANDOM_INT(l, h) \
  (l + std::rand() % (h - l))  // random int in [l, h)
                               // (https://stackoverflow.com/a/7560146)
#define FOR(I, UPPERBND) for (int I = 0; I < int(UPPERBND); ++I)

typedef double Scalar;        // use for calculating Euler
typedef unsigned char uchar;  // local indices with special values
typedef unsigned int uint;
typedef std::array<int, 2> aint2;  //  unique id of powercell edge, medial edge
typedef std::array<int, 3> aint3;  //  medial face
typedef std::array<int, 4> aint4;
typedef std::array<int, 5> aint5;
typedef std::array<double, 3> adouble3;  // for gui

// === Math utilities
constexpr double PI = 3.14159265358979323846;
constexpr double HALF_PI = PI / 180.;

#define SCALAR_ZERO_6 1e-6  // float, 1e-8 if double
#define SCALAR_ZERO_5 1e-5
#define SCALAR_ZERO_4 1e-4
#define SCALAR_ZERO_3 1e-3
#define SCALAR_ZERO_2 1e-2
#define SCALAR_ZERO_1 1e-1
#define SCALAR_ZERO_0 0
#define SCALAR_1 1
#define SCALAR_10 10
#define SCALAR_20 20

//----------------------------------------------------------------------------
constexpr int DIM = 3;
constexpr int UNK_FACE = -1;
constexpr uchar END_OF_LIST = 255;
constexpr int UNK_INT = -1;
constexpr float UNK_FLOAT = -1.f;
constexpr uchar UNK_UCHAR = UCHAR_MAX;
// every tet face is shared by 2 tets
constexpr float F_TET_ADJ_DEFAULT = 2.f;
// every new clip of convex cell is only shared by 1
// bcs we only care about the Euler of a single cell
constexpr float F_CELL_ADJ_DEFAULT = 1.f;
constexpr float UNK_FACE_ID = -1.f;
constexpr double EPS_DEGREE_10 = 10;
constexpr double EPS_DEGREE_20 = 20;
constexpr double EPS_DEGREE_30 = 30;
constexpr double EPS_DEGREE_ = 30;
constexpr double EPS_DEGREE_90 = 90;