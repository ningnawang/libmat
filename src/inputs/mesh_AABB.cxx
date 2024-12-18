/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 */

#include "mesh_AABB.h"

#include <geogram/basic/geometry_nd.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/numerics/predicates.h>

namespace {

using namespace GEO;

/**
 * \brief Finds the nearest point in a mesh facet from a query point.
 * \param[in] M the mesh
 * \param[in] p the query point
 * \param[in] f index of the facet in \p M
 * \param[out] lambda1 barycentric coordinate of the closest point
 *  relative to facet \p f's vertex \p V0
 * \param[out] lambda2 barycentric coordinate of the closest point
 *  relative to facet \p f's vertex \p V1
 * \param[out] lambda3 barycentric coordinate of the closest point
 *  relative to facet \p f's vertex \p V2
 * \param[out] nearest_p the point of facet \p f nearest to \p p
 * \param[out] squared_dist the squared distance between
 *  \p p and \p nearest_p
 * \pre the mesh \p M is triangulated
 */
void get_point_facet_nearest_point_and_lambdas(const Mesh& M, const vec3& p,
                                               index_t f, double& lambda1,
                                               double& lambda2, double& lambda3,
                                               vec3& nearest_p,
                                               double& squared_dist) {
  geo_debug_assert(M.facets.nb_vertices(f) == 3);
  index_t c = M.facets.corners_begin(f);
  const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
  ++c;
  const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
  ++c;
  const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
  squared_dist = Geom::point_triangle_squared_distance(
      p, p1, p2, p3, nearest_p, lambda1, lambda2, lambda3);
}

/****************************************************************************/
/*No Change*/
/**
 * \brief Computes the axis-aligned bounding box of a mesh facet.
 * \param[in] M the mesh
 * \param[out] B the bounding box of the facet
 * \param[in] f the index of the facet in mesh \p M
 */
void get_facet_bbox(const Mesh& M, Box& B, index_t f) {
  index_t c = M.facets.corners_begin(f);
  const double* p = M.vertices.point_ptr(M.facet_corners.vertex(c));
  for (coord_index_t coord = 0; coord < 3; ++coord) {
    B.xyz_min[coord] = p[coord];
    B.xyz_max[coord] = p[coord];
  }
  for (++c; c < M.facets.corners_end(f); ++c) {
    p = M.vertices.point_ptr(M.facet_corners.vertex(c));
    for (coord_index_t coord = 0; coord < 3; ++coord) {
      B.xyz_min[coord] = std::min(B.xyz_min[coord], p[coord]);
      B.xyz_max[coord] = std::max(B.xyz_max[coord], p[coord]);
    }
  }
}

/**
 * \brief Computes the maximum node index in a subtree
 * \param[in] node_index node index of the root of the subtree
 * \param[in] b first facet index in the subtree
 * \param[in] e one position past the last facet index in the subtree
 * \return the maximum node index in the subtree rooted at \p node_index
 */
index_t max_node_index(index_t node_index, index_t b, index_t e) {
  geo_debug_assert(e > b);
  if (b + 1 == e) {
    return node_index;
  }
  index_t m = b + (e - b) / 2;
  index_t childl = 2 * node_index;
  index_t childr = 2 * node_index + 1;
  return std::max(max_node_index(childl, b, m), max_node_index(childr, m, e));
}

/**
 * \brief Computes the hierarchy of bounding boxes recursively.
 * \details This function is generic and can be used to compute
 *  a bbox hierarchy of arbitrary elements.
 * \param[in] M the mesh
 * \param[in] bboxes the array of bounding boxes
 * \param[in] node_index the index of the root of the subtree
 * \param[in] b first element index in the subtree
 * \param[in] e one position past the last element index in the subtree
 * \param[in] get_bbox a function that computes the bbox of an element
 * \tparam GET_BBOX a function (or a functor) with the following arguments:
 *  - mesh: a const reference to the mesh
 *  - box: a reference where the computed bounding box of the element
 *   will be stored
 *  - element: the index of the element
 */
template <class GET_BBOX>
void init_bboxes_recursive(const Mesh& M, vector<Box>& bboxes,
                           index_t node_index, index_t b, index_t e,
                           const GET_BBOX& get_bbox) {
  geo_debug_assert(node_index < bboxes.size());
  geo_debug_assert(b != e);
  if (b + 1 == e) {
    get_bbox(M, bboxes[node_index], b);
    return;
  }
  index_t m = b + (e - b) / 2;
  index_t childl = 2 * node_index;
  index_t childr = 2 * node_index + 1;
  geo_debug_assert(childl < bboxes.size());
  geo_debug_assert(childr < bboxes.size());
  init_bboxes_recursive(M, bboxes, childl, b, m, get_bbox);
  init_bboxes_recursive(M, bboxes, childr, m, e, get_bbox);
  geo_debug_assert(childl < bboxes.size());
  geo_debug_assert(childr < bboxes.size());
  bbox_union(bboxes[node_index], bboxes[childl], bboxes[childr]);
}

/**
 * \brief Finds the nearest point in a mesh facet from a query point.
 * \param[in] M the mesh
 * \param[in] p the query point
 * \param[in] f index of the facet in \p M
 * \param[out] nearest_p the point of facet \p f nearest to \p p
 * \param[out] squared_dist the squared distance between
 *  \p p and \p nearest_p
 * \pre the mesh \p M is triangulated
 */
void get_point_facet_nearest_point(const Mesh& M, const vec3& p, index_t f,
                                   vec3& nearest_p, double& squared_dist) {
  geo_debug_assert(M.facets.nb_vertices(f) == 3);
  index_t c = M.facets.corners_begin(f);
  const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
  ++c;
  const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
  ++c;
  const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
  double lambda1, lambda2, lambda3;  // barycentric coords, not used.
  squared_dist = Geom::point_triangle_squared_distance(
      p, p1, p2, p3, nearest_p, lambda1, lambda2, lambda3);
}

/**
 * \brief Computes the squared distance between a point and a Box.
 * \param[in] p the point
 * \param[in] B the box
 * \return the squared distance between \p p and \p B
 * \pre p is inside B
 */
double inner_point_box_squared_distance(const vec3& p, const Box& B) {
  geo_debug_assert(B.contains(p));
  double result = geo_sqr(p[0] - B.xyz_min[0]);
  result = std::min(result, geo_sqr(p[0] - B.xyz_max[0]));
  for (coord_index_t c = 1; c < 3; ++c) {
    result = std::min(result, geo_sqr(p[c] - B.xyz_min[c]));
    result = std::min(result, geo_sqr(p[c] - B.xyz_max[c]));
  }
  return result;
}

/**
 * \brief Computes the squared distance between a point and a Box
 *  with negative sign if the point is inside the Box.
 * \param[in] p the point
 * \param[in] B the box
 * \return the signed squared distance between \p p and \p B
 */
double point_box_signed_squared_distance(const vec3& p, const Box& B) {
  bool inside = true;
  double result = 0.0;
  for (coord_index_t c = 0; c < 3; c++) {
    if (p[c] < B.xyz_min[c]) {
      inside = false;
      result += geo_sqr(p[c] - B.xyz_min[c]);
    } else if (p[c] > B.xyz_max[c]) {
      inside = false;
      result += geo_sqr(p[c] - B.xyz_max[c]);
    }
  }
  if (inside) {
    result = -inner_point_box_squared_distance(p, B);
  }
  return result;
}

/**
 * \brief Computes the squared distance between a point and the
 *  center of a box.
 * \param[in] p the point
 * \param[in] B the box
 * \return the squared distance between \p p and the center of \p B
 */
double point_box_center_squared_distance(const vec3& p, const Box& B) {
  double result = 0.0;
  for (coord_index_t c = 0; c < 3; ++c) {
    double d = p[c] - 0.5 * (B.xyz_min[c] + B.xyz_max[c]);
    result += geo_sqr(d);
  }
  return result;
}

/**
 * \brief Tests whether a segment intersects a triangle.
 * \param[in] q1 , q2 the two extremities of the segment.
 * \param[in] p1 , p2 , p3 the three vertices of the triangle.
 * \retval true if [q1,q2] has an intersection with (p1, p2, p3).
 * \retval false otherwise.
 */
bool segment_triangle_intersection(const vec3& q1, const vec3& q2,
                                   const vec3& p1, const vec3& p2,
                                   const vec3& p3) {
  //   If the segment does not straddle the supporting plane of the
  // triangle, then there is no intersection.
  vec3 N = cross(p2 - p1, p3 - p1);
  if (dot(q1 - p1, N) * dot(q2 - p1, N) > 0.0) {
    return false;
  }

  //  The three tetrahedra formed by the segment and the three edges
  // of the triangle should have the same sign, else there is no
  // intersection.
  int s1 = geo_sgn(Geom::tetra_signed_volume(q1, q2, p1, p2));
  int s2 = geo_sgn(Geom::tetra_signed_volume(q1, q2, p2, p3));
  if (s1 != s2) {
    return false;
  }
  int s3 = geo_sgn(Geom::tetra_signed_volume(q1, q2, p3, p1));
  return (s2 == s3);
}

/**
 * \brief Tests whether there is an intersection between a segment
 *  and a mesh facet.
 * \param[in] q1 , q2 the extremities of the segment
 * \param[in] M the mesh
 * \param[in] f the facet
 */
bool segment_mesh_facet_intersection(const vec3& q1, const vec3& q2,
                                     const Mesh& M, index_t f) {
  geo_debug_assert(M.facets.nb_vertices(f) == 3);
  index_t c = M.facets.corners_begin(f);
  const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
  ++c;
  const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
  ++c;
  const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
  return segment_triangle_intersection(q1, q2, p1, p2, p3);
}

/**
 * \brief Tests whether a segment intersects a box.
 * \param[in] q1 , q2 the two extremities of the segment.
 * \param[in] box the box.
 * \retval true if [q1,q2] intersects the box.
 * \retval false otherwise.
 */
bool segment_box_intersection(const vec3& q1, const vec3& q2, const Box& box) {
  // Ref:
  // https://www.gamedev.net/forums/topic/338987-aabb---line-segment-intersection-test/
  vec3 d(0.5 * (q2.x - q1.x), 0.5 * (q2.y - q1.y), 0.5 * (q2.z - q1.z));

  vec3 e(0.5 * (box.xyz_max[0] - box.xyz_min[0]),
         0.5 * (box.xyz_max[1] - box.xyz_min[1]),
         0.5 * (box.xyz_max[2] - box.xyz_min[2]));

  vec3 c(q1.x + d.x - 0.5 * (box.xyz_min[0] + box.xyz_max[0]),
         q1.y + d.y - 0.5 * (box.xyz_min[1] + box.xyz_max[1]),
         q1.z + d.z - 0.5 * (box.xyz_min[2] + box.xyz_max[2]));

  vec3 ad(fabs(d.x), fabs(d.y), fabs(d.z));

  if (fabs(c[0]) > e[0] + ad[0]) {
    return false;
  }

  if (fabs(c[1]) > e[1] + ad[1]) {
    return false;
  }

  if (fabs(c[2]) > e[2] + ad[2]) {
    return false;
  }

  if (fabs(d[1] * c[2] - d[2] * c[1]) > e[1] * ad[2] + e[2] * ad[1]) {
    return false;
  }

  if (fabs(d[2] * c[0] - d[0] * c[2]) > e[2] * ad[0] + e[0] * ad[2]) {
    return false;
  }

  if (fabs(d[0] * c[1] - d[1] * c[0]) > e[0] * ad[1] + e[1] * ad[0]) {
    return false;
  }

  return true;
}
}  // namespace

/****************************************************************************/

namespace GEO {

void MeshFacetsAABBLambda::get_nearest_facet_hint_lambdas(
    const vec3& p, index_t& nearest_f, double& lambda0, double& lambda1,
    double& lambda2, vec3& nearest_point, double& sq_dist) const {
  // Find a good initial value for nearest_f by traversing
  // the boxes and selecting the child such that the center
  // of its bounding box is nearer to the query point.
  // For a large mesh (20M facets) this gains up to 10%
  // performance as compared to picking nearest_f randomly.
  index_t b = 0;
  index_t e = mesh_->facets.nb();
  index_t n = 1;
  while (e != b + 1) {
    index_t m = b + (e - b) / 2;
    index_t childl = 2 * n;
    index_t childr = 2 * n + 1;
    if (point_box_center_squared_distance(p, bboxes_[childl]) <
        point_box_center_squared_distance(p, bboxes_[childr])) {
      e = m;
      n = childl;
    } else {
      b = m;
      n = childr;
    }
  }
  nearest_f = b;

  index_t v =
      mesh_->facet_corners.vertex(mesh_->facets.corners_begin(nearest_f));
  lambda0 = 1.0;
  lambda1 = 0.0;
  lambda2 = 0.0;
  nearest_point = Geom::mesh_vertex(*mesh_, v);
  sq_dist = Geom::distance2(p, nearest_point);
}

void MeshFacetsAABBLambda::nearest_facet_recursive_lambdas(
    const vec3& p, index_t& nearest_f, double& lambda0, double& lambda1,
    double& lambda2, vec3& nearest_point, double& sq_dist, index_t n, index_t b,
    index_t e) const {
  geo_debug_assert(e > b);

  // If node is a leaf: compute point-facet distance
  // and replace current if nearer
  if (b + 1 == e) {
    vec3 cur_nearest_point;
    double cur_sq_dist;
    double cur_lambdas[3];
    get_point_facet_nearest_point_and_lambdas(*mesh_, p, b, cur_lambdas[0],
                                              cur_lambdas[1], cur_lambdas[2],
                                              cur_nearest_point, cur_sq_dist);
    if (cur_sq_dist < sq_dist) {
      nearest_f = b;
      nearest_point = cur_nearest_point;
      sq_dist = cur_sq_dist;
      lambda0 = cur_lambdas[0];
      lambda1 = cur_lambdas[1];
      lambda2 = cur_lambdas[2];
    }
    return;
  }
  index_t m = b + (e - b) / 2;
  index_t childl = 2 * n;
  index_t childr = 2 * n + 1;

  double dl = point_box_signed_squared_distance(p, bboxes_[childl]);
  double dr = point_box_signed_squared_distance(p, bboxes_[childr]);

  // Traverse the "nearest" child first, so that it has more chances
  // to prune the traversal of the other child.
  if (dl < dr) {
    if (dl < sq_dist) {
      nearest_facet_recursive_lambdas(p, nearest_f, lambda0, lambda1, lambda2,
                                      nearest_point, sq_dist, childl, b, m);
    }
    if (dr < sq_dist) {
      nearest_facet_recursive_lambdas(p, nearest_f, lambda0, lambda1, lambda2,
                                      nearest_point, sq_dist, childr, m, e);
    }
  } else {
    if (dr < sq_dist) {
      nearest_facet_recursive_lambdas(p, nearest_f, lambda0, lambda1, lambda2,
                                      nearest_point, sq_dist, childr, m, e);
    }
    if (dl < sq_dist) {
      nearest_facet_recursive_lambdas(p, nearest_f, lambda0, lambda1, lambda2,
                                      nearest_point, sq_dist, childl, b, m);
    }
  }
}

index_t MeshFacetsAABBLambda::nearest_facet_lambdas(
    const vec3& p, vec3& nearest_point, double& sq_dist, double& lambda0,
    double& lambda1, double& lambda2) const {
  lambda0 = lambda1 = lambda2 = -1;
  index_t nearest_facet;
  get_nearest_facet_hint_lambdas(p, nearest_facet, lambda0, lambda1, lambda2,
                                 nearest_point, sq_dist);
  nearest_facet_recursive_lambdas(p, nearest_facet, lambda0, lambda1, lambda2,
                                  nearest_point, sq_dist, 1, 0,
                                  mesh_->facets.nb());

  // if (lambda0 < -1e-10 || lambda1 < -1e-10 || lambda2 < -1e-10) {
  //   std::cout << "nearest_facet: " << nearest_facet << ", lambda0: " <<
  //   lambda0 << ", lambda1: " << lambda1 << ", lambda2: " << lambda2 <<
  //   std::endl;
  // }
  geo_assert(lambda0 >= -1e-10 && lambda1 >= -1e-10 && lambda2 >= -1e-10);

  // // check if point computed from lambdas is the same as nearest_point
  // index_t c = this->mesh_->facets.corners_begin(nearest_facet);
  // const vec3& p1 =
  //     Geom::mesh_vertex(*this->mesh_, this->mesh_->facet_corners.vertex(c));
  // ++c;
  // const vec3& p2 =
  //     Geom::mesh_vertex(*this->mesh_, this->mesh_->facet_corners.vertex(c));
  // ++c;
  // const vec3& p3 =
  //     Geom::mesh_vertex(*this->mesh_, this->mesh_->facet_corners.vertex(c));
  // vec3 tmp = p1 * lambda0 + p2 * lambda1 + p3 * lambda2;
  // geo_assert(Geom::distance(nearest_point, tmp) < 1e-2);

  return nearest_facet;
}
}  // namespace GEO