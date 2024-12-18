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

#pragma once
/**
 * \file mesh_AABB.h
 * \brief Axis Aligned Bounding Box trees for accelerating
 *  geometric queries that operate on a Mesh.
 */

#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>

namespace GEO {

/**
 * \brief Axis Aligned Bounding Box tree of mesh facets.
 * \details Used to quickly compute facet intersection and
 *  to locate the nearest facet from 3d query points.
 */
class GEOGRAM_API MeshFacetsAABBLambda : public MeshFacetsAABB {
 private:
  void get_nearest_facet_hint_lambdas(const vec3& p, index_t& nearest_f,
                                      double& lambda0, double& lambda1,
                                      double& lambda2, vec3& nearest_point,
                                      double& sq_dist) const;
  void nearest_facet_recursive_lambdas(const vec3& p, index_t& nearest_f,
                                       double& lambda0, double& lambda1,
                                       double& lambda2, vec3& nearest_point,
                                       double& sq_dist, index_t n, index_t b,
                                       index_t e) const;

 public:
  MeshFacetsAABBLambda(Mesh& M, bool reorder = true)
      : MeshFacetsAABB(M, reorder) {};
  ~MeshFacetsAABBLambda() {};

  /**
   * \brief Finds the nearest facet from an arbitrary 3d query point.
   * \param[in] p query point
   * \param[out] nearest_point nearest point on the surface
   * \param[out] sq_dist squared distance between p and the surface.
   * \return the index of the facet nearest to point p.
   */
  index_t nearest_facet_lambdas(const vec3& p, vec3& nearest_point,
                                double& sq_dist, double& lambda0,
                                double& lambda1, double& lambda2) const;
};

}  // namespace GEO