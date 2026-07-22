//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#ifndef FVMADAPT_PLAN_OPERATIONS_H
#define FVMADAPT_PLAN_OPERATIONS_H

#include <Loci.h>
#include <vector>

#include "hex_defines.h"

void extract_quad_edge(const std::vector<char>& facePlan,
                       std::vector<char>& edgePlan,
                       unsigned int faceEdgeID) ;

void encode_edge_plan(const Loci::SetLong& splitCoordinates,
                      std::vector<char>& edgePlan) ;

/// Return one split-coordinate set per face-local edge.
/// Coordinates follow the stored edge2node orientation.
std::vector<Loci::SetLong> project_face_plan_to_edge_splits(
  const std::vector<char>& facePlan,
  bool isQuadFace,
  const std::vector<bool>& edgeIsReversed) ;

std::vector<char> merge_quad_face(std::vector<char>& facePlan,
                                  char orientCode) ;

std::vector<char> merge_quad_face(std::vector<char>& facePlanL,
                                  char orientCodeL,
                                  std::vector<char>& facePlanR,
                                  char orientCodeR) ;

std::vector<char> extract_hex_face(const std::vector<char>& cellPlan,
                                   DIRECTION direction) ;

std::vector<char> extract_prism_face(const std::vector<char>& cellPlan,
                                     int faceID) ;

std::vector<char> merge_quad_face_p(const std::vector<char>& cellPlan,
                                    int faceID, char orientCode) ;

std::vector<char> merge_quad_face_pp(const std::vector<char>& cellPlanL,
                                     int faceIDL, char orientCodeL,
                                     const std::vector<char>& cellPlanR,
                                     int faceIDR, char orientCodeR) ;

std::vector<char> merge_tri_face_p(const std::vector<char>& cellPlan,
                                   int faceID, char orientCode) ;

std::vector<char> merge_tri_face_pp(const std::vector<char>& cellPlanL,
                                    int faceIDL, char orientCodeL,
                                    const std::vector<char>& cellPlanR,
                                    int faceIDR, char orientCodeR) ;

std::vector<char> extract_general_face(
  const Entity* lower, int lower_size,
  const Entity* upper, int upper_size,
  const Entity* boundary_map, int boundary_map_size,
  const const_multiMap& face2node,
  const const_multiMap& face2edge,
  const const_MapVec<2>& edge2node,
  const std::vector<char>& cellPlan,
  Entity face,
  const const_store<int>& node_remap) ;

std::vector<char> merge_faceplan(std::vector<char>& facePlanL,
                                 std::vector<char>& facePlanR,
                                 int faceNodeCount) ;

std::vector<char> transfer_plan_g2q(std::vector<char>& facePlan) ;

std::vector<char> transfer_plan_q2g(const std::vector<char>& facePlan) ;

#endif
