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
#ifndef TABLES_H
#define TABLES_H
#include <vector>

/**
 * @file tables.h
 * @ingroup fvmadapt_plans
 * @brief Cell-to-face and face-to-edge split-code extraction tables.
 *
 * These tables are part of the low-level refinement-plan encoding. Code `8`
 * is a propagation case, not a synonym for `0`.
 */

/// faceCodeTable is the transfer table from cell code to face code
/// faceCode = faceCodeTable[dd*7+(int(cellCode)-int(1))]
/// Here dd is the face direction, it can be RIGHT, LEFT, FRONT, BACK
/// UP, DOWN. and cellCode is 1~7, 0 not included because if cellCode
/// is 0, faceCode is always 0.
/// 8 is a special value, means even cell is split, the face is not, but
/// the face can be split in the next step. So 8 can not be replaced by 0
/// when building a tree, 8 means the face has only one child.
/// The faceCode has the same meaning as cellCode. The binary order is xy,
///  yz, xz respectively instead of xyz in cellCode.
const char faceCodeTable[42]={
  1, 2, 3, 8, 1, 2, 3,
  1, 2, 3, 8, 1, 2, 3,
  1, 8, 1, 2, 3, 2, 3,
  1, 8, 1, 2, 3, 2, 3,
  8, 1, 1, 2, 2, 3, 3,
  8, 1, 1, 2, 2, 3, 3};

/// Build the ID table
const std::vector<bool> v1(2, true);
const std::vector<bool> v2(4, true);
const bool a3[3] = {0, 1, 0};
const bool a4[3] = {1, 0, 0};
const bool a5[5] = {0, 0, 1, 1, 0};
const bool a6[5] = {1, 1, 0, 0, 0};
const bool a7[5] = {1, 0, 1, 0, 0};
const bool a8[5] = {0, 1, 0, 1, 0};
const bool a9[9] = {0, 0, 0, 0, 1, 1, 1, 1, 0};
const bool a10[9] = {1, 1, 1, 1, 0, 0, 0, 0, 0};
const bool a11[9] = {0, 0, 1, 1, 0, 0, 1, 1, 0};
const bool a12[9] = {1, 1, 0, 0, 1, 1, 0, 0, 0};
const bool a13[9] = {0, 1, 0, 1, 0, 1, 0, 1, 0};
const bool a14[9] = {1, 0, 1, 0, 1, 0, 1, 0, 0};

const std::vector<bool> v3(a3, &a3[2]);
const std::vector<bool> v4(a4, &a4[2]);
const std::vector<bool> v5(a5, &a5[4]);
const std::vector<bool> v6(a6, &a6[4]);
const std::vector<bool> v7(a7, &a7[4]);
const std::vector<bool> v8(a8, &a8[4]);
const std::vector<bool> v9(a9, &a9[8]);
const std::vector<bool> v10(a10, &a10[8]);
const std::vector<bool> v11(a11, &a11[8]);
const std::vector<bool> v12(a12, &a12[8]);
const std::vector<bool> v13(a13, &a13[8]);
const std::vector<bool> v14(a14, &a14[8]);

/// faceIDTable shows the children cells whose codes need to be transferd to face code.
/// std::vector<bool> childrenID = faceIDTable[dd*7+(int(cellCode)-int(1))];
/// Here dd is the face direction, it can be RIGHT, LEFT, FRONT, BACK UP, DOWN.
// cellCode is 1~7, for example, if childrenID = {1, 1, 0, 0}, the cell has 4
// children, child 0 and child 1's codes need to be transfered to face code.
const std::vector<bool> faceIDTable[42] = {
  v1, v1, v2, v3, v5, v5, v9,
  v1, v1, v2, v4, v6, v6, v10,
  v1, v3, v5, v1, v2, v8, v11,
  v1, v4, v6, v1, v2, v7, v12,
  v3, v1, v8, v1, v8, v2, v13,
  v4, v1, v7, v1, v7, v2, v14};

/* this is the IDtable
{
  {1, 1}, {1, 1}, {1, 1, 1, 1}, {0, 1}, {0, 0, 1, 1}, {0, 0, 1, 1}, {0, 0, 0, 0, 1, 1, 1, 1},
  {1, 1}, {1, 1}, {1, 1, 1, 1}, {1, 0}, {1, 1, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 1, 0, 0, 0, 0},
  {1, 1}, {0, 1}, {0, 0, 1, 1}, {1, 1}, {1, 1, 1, 1}, {0, 1, 0, 1}, {0, 0, 1, 1, 0, 0, 1, 1},
  {1, 1}, {1, 0}, {1, 1, 0, 0}, {1, 1}, {1, 1, 1, 1}, {1, 0, 1, 0}, {1, 1, 0, 0, 1, 1, 0, 0},
  {0, 1}, {1, 1}, {0, 1, 0, 1}, {1, 1}, {0, 1, 0, 1}, {1, 1, 1, 1}, {0, 1, 0, 1, 0, 1, 0, 1},
  {1, 0}, {1, 1}, {1, 0, 1, 0}, {1, 1}, {1, 0, 1, 0}, {1, 1, 1, 1}, {1, 0, 1, 0, 1, 0, 1, 0}};
*/

/// The local coordinates system of set up as the following: the origin is at node 0
/// x positive: node 0 -> node 1,  y positive: node 0-> node 3
/// the definition of edges:
///   edge0(node0->node1,y=0),
///   edge1(node1->node2, x=1),
///   edge2(node3->node2, y=1),
///   edge3(node0->node3, x=0)
/// The edgeCodeTable is the transfer table from face code to edge code
/// edgeCode = edgeCodeTable[dd*3+faceCode-1]
/// Here, dd is the edgeID, values are 0~3 and faceCode is 1~3.
/// 0 is not included because if faceCode
/// is 0, edgeCode is always 0
/// 8 is a special value, means even face is split, the edge is not, but
/// the edge can be split in the next step.
/// When building a tree, 8 means the edge has only one child.
const char edgeCodeTable[12]={
  8, 1, 1,
  1, 8, 1,
  8, 1, 1,
  1, 8, 1};

/// Build the ID table.
/// edgeIDTable shows the children faces whose codes need to be transferd to
/// edge code.
/// std::vector<bool> childrenID = edgeIDTable[dd*7+faceCode -1]
/// Here, dd is the edge ID, for example, if childrenID = {1, 1, 0, 0}, the
/// face has 4 children, child 0 and child 1's codes need to be transfered to
/// edge code.
const std::vector<bool> edgeIDTable[12] = {
  v4, v1, v7,
  v1, v3, v5,
  v3, v1, v8,
  v1, v4, v6};

/* this is the IDtable
   {
   {1, 0}, {1, 1}, {1, 0, 1, 0},
   {0, 1}, {1, 1}, {0, 1, 0, 1},
   {1, 1}, {1, 0}, {1, 1, 0, 0},
   {1, 1}, {0, 1}, {0, 0, 1, 1}};
*/
#endif

