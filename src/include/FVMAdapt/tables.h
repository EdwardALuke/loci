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
 *
 * Extract face split codes from cell plans and edge split codes from face
 * plans. Intermediate code 8 continues extraction at the next level without
 * splitting the current face or edge.
 */

/// Face split codes indexed by dd*7+cellCode-1, where dd is RIGHT, LEFT,
/// FRONT, BACK, UP, or DOWN and cellCode is in [1, 7]. Cell code 0 gives face
/// code 0 without a table lookup.
///
/// The two face directions are yz for RIGHT/LEFT, xz for FRONT/BACK, and xy
/// for UP/DOWN. Code 8 means the cell split leaves this face unsplit; continue
/// with the same face at the next extraction level.
const char faceCodeTable[42]={
  1, 2, 3, 8, 1, 2, 3,
  1, 2, 3, 8, 1, 2, 3,
  1, 8, 1, 2, 3, 2, 3,
  1, 8, 1, 2, 3, 2, 3,
  8, 1, 1, 2, 2, 3, 3,
  8, 1, 1, 2, 2, 3, 3} ;

const std::vector<bool> v1(2, true) ;
const std::vector<bool> v2(4, true) ;
const bool a3[3] = {0, 1, 0} ;
const bool a4[3] = {1, 0, 0} ;
const bool a5[5] = {0, 0, 1, 1, 0} ;
const bool a6[5] = {1, 1, 0, 0, 0} ;
const bool a7[5] = {1, 0, 1, 0, 0} ;
const bool a8[5] = {0, 1, 0, 1, 0} ;
const bool a9[9] = {0, 0, 0, 0, 1, 1, 1, 1, 0} ;
const bool a10[9] = {1, 1, 1, 1, 0, 0, 0, 0, 0} ;
const bool a11[9] = {0, 0, 1, 1, 0, 0, 1, 1, 0} ;
const bool a12[9] = {1, 1, 0, 0, 1, 1, 0, 0, 0} ;
const bool a13[9] = {0, 1, 0, 1, 0, 1, 0, 1, 0} ;
const bool a14[9] = {1, 0, 1, 0, 1, 0, 1, 0, 0} ;

const std::vector<bool> v3(a3, &a3[2]) ;
const std::vector<bool> v4(a4, &a4[2]) ;
const std::vector<bool> v5(a5, &a5[4]) ;
const std::vector<bool> v6(a6, &a6[4]) ;
const std::vector<bool> v7(a7, &a7[4]) ;
const std::vector<bool> v8(a8, &a8[4]) ;
const std::vector<bool> v9(a9, &a9[8]) ;
const std::vector<bool> v10(a10, &a10[8]) ;
const std::vector<bool> v11(a11, &a11[8]) ;
const std::vector<bool> v12(a12, &a12[8]) ;
const std::vector<bool> v13(a13, &a13[8]) ;
const std::vector<bool> v14(a14, &a14[8]) ;

/// Select the child cells that meet face dd, indexed by dd*7+cellCode-1 for
/// cellCode in [1, 7]. A true entry selects that child for face-plan
/// extraction. For example, {1, 1, 0, 0} selects children 0 and 1 from a
/// four-child cell.
const std::vector<bool> faceIDTable[42] = {
  v1, v1, v2, v3, v5, v5, v9,
  v1, v1, v2, v4, v6, v6, v10,
  v1, v3, v5, v1, v2, v8, v11,
  v1, v4, v6, v1, v2, v7, v12,
  v3, v1, v8, v1, v8, v2, v13,
  v4, v1, v7, v1, v7, v2, v14} ;

/* this is the IDtable
{
  {1, 1}, {1, 1}, {1, 1, 1, 1}, {0, 1}, {0, 0, 1, 1}, {0, 0, 1, 1}, {0, 0, 0, 0, 1, 1, 1, 1},
  {1, 1}, {1, 1}, {1, 1, 1, 1}, {1, 0}, {1, 1, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 1, 0, 0, 0, 0},
  {1, 1}, {0, 1}, {0, 0, 1, 1}, {1, 1}, {1, 1, 1, 1}, {0, 1, 0, 1}, {0, 0, 1, 1, 0, 0, 1, 1},
  {1, 1}, {1, 0}, {1, 1, 0, 0}, {1, 1}, {1, 1, 1, 1}, {1, 0, 1, 0}, {1, 1, 0, 0, 1, 1, 0, 0},
  {0, 1}, {1, 1}, {0, 1, 0, 1}, {1, 1}, {0, 1, 0, 1}, {1, 1, 1, 1}, {0, 1, 0, 1, 0, 1, 0, 1},
  {1, 0}, {1, 1}, {1, 0, 1, 0}, {1, 1}, {1, 0, 1, 0}, {1, 1, 1, 1}, {1, 0, 1, 0, 1, 0, 1, 0}};
*/

/// Edge split codes indexed by dd*3+faceCode-1, where dd is in [0, 4) and
/// faceCode is in [1, 3]. Face code 0 gives edge code 0 without a table
/// lookup.
///
/// The origin is node 0, x runs from node 0 to node 1, and y runs from node 0
/// to node 3. Edge directions are:
///
/// @verbatim
/// edge 0: node 0 -> node 1, y = 0
/// edge 1: node 1 -> node 2, x = 1
/// edge 2: node 3 -> node 2, y = 1
/// edge 3: node 0 -> node 3, x = 0
/// @endverbatim
///
/// Code 8 means the face split leaves this edge unsplit; continue with the
/// same edge at the next extraction level.
const char edgeCodeTable[12]={
  8, 1, 1,
  1, 8, 1,
  8, 1, 1,
  1, 8, 1} ;

/// Select the child faces that meet edge dd, indexed by dd*3+faceCode-1 for dd
/// in [0, 4) and faceCode in [1, 3]. A true entry selects that child for
/// edge-plan extraction.
const std::vector<bool> edgeIDTable[12] = {
  v4, v1, v7,
  v1, v3, v5,
  v3, v1, v8,
  v1, v4, v6} ;

/* this is the IDtable
   {
   {1, 0}, {1, 1}, {1, 0, 1, 0},
   {0, 1}, {1, 1}, {0, 1, 0, 1},
   {1, 1}, {1, 0}, {1, 1, 0, 0},
   {1, 1}, {0, 1}, {0, 0, 1, 1}};
*/
#endif

