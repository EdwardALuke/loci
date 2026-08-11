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
#ifndef TRANSFER_PLAN_H
#define TRANSFER_PLAN_H

#include <vector>

/**
 * @file transfer_plan.h
 * @brief Converts refinement plans between general-face and quad-face child
 *        ordering conventions.
 *
 * FVMAdapt uses different child orderings for polygonal Face trees and
 * QuadFace trees. These helpers rebuild an empty tree from one representation,
 * traverse it breadth-first, and emit the corresponding plan for the other
 * representation.
 */


 /**
 * Converts a general-face refinement plan to a quad-face plan.
 *
 * The input is first replayed into a temporary Face tree. The output plan is
 * then written in breadth-first order using QuadFace child ordering. Trailing
 * no-split entries are removed before the result is returned.
 *
 * @param facePlan General-face plan to convert. The function does not modify
 *        the plan contents, but the parameter is non-const for compatibility
 *        with existing callers.
 * @return Equivalent plan using the quad-face convention.
 */
std::vector<char> transfer_plan_g2q(std::vector<char>& facePlan) ;


/**
 * Converts a quad-face refinement plan to a general-face plan.
 *
 * The input is replayed into a temporary QuadFace tree with orientation code
 * zero. The output plan is then written in breadth-first order using the
 * general Face child ordering. Trailing no-split entries are removed before the
 * result is returned.
 *
 * @param facePlan Quad-face plan to convert.
 * @return Equivalent plan using the general-face convention.
 */
std::vector<char> transfer_plan_q2g(const std::vector<char>& facePlan) ;

#endif
