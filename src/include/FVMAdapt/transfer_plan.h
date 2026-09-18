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
 *
 * Convert plans between four-edge Face trees and QuadFace trees. Conversion
 * to Face retains only four-way splits.
 */


/**
 * Convert a four-edge Face plan to QuadFace order. Map split code 1 to code
 * 3, reorder the children, and remove trailing zero entries. The contents of
 * facePlan are not modified.
 */
std::vector<char> transfer_plan_g2q(std::vector<char>& facePlan) ;


/**
 * Convert the four-way splits of a QuadFace plan to four-edge Face order. Map
 * code 3 to code 1 and reorder the children. Codes 1 and 2 become unsplit
 * entries, and their descendants are omitted. Remove trailing zero entries
 * before returning the plan.
 */
std::vector<char> transfer_plan_q2g(const std::vector<char>& facePlan) ;

#endif
