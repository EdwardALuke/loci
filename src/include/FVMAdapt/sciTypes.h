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
/*-*- C++ -*-*/
#ifndef LOCISCITYPES_H
#define LOCISCITYPES_H

#include <Loci_types.h>

const double EPSILON = 1e-30 ;

typedef double real;

typedef Loci::vector3d<unsigned int> vect3u ;
typedef Loci::vector3d<double> vect3d ;
typedef Loci::vector3d<vect3d> tens3d ;

// Inner product of two 3D vectors is the dot product
// dot(vect3d, vect3d) is already defined
// double dot(const vect3d lhs, const vect3d rhs);

/// Inner product of a 3D tensor of rank 2 and a 3D vector
vect3d inner(const tens3d & lhs, const vect3d & rhs);

/// Inner product of a 3D vector and a 3D tensor of rank 2
vect3d inner(const vect3d & v, const tens3d & t);

/// Outer product of two 3d vectors
tens3d outer(const vect3d & lhs, const vect3d & rhs);

/// This operation is similar to 3x3 Matrix-Matrix multiplication
tens3d product(const tens3d & t1, const tens3d & t2);

/// Compute the determinant of the 3x3 Matrix
double determinant(const tens3d & M);

/// Compute the transpose of the 3x3 Matrix
tens3d transpose(const tens3d & M);

/// Compute the inverse of the 3x3 Matrix
tens3d inverse(const tens3d & M);

#endif
